# Author: Anna Lyubetskaya. Date: 21-09-03
# Perform DEA of spots with specific pathology categories
# In this setting, the more data the better; integration (common clusters) or even SCT normalization are not necessary


## ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_signatures.R")

source("code/R/b_Cluster_Biomarkers/Seurat_cluster_utils.R")


## PARAMETERS ----


# Path to processed Seurat data
sample_name <- "Syng_merge"

# User-defined score threshold
group_threshold <- 100

#PCT threshold - percent of spots containing top biomarkers
pct1_threshold <- 90


## MC38 mixed v pure

# Name for the pathology class to define groups
#pathology_select_gr1 <- "Pathology.MC38.percent"
#pathology_select_gr2 <- "Pathology.MC38.percent"

# Sample name regex to apply to each group
#sample_list_gr1 <- "_B16-MC38_"
#sample_list_gr2 <- "_MC38_"

# Name to use in folder and file names
#comparison_name <- "MC38m_v_MC38p"
#gr1_name <- "MC38m"
#gr2_name <- "MC38p"

# User defined color scheme
#cols <- list("MC38p" = "#056db5", "MC38m" = "#153c65", "Exclude" = "lightgrey")
# Add user defined genes to the final dot plot if they are DEA
#genes_user <- c()


## B16 mixed v pure

# Name for the pathology class to define groups
pathology_select_gr1 <- "Pathology.B16.percent"
pathology_select_gr2 <- "Pathology.B16.percent"

# Sample name regex to apply to each group
sample_list_gr1 <- "_B16-MC38_"
sample_list_gr2 <- "_B16_"

# Name to use in folder and file names
comparison_name <- "B16m_v_B16p"
gr1_name <- "B16m"
gr2_name <- "B16p"

# User defined color scheme
cols <- list("B16p" = "#056db5", "B16m" = "#153c65", "Exclude" = "lightgrey")

# Add user defined genes to the final dot plot if they are DEA
genes_user <- c()


## B16 pure v MC38 pure

# Name for the pathology class to define groups
#pathology_select_gr1 <- "Pathology.B16.percent"
#pathology_select_gr2 <- "Pathology.MC38.percent"

# Sample name regex to apply to each group
#sample_list_gr1 <- "_B16_"
#sample_list_gr2 <- "_MC38_"

# Name to use in folder and file names
#comparison_name <- "B16p_v_MC38p"
#gr1_name <- "B16p"
#gr2_name <- "MC38p"

# User defined color scheme
#cols <- list("B16p" = "#056db5", "MC38p" = "#db5a7c", "Exclude" = "lightgrey")

# Add user defined genes to the final dot plot if they are DEA
#genes_user <- c("TYRP1", "DCT", "PMEL", "MLANA", "SHISAL2B", "RHOX5", "CD74", "CSF1R", "PTPRC", "COL3A1", "DCN", "COL6A1")


## B16 pure v MC38 pure

# Name for the pathology class to define groups
pathology_select_gr1 <- "Pathology.B16.percent"
pathology_select_gr2 <- "Pathology.MC38.percent"

# Sample name regex to apply to each group
sample_list_gr1 <- "_B16-MC38_"
sample_list_gr2 <- "_B16-MC38_"

# Name to use in folder and file names
comparison_name <- "B16m_v_MC38m"
gr1_name <- "B16m"
gr2_name <- "MC38m"

# User defined color scheme
cols <- list("B16m" = "#056db5", "MC38m" = "#db5a7c", "Exclude" = "lightgrey")

# Add user defined genes to the final dot plot if they are DEA
#genes_user <- NULL


## HumanPanc: DonorA FF-polyA v FFPE-polyA

# Name for the pathology class to define groups
#pathology_select_gr1 <- "Pathology.Tissue"
#pathology_select_gr2 <- "Pathology.Tissue"

# Sample name regex to apply to each group
#sample_list_gr1 <- "_FF_"
#sample_list_gr2 <- "ROI1_FFPE_.*Jun20|ROI1_FFPE_.*Dec20"

# Name to use in folder and file names
#comparison_name <- "FF-polyA_v_FFPE-polyA"
#gr1_name <- "FF-polyA"
#gr2_name <- "FFPE-polyA"

# User defined color scheme
#cols <- list("FF-polyA" = "#056db5", "FFPE-polyA" = "#db5a7c", "Exclude" = "lightgrey")

# Add user defined genes to the final dot plot if they are DEA
#genes_user <- NULL

  
## PATHS ----


# Location of pre-processed data
input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, sample_name, "_", comparison_name, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)

# Path to a signature file
sig_path <- "data/import/Signatures/signatures_syngeneics_t100_Aug21.txt"
#sig_path <- "data/import/Signatures/signatures_pdac_t1000_May21.txt"


## INGEST DATA ----


# Seurat data
data_seurat <- readRDS(input_path)


## INGEST SIGNATURES ----


# Find genes abundant in this sample
gene_list <- seurat_select_abundant_genes_my(data_seurat, sct_threshold=0.5, spot_threshold=10, assay="SCT", slot="data", split_by="user.Sample_Name")

# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_length_min=1, sig_length_max=1000, ratio_threshold=0)

# Find all genes in signatures and their signature assignations
sig_invert_df <- invert_list_my(signature_list)


## WRANGLE DATA ----


# Extract coordinate meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

# Add user.Sample_Name to meta data if absent: happens in stand-alone samples
if(!"user.Sample_Name" %in% colnames(meta_df)){
  meta_df["user.Sample_Name"] = data_seurat@misc$user.Sample_Name
}

# Select spots for group1
barcodes_gr1 <- meta_df %>%
  dplyr::filter(!!rlang::sym(pathology_select_gr1) >= group_threshold &
                  grepl(sample_list_gr1, user.Sample_Name)) %>%
  dplyr::pull(Coordinate)

# Select spots for group1
barcodes_gr2 <- meta_df %>%
  dplyr::filter(!!rlang::sym(pathology_select_gr2) >= group_threshold &
                  grepl(sample_list_gr2, user.Sample_Name)) %>%
  dplyr::pull(Coordinate)

# Add group allocation back to the Seurat object
data_seurat@meta.data[which(rownames(data_seurat@meta.data) %in% barcodes_gr1), "DEA_Group"] <- gr1_name
data_seurat@meta.data[which(rownames(data_seurat@meta.data) %in% barcodes_gr2), "DEA_Group"] <- gr2_name
data_seurat@meta.data[which(is.na(data_seurat@meta.data[["DEA_Group"]])), "DEA_Group"] <- "Exclude"


## PLOT GROUP DISTRIBUTIONS ----


# Establish color scheme
if(is.null(cols)){
  cols <- c("darkblue", "red", "lightgrey")
  names(cols) <- c(gr1_name, gr2_name, "Exclude")
}

# Assign the largest spot size to all samples for better vis
data_seurat@misc$user.pt.size.factor <- max(data_seurat@meta.data$user.pt.size.factor)

## This logic below doesn't work because I screwed up syngeenics
# Assume an integrated file and find those samples that have one of the two categories of interest
#sample_list <- data_seurat@meta.data %>%
#  dplyr::filter(DEA_Group != "Exclude") %>%
#  dplyr::pull(user.Sample_Name) %>%
#  unique

# Plot spots selected as signature-high
#p <- spatial_dim_plot_my(data_seurat, group.by="DEA_Group", title=sample_name, cols=cols, images=sample_list)

p <- spatial_dim_plot_my(data_seurat, group.by="DEA_Group", title=sample_name, cols=cols)

# Write combo plot to file
filename <- paste0(output_path, "/spatial_", sample_name, "_", comparison_name, "_spotstatus")
write_plot2file_my(p, filename, num_row=1, num_col=length(names(data_seurat@images)))


## FIND MARKERS ----


# Subset Seurat data to chosen barcodes
data_seurat <- subset(data_seurat, cells=c(barcodes_gr1, barcodes_gr2))

# Define DEA parameters
params <- cluster_params_my()

# Minimum % in spots
params[["pct_min"]] <- 0.1
# Minimum FC difference to test
params[["logfc_threshold"]] <- 0.5
# FC signficance threshold
params[["sign_fc_threshold"]] <- 1

# Add latent parameters
params[["latent_vars"]] <- c("user.Sample_Name")
  
data_seurat@meta.data$DEA_Group <- factor(data_seurat@meta.data$DEA_Group, 
                                          levels=unique(data_seurat@meta.data$DEA_Group))

# Put the clustering of interest into seurat clusters
data_seurat <- Seurat::SetIdent(data_seurat, value="DEA_Group")

# Define output file name
filename_prefix <- paste0("dea_", sample_name, "_", comparison_name)
filename <- paste0(output_path, filename_prefix, "_markers_significant.txt")

# Find markers of each cluster against the rest
# Only two groups defined so find_all=TRUE still works
if(!file.exists(filename)){
  markers_df <- seurat_find_markers_my(data_seurat, assay=params[["assay"]], 
                                       find_all=TRUE, group.by="DEA_Group", 
                                       min.pct=params[["pct_min"]], 
                                       logfc.threshold=params[["logfc_threshold"]],
                                       test.use=params[["test_use"]], 
                                       latent.vars=params[["latent_vars"]])
  # Find signficant markers
  markers_filt_df <- marker_analysis_my(markers_df, params, "DEA_Group", 
                                        paste0(sample_name, "_", comparison_name), output_path, filename_prefix)
} else{
  markers_filt_df <- readr::read_delim(filename, delim="\t")
}


# Select top biomarkers for the group
markers_top_df <- markers_filt_df %>%
  dplyr::left_join(sig_invert_df, by="Symbol") %>%
  dplyr::filter(pct_1 >= pct1_threshold & 
                  direction == "UP") %>%
  dplyr::filter(!grepl("^MT-|^RP[SL]", Symbol))

readr::write_delim(markers_top_df, paste0(output_path, filename_prefix, "_mostrelevant.txt"), delim="\t")


# Select markers for plotting
markers_plot_df <- markers_top_df %>%
  # dplyr::filter(InSignature == TRUE) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_head(n=10)

# Merge top 10 list with the user defined gene list
gene_plot_list <- union(markers_plot_df$Symbol, intersect(markers_filt_df$Symbol, genes_user))

# Add a meta data column merging DEA_Group with user.Sample_Name
data_seurat@meta.data["DEAGroup_Sample"] <- paste(data_seurat@meta.data[["DEA_Group"]], data_seurat@meta.data[["user.Sample_Name"]], sep=": ")

# Create a dot plot of top biomarkers
p <- Seurat::DotPlot(data_seurat, features = gene_plot_list, group.by = "DEAGroup_Sample",
                     assay = "SCT", col.min = 0, col.max = 1, scale=TRUE, cluster.idents = FALSE) +
  Seurat::RotatedAxis()

# Identify the width of the plot
num_col <- round(log2(length(gene_plot_list)) / 2) + 2
if(num_col < 1){
  num_col <- 1
}

# Write to file
filename <- paste0(output_path, "/dea_dot_", sample_name, "_", comparison_name)
write_plot2file_my(p, filename, num_row=2, num_col=num_col)
