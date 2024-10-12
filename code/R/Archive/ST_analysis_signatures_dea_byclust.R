# Author: Anna Lyubetskaya. Date: 21-07-21
# Perform DEA between different clusters that share a signature
# In this setting, integration is necessary for calling clusters; SCT is necessary for signatures; and DEA is done on raw data


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
# Integrated PDAC: HumanPanc_ROI1_FFPEv1_cca_sct; HumanPanc_FF_cca_sct
# Stand-alone PDAC : HumanPanc_ROI2_FFPE_D_Dec20_T; HumanPanc_ROI1_FFPE_A_Apr21; 
#                    HumanPanc_S001890_FFPE_C_Apr21; HumanPanc_S001891_FFPE_D_Apr21
sample_name <- "HumanPanc_S001891_FFPE_D_Apr21"
suffix <- "_annotated_pathology"  # "_annotated_pathology"

# Name for the signature group being plotted
sig_select <- "collisson.pdac.classical"

# SCT and spot number thresholds for gene filtering
sct_threshold <- 1
spot_threshold <- 5

# Assay and slot to use to calculate signature score
assay <- "SCT"
slot <- "data"


## PATHS ----


# Path to a signature file
sig_path <- "data/import/Signatures/signatures_pdac_t1000_May21.txt"

# Location of pre-processed data
input_path <- paste0("XXXX")
#input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, sample_name, "_", sig_select, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Seurat data
data_seurat <- readRDS(input_path)


## CALCULATE SIGNATURE SCORES ----


## Calculate signature scores and add them to Seurat meta data

# Find genes abundant in this sample
gene_list <- seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold,
                                             assay=assay, slot=slot)

# Load signatures and filter them down to only well represented genes
signature_init_list <- read_filter_signatures_my(sig_path, gene_list,
                                                 sig_length_min=1, sig_length_max=1000, ratio_threshold=0)

# Find all genes in signatures and their signature assignations
sig_invert_df <- invert_list_my(signature_init_list)


# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=c(sig_select),
                                            sig_length_min=1, sig_length_max=1000, ratio_threshold=0)

# Add signature scores to a seurat object
data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.", assay=assay)

# Find column names in the Seurat object - Seurat can rename original signatures
sig_names <- colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))])


## WRANGLE DATA ----


# User-defined resolution
resolution <- data_seurat@misc$user.Clustering

# Extract coordinate meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

# Add user.Sample_Name to meta data if absent: happens in stand-alone samples
if(!"user.Sample_Name" %in% colnames(meta_df)){
  meta_df["user.Sample_Name"] = data_seurat@misc$user.Sample_Name
}

# Global signature stats
sig_mean <- mean(meta_df[[sig_names]]) 
sig_sd <- sd(meta_df[[sig_names]]) 
sig_threshold <- sig_mean + sig_sd

# Calculate signature mean and SD for each sample in the cohort
meta_sig_df <- meta_df %>%
  dplyr::select(dplyr::all_of(c(resolution, sig_names))) %>%
  dplyr::group_by(!!rlang::sym(resolution)) %>%
  dplyr::summarise(Sig_Mean = mean(!!rlang::sym(sig_names)),
                   Sig_SD = sd(!!rlang::sym(sig_names))) %>%
  dplyr::mutate(Sig_IsHigh = Sig_Mean >= sig_threshold,
                Sig_IsHigh = ifelse(Sig_IsHigh == TRUE, paste0(Sig_IsHigh, "-", !!rlang::sym(resolution)), Sig_IsHigh))

# Add sample-specific thresholds to the meta data and select spots with high signature using the thresholds
meta_df <- meta_df %>%
  dplyr::inner_join(meta_sig_df, by=resolution)

print("Number of spots that have high (TRUE) signature vs low (FALSE):")
table(meta_df$Sig_IsHigh)

# Add updated meta data back to the Seurat object
data_seurat@meta.data <- meta_df %>%
  tibble::column_to_rownames("Coordinate")


## PLOT SIGNATURE DISTRIBUTIONS ----


# Plot signature score split by sample with global thresholds
filename <- paste0(output_path, "/bar_", sample_name, "_", sig_names, "_scores")
create_bar_plot_my(meta_sig_df, x_label=resolution, y_label="Sig_Mean", fill_label="Sig_IsHigh",
                   error_label="Sig_SD", reorder_x=TRUE, filename=filename)


# Plot spots selected as signature-high
p <- spatial_dim_plot_my(data_seurat, group.by="Sig_IsHigh", title=sample_name)

# Write combo plot to file
filename <- paste0(output_path, "/spatial_", sample_name, "_", sig_names, "_spotstatus")
write_plot2file_my(p, filename, num_row=1, num_col=length(names(data_seurat@images)))


## FIND MARKERS ----


# Find clusters with high signature scores
clust_list <- meta_sig_df %>%
  dplyr::filter(grepl("TRUE", Sig_IsHigh)) %>%
  dplyr::pull(resolution) %>%
  as.character()

if(length(clust_list) >= 2){
  # Find barcodes representing clusters of interest
  barcode_list <- meta_df %>%
    dplyr::filter(!!rlang::sym(resolution) %in% clust_list) %>%
    dplyr::pull(Coordinate)
  
  # Subset Seurat object to only barcodes in clusters of interest
  data_subset_seurat <- subset(data_seurat, cell=barcode_list)
  
  
  # Define DEA parameters
  params <- cluster_params_my()
  
  # Add latent parameters
  params[["latent_vars"]] <- c("user.Sample_Name")
  
  # Put the clustering of interest into seurat clusters
  data_subset_seurat <- Seurat::SetIdent(data_subset_seurat, value=resolution)
  
  # Find markers of each cluster against the rest
  # Only two groups defined so find_all=TRUE still works
  markers_df <- seurat_find_markers_my(data_subset_seurat, assay=params[["assay"]], 
                                       find_all=TRUE, group.by=resolution, 
                                       min.pct=params[["pct_min"]], 
                                       logfc.threshold=params[["logfc_threshold"]],
                                       test.use=params[["test_use"]], 
                                       latent.vars=params[["latent_vars"]])
  
  # Find signficant markers
  filename_prefix <- paste0("dea_", sample_name, "_", sig_names)
  markers_filt_df <- marker_analysis_my(markers_df, params, "Sig_IsHigh", 
                                        paste0(sample_name, "_", sig_names), output_path, filename_prefix)
}
