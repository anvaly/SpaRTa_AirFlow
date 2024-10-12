# Author: Anna Lyubetskaya. Date: 22-08-23
# Investigate various annotation factors relative to each other
# Perform DEA of spots with a set of properties to the rest
# Spots are selected using signature scores, pathology, and clusters
# In this setting, the more data the better


## ENVIRONMENT ----


# Install MAST
BiocManager::install("MAST")

library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_signatures.R")

source("code/R/b_Cluster_Biomarkers/Seurat_cluster_utils.R")


## PARAMETERS ----


# Name of the analysis to use in folder/file names
run_name <- "T_v_NT"

# Path to processed Seurat data
sample_name <- "PDAC84_path12_merge"
sample_exclude <- NULL

# Name for the signature group being plotted
sig_select <- c("PDAC.collisson.classical", "PDAC.moffitt.basal")

# Name of the pathology field to filter by groups
pathology_select <- c("Pathology.Epithelium.percent", "Pathology.Benign Glands.percent", "Pathology.Luminal Debris.percent")
pathology_name <- "Pathology.Tumor"

# SCT and spot number thresholds for gene filtering
sct_threshold <- 0.5
spot_threshold <- 100

# Assay and slot to use to calculate signature score
assay <- "SCT"
slot <- "data"

# Feature # in a spot threshold
feature_threshold <- NULL

# Subsample 2nd group?
sample_group2 <- TRUE

# Groups to compare
group1_name <- "tumor"
group1 <- c(paste0("sig.", sig_select[1], "_", pathology_name, "_TumorTrue"), paste0("sig.", sig_select[2], "_", pathology_name, "_TumorTrue"))
print(group1)

group1_exclude_name <- "other"
group1_exclude <- NULL  # paste0("sig.", sig_select[1], "_", pathology_name, "_TumorTrue")
print(group1_exclude)

group2_name <- "non_tumor"
group2 <- c(paste0("sig.", sig_select[1], "_", pathology_name, "_TumorTrue"), paste0("sig.", sig_select[2], "_", pathology_name, "_TumorFalse"))
print(group2)


## PATHS ----


# Path to a signature file
sig_path <- "data/import/Signatures/signatures_220210.txt"

# Location of pre-processed data
input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, run_name, "_", sample_name, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Seurat data
data_seurat <- readRDS(input_path)

# Remove samples if necessary
if(!is.null(sample_exclude)){
  barcode_list <- rownames(data_seurat@meta.data[which(!data_seurat@meta.data[["user.Sample_Name"]] %in% sample_exclude),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# Remove spots if necessary
if(!is.null(feature_threshold)){
  barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data$nFeature_SCT >= feature_threshold),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}


## CALCULATE TARGET SIGNATURE SCORES ----


## Calculate signature scores and add them to Seurat meta data

# Find genes abundant in this sample
gene_list <- unique(seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold, 
                                                    assay=assay, slot=slot))

# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select, 
                                            sig_length_min=1, sig_length_max=1000, ratio_threshold=0)

# Load signatures for annotation
signature_ann_list <- read_filter_signatures_my(sig_path, gene_list, 
                                                sig_length_min=1, sig_length_max=1000, ratio_threshold=0)

# Find all genes in signatures and their signature assignations
sig_invert_df <- invert_list_my(signature_ann_list)

# Add signature scores to a Seurat object
data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.", assay=assay)

# Find column names in the Seurat object - Seurat can rename original signatures
sig_names <- colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))])


## WRANGLE META DATA ----


# Extract coordinate meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

# Create a joint pathology label
meta_df[[pathology_name]] <- matrixStats::rowMaxs(as.matrix(meta_df[pathology_select]))

# Add user.Sample_Name to meta data if absent: happens in stand-alone samples
if(!"user.Sample_Name" %in% colnames(meta_df)){
  meta_df["user.Sample_Name"] = data_seurat@misc$user.Sample_Name
}

# Threshold signatures
for(s in sig_names){
  # Global signature stats
  sig_median_global <- median(meta_df[[s]]) 
  sig_sd_global <- sd(meta_df[[s]]) 
  sig_threshold_h <- sig_median_global + sig_sd_global * 2
  sig_threshold_l <- sig_median_global
  
  cat(s, sig_median_global, sig_sd_global, sig_threshold_h, sig_threshold_l, "\n")
  
  # New column names
  name1 <- paste0(s, "_TumorTrue")
  name2 <- paste0(s, "_TumorFalse")
  
  # ID tumor categories
  meta_df <- meta_df %>%
    dplyr::mutate(PathSig1 = !!rlang::sym(pathology_name) >= 50 & !!rlang::sym(s) >= sig_threshold_h,
                  PathSig2 = !!rlang::sym(pathology_name) <= 25 & !!rlang::sym(s) <= sig_threshold_l)
  
  colnames(meta_df) <- gsub("^PathSig2$", paste0(s, "_", pathology_name, "_TumorFalse"), colnames(meta_df))
  colnames(meta_df) <- gsub("^PathSig1$", paste0(s, "_", pathology_name, "_TumorTrue"), colnames(meta_df))
}


## DEFINE GROUPS ----


# Define categorical labels of groups to compare
group1_indx <- which(rowSums(meta_df[group1]) > 0)

# Exclude a set if user-defined
if(!is.null(group1_exclude)){
  group1_indx <- setdiff(which(rowSums(meta_df[group1]) > 0), which(rowSums(meta_df[group1]) > 0))
}

group2_indx <- which(rowSums(meta_df[group2]) > 0)

# Subsample group2 if it's too large
if(sample_group2 == TRUE){
  group2_indx <- sample(group2_indx, length(group1_indx), replace=FALSE)
}

meta_df["Comparison"] <- NA
meta_df[group1_indx, "Comparison"] <- group1_name
meta_df[group2_indx, "Comparison"] <- group2_name

table(meta_df["Comparison"])

# Update the Seurat object meta data
data_seurat@meta.data <- meta_df %>%
  tibble::column_to_rownames("Coordinate")

# Subset the Seurat object to groups of interest only
data_seurat <- subset(data_seurat, cells=meta_df %>%
                        dplyr::filter(!is.na(Comparison)) %>%
                        dplyr::pull(Coordinate))


## FIND MARKERS ----


# Specific comparison to test
cluster_pairs <- as.matrix(c(group1_name, group2_name))

# Define DEA parameters
params <- cluster_params_my()

# Add latent parameters
params[["latent_vars"]] <- c("user.Sample_Name")

# Biomarker method and data source
params[["test_use"]] <- "MAST"
params[["assay"]] <- "Spatial"

# Minimum % in spots
params[["pct_min"]] <- 0.25
# Minimum FC difference to test
params[["logfc_threshold"]] <- 0.5

# Put the clustering of interest into seurat clusters
data_seurat <- Seurat::SetIdent(data_seurat, value="Comparison")

# Find markers of each cluster against the rest
# Only two groups defined so find_all=TRUE still works
markers_df <- seurat_find_markers_my(data_seurat, assay=params[["assay"]], 
                                     find_all=FALSE, group.by="Comparison", 
                                     min.pct=params[["pct_min"]], 
                                     logfc.threshold=params[["logfc_threshold"]],
                                     test.use=params[["test_use"]], 
                                     latent.vars=params[["latent_vars"]],
                                     cluster_pairs=cluster_pairs)

# Find significant markers
filename_prefix <- paste0("dea_", sample_name, "_", run_name)
markers_filt_df <- marker_analysis_my(markers_df, params, "Comparison", 
                                      paste0(sample_name, "_", run_name), output_path, filename_prefix) %>%
  dplyr::left_join(sig_invert_df, by="Symbol") %>%
  dplyr::filter(!grepl("^MT-|^RP[SL]", Symbol)) %>%
  tidyr::replace_na(list(InSignature = FALSE))

readr::write_delim(markers_filt_df, paste0(output_path, filename_prefix, "_markers_significant_ann.txt"), delim="\t")

# Select top biomarkers for the group
markers_top_df <- markers_filt_df %>%
  dplyr::filter(direction == "UP" & !grepl("^A[LPC]\\d+{6}", Symbol) & 
                  pct_1 >= 75 & pct_2 <= 40 & avg_logFC >= 100 & p_val_adj_neg_log10 == Inf)

print(markers_top_df$Symbol)

readr::write_delim(markers_top_df, paste0(output_path, filename_prefix, "_mostrelevant.txt"), delim="\t")


## VISUALIZATIONS ----


# Dot plot of the most relevant genes
# Before: sapply(1:nrow(meta_df), function(x) paste0(data_seurat@meta.data[x, c("Comparison", "user.Sample_Name")], collapse=" "))
data_seurat@meta.data[["Group"]] <- data_seurat@meta.data[["user.Sample_Name"]]  

# Identify the width of the plot
num_col <- round(log2(length(markers_top_df$Symbol)) / 2) + 2
if(num_col < 1){
  num_col <- 1
}

# Go through every sample and create a dotplot of most significant genes
for(s in unique(data_seurat@meta.data$user.Sample_Name)){
  # Identify barcodes of interest
  cells_select <- rownames(data_seurat@meta.data[which(data_seurat@meta.data$user.Sample_Name == s),])
  data_seurat_loc <- subset(data_seurat, cells=cells_select)
  
  # Create a dot plot for a sample
  p <- Seurat::DotPlot(data_seurat_loc, features = markers_top_df$Symbol, group.by = "Group",
                       assay = assay, col.min = 0, col.max = 0.5, scale=TRUE) +
    Seurat::RotatedAxis()
  
  filename <- paste0(output_path, "/dea_dot_", s, "_", run_name)
  write_plot2file_my(p, filename, num_row=1, num_col=num_col)
}

# Create a dot plot for the integrated cohort
p1 <- Seurat::DotPlot(data_seurat, features = markers_top_df$Symbol, group.by = "Comparison",
                      assay = assay, col.min = 0, col.max = 0.5, scale=TRUE) +
  Seurat::RotatedAxis()

filename <- paste0(output_path, "/dea_dot_", run_name)
write_plot2file_my(p1, filename, num_row=1, num_col=num_col)

# Format the filtered data for plotting
markers_filt_plot_df <- markers_filt_df %>%
  dplyr::mutate(avg_logFC = log10(2^avg_logFC),
                random = round(runif(nrow(markers_filt_df), 1, 60)),
                p_val_adj_neg_log10 = ifelse(p_val_adj_neg_log10 == Inf, 340 + random, p_val_adj_neg_log10),
                labels = ifelse(Symbol %in% markers_top_df$Symbol, Symbol, NA),
                size = round(pct_1/20, 1),
                TopHit = Symbol %in% markers_top_df$Symbol) %>%
  dplyr::filter(direction == "UP" & pct_1 >= 60 & avg_logFC >= 1)

# Create a scatter plot illustrating significant gene distribution
p2 <- create_scatter_plot_my(markers_filt_plot_df, 
                             x_label="avg_logFC", y_label="p_val_adj_neg_log10", fill_label="TopHit",
                             shape=21, size=6, filename=filename, dot_labels="labels", do_fit=NULL,
                             labels=c("Fold change, log10", "Negative adjusted p-value, log10", "Tumor v Non-tumor in PDAC")) + 
  ggplot2::geom_hline(yintercept=320, color="black", linetype="dashed", size=0.2)

filename <- paste0(output_path, "/dea_scatter_pval_", run_name)
write_plot2file_my(p2, filename, num_row=2, num_col=2)


# Format the filtered data for plotting
markers_filt_plot_df <- markers_filt_df %>%
  dplyr::mutate(labels = ifelse(Symbol %in% markers_top_df$Symbol, Symbol, NA),
                TopHit = Symbol %in% markers_top_df$Symbol) %>%
  dplyr::filter(direction == "UP" & p_val_adj_neg_log10 >= 300 & avg_logFC >= 2 & pct_1 >= 75)

# Create a scatter plot illustrating significant gene distribution
p3 <- create_scatter_plot_my(markers_filt_plot_df, 
                             x_label="avg_logFC", y_label="pct_1", fill_label="InSignature",
                             shape=21, size=6, filename=filename, dot_labels="labels", do_fit=NULL,
                             labels=c("Fold change, log2", "Percent spots expressing the gene in tumor", "Tumor v Non-tumor in PDAC")) + 
  ggplot2::geom_hline(yintercept=90, color="black", linetype="dashed", size=0.2)

filename <- paste0(output_path, "/dea_scatter_pct1_", run_name)
write_plot2file_my(p3, filename, num_row=2, num_col=2)
