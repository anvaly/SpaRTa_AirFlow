# Author: Anna Lyubetskaya. Date: 23-10-02
# Subset a cohort using pathology and cluster labels and do a DEA of leftover clusters


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
# TLS_LN, Mac_in_Fib, Tcell_in_Fib, Bcell_in_Fib
run_name <- "Bcell_in_Fib"

# Path to processed Seurat data
sample_name <- "PDAC108_path14_5K_harmony"
sample_exclude <- NULL

# Name of the pathology field to filter by groups
# c("LymphNode", "TLSAggregate", "TLSImmature", "TLSMature") or c("NonEpi")
pathology_select <- c("NonEpi")

# Select specific clusters
cluster_select <- NULL

# Feature # in a spot threshold
feature_threshold <- NULL

# Cluster by the following columns
# Pathology.Group, sig.PDAC.P19.Macrophage, sig.PDAC.P19.Tcell, sig.PDAC.P19.Bcell
cluster_by_col <- "sig.PDAC.P19.Bcell"


## PATHS ----


# Path to a signature file
sig_path <- "data/import/Signatures/signatures_pdac_20230720.txt"

# Location of pre-processed data
input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, run_name, "_", sample_name, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


# A file of additional barcode labels
meta_data_labels <- "XXXX"

# Log output file
log_path <- paste0(output_path, "log.txt")
write("Cluster\tSpot_num\tSample_num\tRegion_num\tBlock_num", log_path, append=FALSE)


## INGEST DATA ----


# Seurat data
data_seurat <- readRDS(input_path)

# Add additional meta data labels and do bespoke subsetting
if(!is.null(meta_data_labels)){
  meta_df <- readr::read_delim(meta_data_labels, delim="\t") %>%
    dplyr::filter(sig.PDAC.P19.Fibroblast == 1 &
                    sig.PDAC.P19.Macrophage == 0 &
                    sig.PDAC.P19.Tcell == 0 &
                    # sig.PDAC.P19.Bcell == 0 &
                    sig.PDAC.collisson.classical == 0 &
                    sig.PDAC.moffitt.basal == 0 &
                    sig.PDAC.P19.Acinar == 0 &
                    sig.PDAC.P19.Ductal_1 == 0 &
                    sig.PDAC.P19.Endocrine == 0 &
                    sig.PDAC.P19.Endothelial == 0 &
                    sig.PDAC.P19.Stellate == 0)
  
  data_seurat <- subset(data_seurat, cells=meta_df$Coordinate)
  
  data_seurat@meta.data <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate") %>%
    dplyr::inner_join(meta_df, by="Coordinate") %>%
    tibble::column_to_rownames("Coordinate")
}

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

# Remove clusters if necessary
if(!is.null(cluster_select)){
  barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data[[resolution]] %in% cluster_select),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# Remove pathology groups if necessary
if(!is.null(pathology_select)){
  barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data[["Pathology.Group"]] %in% pathology_select),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

gc()


# Print some stats
for(p in unique(data_seurat@meta.data[[cluster_by_col]])){
  data_loc_df <- data_seurat@meta.data %>%
    dplyr::filter(!!rlang::sym(cluster_by_col) == p)
  
  write(paste0(c(p,
                 length(data_loc_df$user.Sample_Name),
                 length(unique(data_loc_df$user.Sample_Name)),
                 length(unique(data_loc_df$user.Region_ID)),
                 length(unique(data_loc_df$user.Block_ID))), collapse="\t"), 
        log_path, append=TRUE)
}


## EXTRACT META DATA ----


meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")
readr::write_delim(meta_df, paste0(output_path, "meta_data.txt"), delim="\t")


## INGEST SIGNATURES ----


# Load signatures for annotation
signature_list <- read_filter_signatures_my(sig_path)

# Find all genes in signatures and their signature assignations
sig_invert_df <- invert_list_my(signature_list)


## FIND MARKERS ----


# Define DEA parameters
params <- cluster_params_my()

# Add latent parameters
params[["latent_vars"]] <- c("user.Sample_Name")

# Biomarker method and data source
params[["test_use"]] <- "wilcox"
params[["assay"]] <- "SCT"

# Minimum % in spots
params[["pct_min"]] <- 0.1
# Minimum FC difference to test
params[["logfc_threshold"]] <- 0.25

# Put the clustering of interest into seurat clusters
data_seurat <- Seurat::SetIdent(data_seurat, value=cluster_by_col)


# Find markers of each cluster against the rest
markers_df <- seurat_find_markers_my(data_seurat, 
                                     assay=params[["assay"]], 
                                     group.by=cluster_by_col, 
                                     min.pct=params[["pct_min"]], 
                                     logfc.threshold=params[["logfc_threshold"]],
                                     test.use=params[["test_use"]], 
                                     latent.vars=params[["latent_vars"]])

# Find significant markers
filename_prefix <- paste0("markers_signficant_dea_", sample_name, "_", run_name)
markers_filt_df <- marker_analysis_my(markers_df, params, cluster_by_col, 
                                      paste0(sample_name, "_", run_name), output_path, filename_prefix) %>%
  dplyr::left_join(sig_invert_df, by="Symbol") %>%
  dplyr::filter(!grepl("^MT-|^RP[SL]", Symbol)) %>%
  tidyr::replace_na(list(InSignature = FALSE))

readr::write_delim(markers_filt_df, paste0(output_path, filename_prefix, "_markers_significant_ann.txt"), delim="\t")

# Select top biomarkers for the group
markers_top_df <- markers_filt_df %>%
  dplyr::filter(direction == "UP" & !grepl("^A[LPC]\\d+{6}", Symbol))

# Count re-curring markers
marker_count_df <- markers_top_df %>% 
  dplyr::group_by(Symbol) %>% 
  dplyr::summarise(clust_count = dplyr::n_distinct(cluster)) %>% 
  dplyr::arrange(desc(clust_count))

# Flag repetitive top genes
markers_top_df <- markers_top_df %>%
  dplyr::inner_join(marker_count_df, by="Symbol")

print(markers_top_df$Symbol)

readr::write_delim(markers_top_df, paste0(output_path, filename_prefix, "_mostrelevant.txt"), delim="\t")


## VISUALIZATIONS ----


## Dot plot of shared biomarkers


# Select data for a cluster
markers_top_loc_df <- markers_top_df %>%
  dplyr::filter(Symbol %in% (marker_count_df %>%
                               dplyr::filter(clust_count > 1) %>%
                               dplyr::pull(Symbol)))

if(nrow(markers_top_loc_df) > 0){
  # Identify the width of the plot
  num_col <- round(log2(length(markers_top_loc_df$Symbol)))
  if(num_col < 1){
    num_col <- 1
  }
  
  # Create a dot plot for the integrated cohort
  p1 <- Seurat::DotPlot(data_seurat, features=unique(markers_top_loc_df$Symbol), group.by=cluster_by_col,
                        assay="SCT", col.min=0, col.max=0.5, scale=TRUE) +
    Seurat::RotatedAxis()
  
  filename <- paste0(output_path, "/dea_dot_shared")
  write_plot2file_my(p1, filename, num_row=1, num_col=num_col)
}


## Dot plot of unique biomarkers


# Create dotplots for each cluster
for(cl in unique(markers_top_df$cluster)){
  
  # Select data for a cluster
  markers_top_loc_df <- markers_top_df %>%
    dplyr::filter(cluster == cl,
                  Symbol %in% (marker_count_df %>%
                                 dplyr::filter(clust_count <= 1) %>%
                                 dplyr::pull(Symbol)))
  
  if(nrow(markers_top_loc_df) > 0){
    # Identify the width of the plot
    num_col <- round(log2(length(markers_top_loc_df$Symbol)))
    if(num_col < 1){
      num_col <- 1
    }
    
    # Create a dot plot for the integrated cohort
    p1 <- Seurat::DotPlot(data_seurat, features=unique(markers_top_loc_df$Symbol), group.by=cluster_by_col,
                          assay="SCT", col.min=0, col.max=0.5, scale=TRUE) +
      Seurat::RotatedAxis()
    
    filename <- paste0(output_path, "/dea_dot_", cl)
    write_plot2file_my(p1, filename, num_row=1, num_col=num_col)
  }
}


# Cleanup images for plotting
data_seurat@images <- data_seurat@images[unique(data_seurat@meta.data$user.Sample_Name)]

# Create a PCA, UMAP, and spatial plots
output_name <- paste0(output_path, "/cluster_", filename_prefix)
Seurat_pca_umap_spatial_my(data_seurat, cluster_by_col, output_name)
