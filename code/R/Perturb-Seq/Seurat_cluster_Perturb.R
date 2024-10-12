# Author: Anna Lyubetskaya. Date: 20-12-30

# Add various data to the Seurat objects
# 1. Pick the number of PCs to use based on variance explained
# 2. Cluster at a few resolutions and visualize the result via clustree
# 3. Perform FindMarkers for each resolution that produces a unique number of clusters


## INSTALL PACKAGES ABSENT IN DOMINO ----


library(clustree)


## SETUP ENVIRONMENT ----


# The library doesn't work without loading because it has dependency on ggraph
library(clustree)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_dea.R")
source("code/R/Utils/utils_10X_matrix.R")


## PARAMETERS ----


# PCA, minimum variance explained
pca_var_threshold <- 0.2
# List of resolutions to test in clustering
resolution_list <- c(0.8)
# All cluster resolutions tested
cluster_names <- c(paste0("SCT_snn_res.", resolution_list), "CRISPRGeneCategory")
de_reference <- c("NON-TARGETING", "background")

# Minimum cluster size
cluster_num_min <- 3
# Maximun cluster size
cluster_num_max <- 200
# Minimum number of spots in a cluster
cluster_size_min <- 5

# Minimum % in spots
pct_min <- 0.5
# Minimum FC difference to test
logfc_threshold <- 0.25

# FC signficance threshold
sign_fc_threshold <- 0.5
# Signficance p-value threshold
sign_pval_adj_neglog10_threshold <- 2

# Color scheme
col_type <- "jet"


## PATHS ----


# Input folder
input_path <- "XXXX"

# Output folder
output_path <- "XXXX"
output_figs <- "XXXX"

# Seurat RDS files are tagged as follows
data_regex_files <- "_all.rds"

# Create output folders
dir.create(output_figs, showWarnings = FALSE)


## INGEST DATA ----


# Find all files following the regex
file_list <- dir(input_path, pattern=data_regex_files, full.names=TRUE)


## WRANGLE DATA ----


for(rds_file in file_list){
  
  name <- gsub("_all", "", gsub("^.+/|\\.[^\\.]+$", "", rds_file))
  filename_rds <- paste0(output_path, name, "_annotated.rds")
  
  print(name)
  
  if(!file.exists(filename_rds)){
    
    
    ## Load data ----
    
    data_seurat <- readRDS(rds_file)
    
    data_seurat@meta.data$CRISPRGeneCategory <- gsub("\\-\\d+$", "", data_seurat@meta.data$CRISPRGeneCategory)
    
    # Create a sample-specific output folder
    output_figs_sample <- paste0(output_figs, name)
    dir.create(output_figs_sample, showWarnings = FALSE)
    
    
    ## PCA, UMAP, FindNeighbors, FindClusters ----
    
    # Perform PCA, visualize variance explained by each PC, and report the number of useful PCs
    filename <- paste0(output_figs_sample, "/pca_barplot_", name, "")
    seurat_pca <- seurat_pca_my(data_seurat, var_threshold=pca_var_threshold, output_file=filename)
    
    # Perform UMAP and find neighbors on previously PCA Seurat data
    data_seurat <- seurat_umap_nb_my(seurat_pca$data, num_dimensions=seurat_pca$num_pcs)
    
}
