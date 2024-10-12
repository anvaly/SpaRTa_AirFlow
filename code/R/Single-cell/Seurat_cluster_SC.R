# Author: Anna Lyubetskaya. Date: 20-12-30

# Add various data to the Seurat objects
# 1. Pick the number of PCs to use based on variance explained
# 2. Cluster at a few resolutions and visualize the result via clustree
# 3. Perform FindMarkers for each resolution that produces a unique number of clusters


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_dea.R")
source("code/R/Utils/utils_10X_matrix.R")

source("code/R/b_Cluster_Biomarkers/Seurat_cluster_utils.R")


## PARAMETERS ----


# Setup parameters
params <- cluster_params_my()


## PATHS ----


# Input folder
input_path <- "XXXX"

# Output folder
output_path <- "XXXX"
output_figs <- "XXXX"

# Seurat RDS files are tagged as follows
data_regex_files <- "_all"

# Create output folders
dir.create(output_path, showWarnings = FALSE)
dir.create(output_figs, showWarnings = FALSE)


## INGEST DATA ----


# Find all files following the regex
file_list <- dir(input_path, pattern=paste0(data_regex_files, ".rds"), full.names=TRUE)


## WRANGLE DATA ----


for(rds_file in file_list){
  
  name <- gsub(data_regex_files, "", gsub("^.+/|\\.[^\\.]+$", "", rds_file))
  filename_rds <- paste0(output_path, name, "_annotated.rds")
  
  print(name)
  
  if(!file.exists(filename_rds)){
    
    
    ## Load data ----
    
    data_seurat <- readRDS(rds_file)
    
    # Create a sample-specific output folder
    output_figs_sample <- paste0(output_figs, name)
    dir.create(output_figs_sample, showWarnings = FALSE)
    
    
    ## Perform clustering analysis ----
    
    # Perform PCA, UMAP, FindNeighbors, FindClusters, and FindMarkers
    data_seurat <- cluster_analysis_my(data_seurat, params, name, output_figs_sample)
    
    
    ## Save annotated Seurat object ----
    
    # Save the full dataset to RDS
    saveRDS(data_seurat, file = filename_rds)
  }
}
