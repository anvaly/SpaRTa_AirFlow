# Author: Anna Lyubetskaya. Date: 20-12-30

# Cluster and find markers for unannotated Seurat objects
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

# FindMarkers expects you to define your favorite resolution in the data_seurat@misc$user.Clustering slot
do_findmarkers <- FALSE

# If no covariates present, default to Wilcoxon for DEA
params[["test_use"]] <- "wilcox"  # wilcox; "LR", "negbinom", "poisson", or "MAST"
params[["assay"]] <- "SCT"  # Spatial or SCT

# Identify the best clustering resolution as the smallest resolution corresponding to the mean number of clusters
cluster_id <- TRUE

set.seed(531)


## PATHS ----


# Input folder
input_path <- "XXXX"

# Output folder
output_path <- "XXXX"
output_figs <- "XXXX"

# Seurat RDS files are tagged as follows
group_name <- "PDAC"
cohort_name <- "PDAC"

# Create output folders
dir.create(output_path, showWarnings = FALSE)
dir.create(output_figs, showWarnings = FALSE)


## INGEST DATA ----


# Find all files following the regex
file_list <- dir(input_path, pattern=paste0(cohort_name, ".*.rds"), full.names=TRUE)


## WRANGLE DATA ----


for(rds_file in file_list){
  
  print(rds_file)
  gc()
  
  
  ## Load data ----
  
  
  data_seurat <- readRDS(rds_file)
  
  # Define output file name
  sample_name <- data_seurat@misc$user.Sample_Name
  filename_rds <- paste0(output_path, sample_name, "_annotated.rds")
  
  if(!file.exists(filename_rds)){
    
    # Create a sample-specific output folder
    output_figs_sample <- paste0(output_figs, sample_name)
    dir.create(output_figs_sample, showWarnings = FALSE)
    
    
    ## Perform clustering analysis ----
    
    
    # Perform PCA, UMAP, FindNeighbors, FindClusters, and FindMarkers
    data_seurat <- cluster_analysis_my(data_seurat, params, sample_name, output_figs_sample)
    
    # Identify the best clustering resolution as the smallest resolution corresponding to the mean number of clusters
    if(cluster_id){
      clust_num <- sapply(params$resolution_list, function(x) length(unique(data_seurat@meta.data[[paste0("SCT_snn_res.", x)]])))
      data_seurat@misc$user.Clustering <- paste0("SCT_snn_res.", params$resolution_list[which.min(abs(clust_num - floor(mean(clust_num))))])
    }
    
    
    ## Find markers ----
    
    
    if(do_findmarkers == TRUE){
      data_seurat <- cluster_analysis_markers_my(data_seurat, params, sample_name, output_figs_sample)
    }
    
    
    ## Save annotated Seurat object ----
    
    
    # Save the full dataset to RDS
    saveRDS(data_seurat, file = filename_rds)
  }
}
