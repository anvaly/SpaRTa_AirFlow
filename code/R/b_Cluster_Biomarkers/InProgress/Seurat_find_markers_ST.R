# Author: Andy Kavran. Date: 2021-Dec-01

# Find markers for annotated Seurat objects
# 1. 
# 2. Cluster at a few resolutions and visualize the result via clustree
# 3. Perform FindMarkers for each resolution that produces a unique number of clusters


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

library(Seurat)
library(future)

setwd("/repos/P02567_TBIO-3021_10X_pilot/")

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
do_findmarkers <- TRUE

# If no covariates present, default to Wilcoxon for DEA
params[["test_use"]] <- "MAST"  # wilcox; "LR", "negbinom", "poisson", or "MAST"
params[["assay"]] <- "Spatial"  # Spatial or SCT

plan("multicore", workers = 20)
plan()
options(future.globals.maxSize = 10000 * 1024^2)

set.seed(823)


## PATHS ----


# Input folder
input_path <- "XXXX"

# Output folder
output_figs <- "XXXX"

# Seurat RDS files are tagged as follows
cohort_name <- "pdac_cohort_spotclean_harmony_integrated"


# Create output folders
dir.create(output_figs, showWarnings = FALSE)


## INGEST DATA ----


# Find all annotated files following the regex
file_list <- dir(input_path, pattern=paste0(cohort_name, ".*rds"), full.names=TRUE)


## WRANGLE DATA ----


for(rds_file in file_list){
  
  print(rds_file)
  ## Load data ----

  data_seurat <- readRDS(rds_file)
  
  # Define output file name
  sample_name <- data_seurat@misc$user.Sample_Name
  data_seurat@misc$user.Clustering <- "SCT_snn_res.0.5"
  # Create a sample-specific output folder if it doesn't exist
  output_figs_sample <- paste0(output_figs, sample_name, "/")
  dir.create(output_figs_sample, showWarnings = FALSE)
  
  
  ## Perform clustering analysis ----
  
  # Perform FindMarkers
  data_seurat <- cluster_analysis_markers_my(data_seurat, params, sample_name, output_figs_sample)
  
  ## Save annotated Seurat object ----
  
  # Save the full dataset to RDS
  saveRDS(data_seurat, file = rds_file)
}