# Author: Anna Lyubetskaya. Date: 20-12-30

# decontX implementation


## SETUP ENVIRONMENT ----


if(!"celda" %in% rownames(installed.packages())){
  BiocManager::install("celda")
}

library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_signatures.R")

source("code/R/Utils/utils_10X_matrix.R")


## PARAMETERS ----


# Cohort name
cohort_name <- "HumanPanc"
# Seurat RDS files are tagged as follows
cohort_expr <- "HumanPanc.*"
# User-defined resolution
resolution <- "SCT_snn_res.0.6"
resolutions <- c("SCT_snn_res.0.1", "SCT_snn_res.0.2", "SCT_snn_res.0.3", "SCT_snn_res.0.4", "SCT_snn_res.0.5", 
                 "SCT_snn_res.0.6", "SCT_snn_res.0.7", "SCT_snn_res.0.8", "SCT_snn_res.0.9")

## PATHS ----


# Input folder
input_path <- "XXXX"

# Output folder
output_path_init <- "XXXX"

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)


## INGEST DATA ----


# Find all files following the regex
file_list <- dir(input_path, pattern=paste0(cohort_expr, ".rds"), full.names=TRUE)


## WRANGLE DATA ----


for(rds_file in file_list){
  
  # Read the Seurat object
  data_seurat <- readRDS(rds_file)
  
  # Sample name
  sample_name <- data_seurat@misc$user.Sample_Name
  
  # Create a sample-specific output folder
  output_path <- paste0(output_path_init, sample_name, "/")
  dir.create(output_path, showWarnings = FALSE)
  
  # Run or read decontX data
  filename <- paste0(output_path, "decontX_result.RDS")
  if(!file.exists(filename)){
    # Run decontX on the data and save to RDS
    decont_list <- list()
    for(res in resolutions){
      decont_list[[res]] <- celda::decontX(data_seurat@assays$Spatial@counts, z=data_seurat@meta.data[[res]])
    }
  
    saveRDS(decont_list, filename)
  } else{
    decont_list <- readRDS(filename)
  }
  
  # Create a tibble of all contamination values
  contamination_wide_df <- tibble::as_tibble(data.frame(lapply(decont_list, function(x) x$contamination)))
  colnames(contamination_wide_df) <- resolutions
  
  contamination_df <- contamination_wide_df %>%
    dplyr::mutate(Coordinate = colnames(data_seurat)) %>%
    dplyr::relocate(Coordinate) %>%
    df_wide2long_my(key="Resolution", val="ContaminationPercent") %>%
    dplyr::mutate(ContaminationPercent = round(ContaminationPercent * 100),
                  Negligible = ContaminationPercent <= 5)
  
  # Plot contamination values across resolutions
  filename <- paste0(output_path, "hist_decont")
  create_hist_plot_my(contamination_df, x_label="ContaminationPercent", fill_label="Negligible", 
                      facet_var=c("Resolution", "fixed"), intercept=c(0, 100), binwidth=1, filename=filename, labels=NULL, 
                      add_density=FALSE, log_scale=TRUE)
  
  # Plot contamination on UMAP
  # celda::plotDecontXContamination(decont_list[[1]])
}
