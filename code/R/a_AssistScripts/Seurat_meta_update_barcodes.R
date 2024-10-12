# Author: Anna Lyubetskaya. Date: 21-09-13

# Assist script for updating and augmenting data
# Update barcode meta data information from other Seurat objects in another Seurat object
# Gather a list of Seurat RDS objects, read a specific feature from them, update barcode meta data in another Seurat object, output to RDS

# Note: Proceed with caution! This script requires hand-holding because barcode names might require bespoke wrangling

# Another warning: Be careful before you overwrite this column and check for accidentally created .x and .y copied columns


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`


## PARAMETERS ----


# Sample / Cohort name
cohort_name <- "PDAC84_path12_epi_cca_sct"

# List of data_seurat@misc[["user.X"]] to update
# Note: user.Clustering is a flag to extract the real clustering resolution from the misc slot under user.Clustering
parameter_list <- c("integrated_snn_res.0.3")
# Assign different names to the selected parameters in the target object
parameter_rename <- c("epi_cca_sct.integrated_snn_res.0.3")

# Don't write the final RDS object
no_rds_output <- FALSE

# Filename suffix to remove
rds_suffix <- ""


## PATHS ----


# Input folder
input_path <- "XXXX"

# Output folder
output_path <- "XXXX"


## INGEST DATA ----


# Find all files following the regex
file_list <- dir(input_path, pattern=paste0(cohort_name, ".*.rds"), full.names=TRUE)

# Ingest the target Seurat object
data_target_seurat <- readRDS(output_path)


## WRANGLE DATA ----


data_list <- list()

# Cycle through input Seurat objects and collect barcode meta-data
for(filename in file_list){
  
    # Ingest the Seurat object
    data_seurat <- readRDS(filename)
    
    # Meta data parameter list copy
    parameter_list_loc <- parameter_list
    
    # user.Clustering is a flag to extract the real clustering resolution from the misc slot under user.Clustering
    if("user.Clustering" %in% parameter_list){
      parameter_list_loc[which(parameter_list_loc == "user.Clustering")] <- data_seurat@misc$user.Clustering
    }
    
    # Re-discover sample name
    sample_name <- data_seurat@misc$user.Sample_Name
    if(is.null(sample_name)){
      sample_name <- gsub(rds_suffix, "", gsub("^.+/", "", filename))
    }
    
    # Extract relevant meta data information
    data_list[[filename]] <- data_seurat@meta.data[parameter_list_loc] %>%
      tibble::rownames_to_column(var="Coordinate")
    
    # Rename the clustering resolution into a standard variable name
    if("user.Clustering" %in% parameter_list){
      colnames(data_list[[filename]]) <- gsub(data_seurat@misc$user.Clustering, "Clustering_Preferred", colnames(data_list[[filename]]))
    }
}

# Combine barcode meta data into a single tibble
data_df <- do.call(rbind, data_list)

# Rename parameters according to the user-defined list
if(!is.null(parameter_rename)){
  for(i in 1:length(parameter_list)){
    colnames(data_df) <- gsub(parameter_list[i], parameter_rename[i], colnames(data_df))
  }
}

# Are there any missing coordinates?
setdiff(rownames(data_target_seurat@meta.data), data_df$Coordinate)
intersect(rownames(data_target_seurat@meta.data), data_df$Coordinate)

# Add new meta data to the target Seurat object
data_target_seurat@meta.data <- data_target_seurat@meta.data %>%
  tibble::rownames_to_column(var="Coordinate") %>%
  dplyr::left_join(data_df, by="Coordinate") %>%
  tibble::column_to_rownames("Coordinate")

# Write the updated Seurat object to file
if(no_rds_output == FALSE){    
  saveRDS(data_target_seurat, file=output_path)
}
