# Andy Kavran, 2022-August-18
# Re-normalize a set of Seurat objects
# This script over-writes original input RDS objects!


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

library(Seurat)


## DEFINE PATHS AND PARAMETERS ----


# Regex to select samples
cohort_regex <- "PDAC"

# Input and output path
input_objects <- "XXXX"
output_path <- input_objects

# Read all files in the input folder
file_list <- dir(input_objects, full.names = TRUE)[grepl(cohort_regex, dir(input_objects))]

vars_to_regress <- NULL  # c("mito_percent", "ribo_percent")


## RUN THROUGH LIST OF OBJECTS ----


for(rds_file in file_list){
  
  print(rds_file)
  
  
  ## READ SEURAT ----
  
  
  # Open a connection to the RDS object
  con <- gzfile(rds_file)
  
  # Ingest the Seurat object
  data_seurat <- readRDS(con)
  
  # Close the connection to be able to overwrite
  close(con)
  
  
  ## NORMALIZE SEURAT ----
  
  
  # Perform SCTransform normalization
  data_seurat <- Seurat::SCTransform(data_seurat, assay="Spatial", vars.to.regress=vars_to_regress, variable.features.n=8000,
                                     return.only.var.genes=FALSE, verbose=FALSE, vst.flavor="v2")

  
  ## WRITE SEURAT ----
  
  
  # Write the updated Seurat object
  saveRDS(data_seurat, file=rds_file)
}
