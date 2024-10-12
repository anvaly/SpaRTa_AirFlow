# Simple Merge Script
# Takes in a list of seurat objects and merges them.
# SCTransform all the objects, then writes to disk.

## Environment ----
# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

# Update Seurat
if(packageVersion("Seurat") < 4.1){
  install.packages("Seurat")
}
if(packageVersion("glmGamPoi") <1.6){
  install.packages("glmGamPoi")
}
library(Seurat)


source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_10X_cohort.R")

## Define paths and parms ----

seurat_path <- c("XXXX")

output_path <- "XXXX"
file_out_name <- "full_cohort_spotclean_pathv9_SCT_v2"

dir.create(output_path)

misc_field_list <- c("user.Sample_Name", "user.Group_Name", "user.Tissue",
                     "user.Date", "user.Area", "user.Project_ID", "user.Primary",
                     "user.Protocol", "user.Organism", "user.pt.size.factor",
                     "user.Clustering", "user.Slide", "Pathology.Tissue", "user.KRAS",
                     "user.Block_ID", "user.Concentriq_Repo_ID", "user.Concentriq_Image_ID",
                     "user.Image_Filename")
cohort_regex <- "PDAC"
n_var_features <- 7000
file_list <-  file_list <- dir(seurat_path, pattern = cohort_regex, full.names = TRUE)

## Merge Objects ----

# Ingest a set of RDS Seurat objects
start_time <- proc.time()
seurat_list <- read_rds_list_simple_my(file_list) #full cohort = 14 minutes
proc.time() - start_time

# Remove old SCT model and clustering
for(name in names(seurat_list)){
  sample_name <- seurat_list[[name]]@misc$user.Sample_Name #get sample name
  seurat_list[[name]]@active.assay <- "Spatial"
  seurat_list[[name]][["SCT"]] <- NULL
  names(seurat_list[[name]]@images) <- sample_name # rename image slice with sample name
  seurat_list[[name]]@meta.data <- seurat_list[[name]]@meta.data[which(!grepl("SCT_snn_res.", names(seurat_list[[name]]@meta.data)))] # remove old clustering
}


# Move misc parameters to meta data slot
seurat_list <- seurat_misc_to_meta(seurat_list, misc_field_list=misc_field_list)

# Merge individual objects together
start_time <- proc.time()
merged_seurat_data <- seurat_merge_my(seurat_list)
proc.time() - start_time # whole cohort 35 minutes

# SCTransform the merged dataset
merged_seurat_data <- Seurat::SCTransform(merged_seurat_data, assay="Spatial", variable.features.n = n_var_features,
                                          return.only.var.genes = TRUE, verbose = FALSE, conserve.memory = FALSE, vst.flavor = "v2")
# Write to disk
saveRDS(merged_seurat_data, file = paste0(output_path, file_out_name, ".rds"), compress = FALSE)

