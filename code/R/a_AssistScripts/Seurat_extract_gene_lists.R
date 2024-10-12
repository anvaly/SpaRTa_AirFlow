# Author: Anna Lyubetskaya. Date: 20-08-06

# Extract gene lists from representative experiments for reference
# This is an assist script


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/R/Utils/utils_10X_in_out.R")


## PARAMETERS ----


# Run name
run_name <- "human_FFPE-probes"


## PATHS ----


# Input folder
input_path <- "XXXX"

# Output folder
output_path <- "XXXX"

# Seurat RDS files are tagged as follows
data_regex_files <- ".rds"


## INGEST DATA ----


# Find all files following the regex
file_list <- dir(input_path, pattern=paste0(data_regex_files), full.names=TRUE)

# Ingest a set of RDS Seurat objects
seurat_list <- read_rds_list_simple_my(file_list)


## FIND GENE LISTS ----


gene_list <- list()

for(name in names(seurat_list)){
  
  gene_list[[paste(name, "Spatial_counts")]] <- sort(rownames(seurat_list[[name]]@assays$Spatial@counts))
  gene_list[[paste(name, "Spatial_data")]] <- sort(rownames(seurat_list[[name]]@assays$Spatial@data))
  gene_list[[paste(name, "Spatial_scale.data")]] <- sort(rownames(seurat_list[[name]]@assays$Spatial@scale.data))
  
  gene_list[[paste(name, "SCT_counts")]] <- sort(rownames(seurat_list[[name]]@assays$SCT@counts))
  gene_list[[paste(name, "SCT_data")]] <- sort(rownames(seurat_list[[name]]@assays$SCT@data))
  gene_list[[paste(name, "SCT_scale.data")]] <- sort(rownames(seurat_list[[name]]@assays$SCT@scale.data))

}

# Reduce gene list
gene_unique_list <- unique(gene_list)
