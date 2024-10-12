# Author: Anna Lyubetskaya. Date: 23-07-24
# Split cohort by pathology


## ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")
source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_in_out.R")


## PARAMETERS ----


# Name of the analysis to use in folder/file names
run_name <- "rest"

# Path to processed Seurat data
cohort_name <- "PDAC108_path14_merge"
sample_exclude <- NULL


## COMPARTMENT


# Name of the pathology field to filter by groups
pathology_select <- paste0("Pathology.", c("NonEpi"), ".percent")
pathology_name <- "Pathology.RestAll.percent"

# Feature # in a spot threshold
feature_threshold <- 300

# Subset and write out a Seurat object
subset_def <- "Pathology.RestAll.percent"
subset_threshold <- -100


## PATHS ----


# A list of Seurat objects which spots should be explicitly excluded
seurat_spot_exclude <- c("XXXX",
                         "XXXX")

# Location of pre-processed data
input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, run_name, "_", cohort_name, "/")

# Path to which to write a subset Seurat object
output_path_subset <- "XXXX"

# Log information
log_file <- paste0(output_path, "log.txt")

# Create output folders
dir.create(output_path_subset, showWarnings = FALSE)
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Ingest Seurat objects that contain spots that should be excluded from analyses
spot_exclude_list <- c()
if(!is.null(seurat_spot_exclude)){
  
  # Find all files following the regex
  file_list <- unname(unlist(sapply(seurat_spot_exclude, function(p) dir(p, pattern=".*.rds", full.names=TRUE))))
  
  # Ingest a set of RDS Seurat objects
  spot_exclude_list <- sapply(read_rds_list_simple_my(file_list), function(x) colnames(x))
 
  gc() 
}

# Seurat data
data_seurat <- readRDS(input_path)

# Remove spots specified by the user
if(length(spot_exclude_list) >= 1){
  spot_exclude <- c(unlist(unname(spot_exclude_list)))
  
  cat(ncol(data_seurat), length(spot_exclude), length(unique(spot_exclude)), "\n")
  cat(ncol(data_seurat) - length(unique(spot_exclude)), "\n")
  
  barcode_list <- setdiff(colnames(data_seurat), spot_exclude)
  cat(length(barcode_list))
  
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# Remove samples if necessary
if(!is.null(sample_exclude)){
  barcode_list <- rownames(data_seurat@meta.data[which(!data_seurat@meta.data[["user.Sample_Name"]] %in% sample_exclude),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# Remove spots if necessary
if(!is.null(feature_threshold)){
  barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data$nFeature_Spatial >= feature_threshold),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}


## WRANGLE META DATA ----


# Extract coordinate meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

# Create a joint pathology label
meta_df[[pathology_name]] <- matrixStats::rowMaxs(as.matrix(meta_df[pathology_select]))

# Add user.Sample_Name to meta data if absent: happens in stand-alone samples
if(!"user.Sample_Name" %in% colnames(meta_df)){
  meta_df["user.Sample_Name"] = data_seurat@misc$user.Sample_Name
}


## SUBSET OBJECT ----


# Subset and write out a Seurat object
table(meta_df[subset_def])


# Find samples and corresponding barcodes to subset
subset_df <- meta_df %>%
  dplyr::filter(!!rlang::sym(subset_def) >= subset_threshold) %>% 
  dplyr::select(Coordinate, user.Sample_Name)


# Count spots selected at the current thresholds
print(nrow(subset_df))
print(nrow(meta_df))
print(round(nrow(subset_df) / nrow(meta_df) * 100))

sort(table(subset_df$user.Sample_Name))


# Create objects for each sample
for(s in unique(subset_df$user.Sample_Name)){
  
  # Barcodes for the select samples
  cell_list <- subset_df %>%
    dplyr::filter(user.Sample_Name == s) %>%
    dplyr::pull(Coordinate)
  
  filename_rds <- paste0(output_path_subset, cohort_name, "_", s, "_", run_name, ".rds")
  
  if(!file.exists(filename_rds)){
    data_subset_seurat <- subset(data_seurat, cells=cell_list)
    
    # Only leave the relevant image
    data_subset_seurat@images <- data_subset_seurat@images[unique(data_subset_seurat@meta.data$user.Sample_Name)]
    
    saveRDS(data_subset_seurat, filename_rds)
  }
  
}
