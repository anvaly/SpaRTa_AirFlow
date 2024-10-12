# Author: Anna Lyubetskaya. Date: 20-11-27

# Create a normalized, filtered Seurat object from 10X single-cell or single-nuclei folder data

# For a list of samples from a meta data file, 
# 1. load data
# 2. upper case all gene symbols
# 3. augment meta data
# --- data from the input meta data file, with prefix "user." in the misc slot
# --- data from the 10X metrics file, with prefix "10X." in the misc slot
# 5. Calculate ribo / mito content
# 6. Claculate additional QC metrics
# 7. SCT transform
# 8. save to RDS as an object

# Input meta data file contain te following columns:
# - Sample_Name, and FullPath: Sample_Name is the RDS object output name, FullPath is the location of the 10X spatial report


## ENVIRONMENT ----


# Allow piping throughout the package
`%>%` <- magrittr::`%>%`

source("code/utils/utils_in_out.R")
source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_vis.R")


## PARAMETERS ----


input_file_name <- "cellbender_feature_bc_matrix_filtered.h5"
feature_min <- 500
mito_max <- 20
ribo_max <- 40


## PATHS ----


# Input location containing 10X folders
input_path <- "XXXX"

# Output location for Seurat objects
output_path <- "XXXX"

# Create output subfolders
output_folders <- create_output_subfolders_my(output_path, c("Seurat_object"))

# Input file with sample meta data
meta_file <- paste0(input_path, "meta_data_snRNA.txt")


## INGEST, WRANGLE, OUTPUT ----


# Read sample file
meta_df <- readr::read_delim(file=meta_file, delim="\t") %>%
  tidyr::drop_na()


for(i in 1:nrow(meta_df)){
  
  # Load the 10X single-cell folder
  filename_in <- meta_df[[i, "FullPath"]]
  sample_name <- meta_df[[i, "Sample_Name"]]
  
  # Output filename for the RDS object
  filename_out <- paste0(output_folders$Seurat_object, sample_name, "_all.rds")
  
  if(!file.exists(filename_out)){
    print(filename_in)
    
    
    ## Load data ----
    
    # Read a 10X single-cell folder
    data_seurat_init <- read_10X_H5_my(paste0(filename_in, input_file_name))
    
    # Remove cells with zero UMI counts
    barcodes_keep <- rownames(data_seurat_init@meta.data)[which(data_seurat_init@meta.data$nCount_RNA > 0)]
    data_seurat <- subset(data_seurat_init, cells = barcodes_keep)
    
    # Read a 10X defined metrics file
    qc_df <- read_10X_qc_metrics_my(filename_in)
    
    
    ## Update data with parameters ----
    
    # Add user-defined meta-data to a Seurat object misc slot tagged with "user." prefix
    for(v in colnames(meta_df)){
      data_seurat <- add_misc_to_seurat_object_my(data_seurat, field=paste0("user.", v), dict_value=meta_df[[i, v]])
    }
    
    # Add 10X default QC metrics to a Seurat object tagged with "10X." prefix
    for(col in 1:ncol(qc_df)){
      data_seurat <- add_misc_to_seurat_object_my(data_seurat, field=paste0("10X.", colnames(qc_df)[col]), dict_value=qc_df[[1, col]])
    }
    
    
    ## Calculate ribosomal and mitochondrial content ----
    
    
    # Calculate mitochondrial and ribosomal content
    data_seurat <- calculate_mt_ribo_my(data_seurat)
    
    # List of cells
    cell_list <- data_seurat@meta.data %>% 
      dplyr::filter(mito_percent <= mito_max & 
                      ribo_percent <= ribo_max &
                      nFeature_RNA >= feature_min) %>%
      rownames()
    
    cat(filename_in, "\n",
        "Cell number at ingestion =", ncol(data_seurat_init), "\n",
        "Cell number at pre-filtering =", ncol(data_seurat), "\n",
        "Cell number at filtering =", length(cell_list), "\n")
    
    # Remove cells with high MT or ribo content or low feature content
    data_seurat <- subset(data_seurat, cells = cell_list)
    
    
    ## SCTransform ----
    
    # Perform SCTransform normalization
    data_seurat <- Seurat::SCTransform(data_seurat, assay="RNA", vars.to.regress = c("mito_percent", "ribo_percent"),
                                       return.only.var.genes = FALSE, verbose = FALSE)
    
    
    ## Calculate mean QC metrics on final objects ----
    
    # List of meta data fields to parse
    field_list <- c("nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT", "mito_percent", "ribo_percent")
    
    # Add mean QC metrics to the misc slot
    for(f in field_list){
      val <- mean(data_seurat@meta.data[[f]])
      data_seurat <- add_misc_to_seurat_object_my(data_seurat, field=paste0("qc.mean.", f), dict_value=val)
    }
    
    
    ## Write data to RDS ----
    
    # Save the full dataset to RDS
    saveRDS(data_seurat, file = filename_out)
  }
}
