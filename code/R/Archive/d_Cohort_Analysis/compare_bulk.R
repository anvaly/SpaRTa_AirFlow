# Author: Anna Lyubetskaya. Date: 20-08-06


##_ SETUP ENVIRONMENT _##


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_archs4.R")
source("code/utils/utils_pca.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_clinical_heritage.R")
source("code/utils/utils_rna_diff_expr.R")

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_matrix.R")


## PARAMETERS ----


# The cohort of interest regex ID
cohort_name <- "Syngeneic_FF"
project_id <- "P-20200928-0001"


## PATHS ----


# Seurat input folder
input_folder <- "C:/Users/lyubetsa/Documents/Projects_NewTech/10X_ST_Reports/"

# Seurat RDS files are tagged as follows
data_regex_files <- "_all.Sobj.rds"

# Output folder
output_folder <- "C:/Users/lyubetsa/Documents/Projects_NewTech/10X_ST_BulkComparison/"


## INGEST SEURAT DATA ----


# Find all files following the regex
file_list <- dir(input_folder, pattern=paste0(cohort_name, ".*", data_regex_files), full.names=TRUE)

# Ingest clean Seurat data
seurat_list <- list()  
for(rds_file in file_list){
  seurat_list[[gsub(paste0("^.+/|", data_regex_files), "", rds_file)]] <- as.list(read_seurat_rds_my(rds_file, do_pseudobulk=TRUE))
}

# Pseudo bulk Seurat tibble
seurat_data_df <- pseudo_bulk_vectors_combine_my(seurat_list, gsub(paste0("^.+/|", data_regex_files), "", file_list))

colnames(seurat_data_df) <- gsub(paste0(project_id, "_", cohort_name), "ST", colnames(seurat_data_df))


## JOIN, NORMALIZE, FILTER ESET ----


# Merge all 3 datasets
data_df <- seurat_data_df %>%
  dplyr::inner_join(bms_data_df, by="Symbol") %>%
  dplyr::inner_join(archs4_data_df, by="Symbol") %>%
  dplyr::mutate_if(is.double, as.integer)

# Calculate CPM values
data_full_eset <- Biobase::ExpressionSet(as.matrix(data_df %>% 
                                                     tibble::column_to_rownames("Symbol")))
data_norm_eset <- edgeR_Normalize(data_full_eset)

data_norm_matrix <- Biobase::exprs(data_norm_eset)
genes_represented <- names(which(rowSums(data_norm_matrix >= 1) >= 10))

# Normalized data wide tibble
data_norm_wide_df <- tibble::as_tibble(Biobase::exprs(data_norm_eset), rownames="Symbol") %>%
  dplyr::filter(Symbol %in% genes_represented)

# Normalized data long tibble
data_norm_df <- data_norm_wide_df %>%
  df_wide2long_my(key="Sample_ID", val="Expression")

# Meta data
meta_list <- sapply(unique(data_norm_df$Sample_ID), function(x) strsplit(x, "_"))
meta_df <- tibble::tibble(Sample_ID = names(meta_list),
                          Dataset = sapply(unname(meta_list), function(x) x[1]),
                          CellLine = sapply(unname(meta_list), function(x) x[2]))


## PCA ----


params <- list(sample_value = "Sample_ID", sample_filter = data_norm_df$Sample_ID,
               feature_value = "Symbol", feature_filter = data_norm_df$Symbol,
               cell_value = "Expression", color_by = c("Dataset", "CellLine"))

filename <- paste0(output_folder, "pca_cell_line")
# Perform first round of PCA
pc1_df <- pca_from_long_tibble_my(data_norm_df, meta_df, params, filename=filename, table2file=FALSE, extra_plot=TRUE)

# Sample extreme outlier sample using first PCA
sample_filt_list <- pc1_df %>% 
  dplyr::filter(.fittedPC1 < 0 & .fittedPC2 < 50) %>%
  dplyr::pull(Sample_ID)

params <- list(sample_value = "Sample_ID", sample_filter = sample_filt_list,
               feature_value = "Symbol", feature_filter = data_norm_df$Symbol,
               cell_value = "Expression", color_by = c("Dataset", "CellLine"))

filename <- paste0(output_folder, "pca_cell_line_filt")
pc2_df <- pca_from_long_tibble_my(data_norm_df, meta_df, params, filename=filename, table2file=FALSE, extra_plot=TRUE)


## BATCH CORRECTION ----


combat_input_matrix <- data_norm_wide_df %>% 
  tibble::column_to_rownames("Symbol") %>%
  dplyr::select(dplyr::all_of(sample_filt_list)) %>%
  as.matrix()

meta_filt_df <- meta_df %>%
  dplyr::filter(Sample_ID %in% sample_filt_list) %>%
  dplyr::arrange(Sample_ID)

colnames(combat_input_matrix) == meta_filt_df$Sample_ID

# Combat with parametric adjustment to correct for batches
combat_matrix = sva::ComBat(combat_input_matrix, batch=meta_filt_df$Dataset, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

# Normalized and batch corrected data long tibble
data_batch_df <- tibble::as_tibble(combat_matrix, rownames="Symbol") %>%
  df_wide2long_my(key="Sample_ID", val="Expression")

params <- list(sample_value = "Sample_ID", sample_filter = data_norm_df$Sample_ID,
               feature_value = "Symbol", feature_filter = data_norm_df$Symbol,
               cell_value = "Expression", color_by = c("Dataset", "CellLine"))

filename <- paste0(output_folder, "pca_cell_line_batch_filt")
pc2_df <- pca_from_long_tibble_my(data_batch_df, meta_filt_df, params, filename=filename, table2file=FALSE, extra_plot=TRUE)
