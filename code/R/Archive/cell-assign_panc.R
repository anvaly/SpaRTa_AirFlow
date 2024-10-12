# Author: Anna Lyubetskaya. Date: 20-10-06
# CellAssign: https://irrationone.github.io/cellassign/index.html


## ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

library(Seurat)
source("code/utils/utils_tibble.R")

# Setup for CellAssign
# devtools::install_github("rstudio/keras")
# library(keras)
# install_keras(method = "conda")
# install_keras(tensorflow = "gpu")
# tensorflow::install_tensorflow()
# tensorflow::install_tensorflow(extra_packages='tensorflow-probability', version = "2.1.0")


## PARAMETERS ----


# Sample ID
#sample_id <- "P-20200612-0001_ST_FFPE_HumanPanc_A"
sample_id <- "P-20200625-0002_ST_FF_HumanPanc_A"


## 10X PATHS ----


# Project path
global_path <- "C:/Users/lyubetsa/Documents/Projects_NewTech/"

# Path to the processed RDS data
input_file <- paste0(global_path, "10X_ST_Reports/", sample_id, "_all.Sobj.rds")

# Output folder
output_path <- paste0(global_path, "10X_ST_Annotation/", sample_id, "/")


## ANNOTATION PATHS ----


# Location of references
input_ref_path <- paste0(global_path, "References/Panc/")

# A list of reference files
ref_file_list <- c("PDAC_atlas_Peng_2019/41422_2019_195_MOESM9_ESM.xls",  # Peng 2019 PDAC cell populations
                   "PDAC_atlas_Peng_2019/41422_2019_195_MOESM11_ESM.xls",  # Peng 2019 PDAC cell subpopulations
                   "Panc_atlas_Muraro_2016/mmc3.xlsx",  # Muraro 2016 panc cell DEA genes
                   "Panc_atlas_Muraro_2016/mmc4.xlsx",  # Muraro 2016 panc cell DEA TFs
                   "Panc_atlas_Muraro_2016/mmc7.xlsx"  # Muraro 2016 panc cell DEA cell surface genes
)

# Simplified dataset names
name_list <- c("Peng2019_groups", "Peng2019_subgroups", "Muraro2016_genes", "Muraro2016_tfs", "Muraro2016_cellsurf")
names(name_list) <- ref_file_list

# How many lines to skip when reading datasets
skip_list <- c(1, 1, 0, 0, 0)
names(skip_list) <- ref_file_list

# Column names for individual datasets before and after wrangling
group_list <- list(list(c("Gene", "p value", "avg_logFC", "Cell type"), c("Symbol", "pvalue", "logFC", "Cell_type")),
                   list(c("Gene", "p value", "avg_logFC", "Subpopulation"), c("Symbol", "pvalue", "logFC", "Cell_type")), 
                   list(c("...1", "padj", "log2FoldChange"), c("Symbol", "pvalue", "logFC", "Cell_type")),
                   list(c("...1", "padj", "log2FoldChange"), c("Symbol", "pvalue", "logFC", "Cell_type")),
                   list(c("...1", "padj", "log2FoldChange"), c("Symbol", "pvalue", "logFC", "Cell_type")))
names(group_list) <- ref_file_list


## INGEST DATA ----


# Load Seurat data
data_seurat <- readRDS(input_file)


df_list <- list()

# Loop through XLS files and their sheets to create a list of tibbles
for(f in names(skip_list)){
  input_file <- paste0(input_ref_path, f)
  sheet_list <- readxl::excel_sheets(input_file)
  
  for(s in sheet_list){
    
    # Read an XLS file
    df <- readxl::read_excel(input_file, sheet=s, skip=skip_list[[f]]) %>%
      dplyr::select(group_list[[f]][[1]])
    
    # If the cell type is captured in the sheet name, add it to the tibble
    if(s != "Sheet1"){
      df <- df %>%
        dplyr::mutate(Cell_type = s)
    }
    
    # Standardize column names
    colnames(df) <- c(group_list[[f]][[2]])
    
    # Make cell types characters and add the dataset name
    df <- df %>%
      dplyr::mutate(Cell_type = as.character(Cell_type),
                    Source = name_list[[f]])
    
    # Add a suffix to the cell type that is only a number
    if(f == "PDAC_atlas_Peng_2019/41422_2019_195_MOESM11_ESM.xls"){
      df <- df %>%
        dplyr::mutate(Cell_type = paste0("PDAC_", Cell_type))
    }
    
    df_list[[paste0(f, "_", s)]] <- df
    
  }
}


## WRANGLE DATA ----


# Create a joint tibble and filter for FC
data_df <- dplyr::bind_rows(df_list) %>%
  dplyr::filter(abs(logFC) >= 1) %>%
  dplyr::mutate(Status = 1)

# Update cell type names
data_df$Cell_type <- gsub("Acinar cell", "acinar", data_df$Cell_type)
data_df$Cell_type <- gsub("Ductal cell", "duct", data_df$Cell_type)
data_df$Cell_type <- gsub("Endocrine cell", "endocrine", data_df$Cell_type)
data_df$Cell_type <- gsub("Endothelial cell", "endothelial", data_df$Cell_type)
data_df$Cell_type <- gsub("Stellaate cell", "stellate", data_df$Cell_type)

# Count groups per gene
count_groups_per_gene <- data_df %>%
  dplyr::group_by(Symbol) %>%
  dplyr::summarise(Groups_per_gene = dplyr::n_distinct(Cell_type),
                   Sources_per_gene = dplyr::n_distinct(Source),
                   Group_List = paste(unique(Cell_type), collapse="; ")) %>%
  dplyr::arrange(desc(Groups_per_gene))

# Count genes per group
count_genes_per_group <- data_df %>%
  dplyr::group_by(Cell_type) %>%
  dplyr::summarise(Genes_per_group = dplyr::n_distinct(Symbol)) %>%
  dplyr::arrange(desc(Genes_per_group))


## CREATE MARKER MATRIX ----


# Select only marker - cell type data
data_prep_df <- data_df %>%
  dplyr::select(c("Symbol", "Cell_type", "Status")) %>%
  unique

# Create a wide marker cell type matrix
data_wide_matrix <- df_long2wide_my(data_df, rows="Symbol", cols="Cell_type", value="Status") %>%
  dplyr::mutate_all(list(~tidyr::replace_na(.,0))) %>%
  tibble::column_to_rownames("Symbol") %>%
  as.matrix()

# Marker heatmap
# pheatmap::pheatmap(as.matrix(data_wide_df))


## CREATE MARKER MATRIX ----


# Extract counts
expression_matrix <- as.matrix(Seurat::GetAssayData(data_seurat, assay="RNA"))

# Calculate sum factors
# https://rdrr.io/bioc/scran/man/computeSumFactors.html
sum_factors <- scran::calculateSumFactors(expression_matrix)

# Make sure that all biomarkers are present in 10X dataset
data_wide_matrix_filt <- data_wide_matrix[intersect(rownames(data_seurat), rownames(data_wide_matrix)),]
expression_matrix_filt <- t(expression_matrix[rownames(data_wide_matrix_filt),])

# Annotate data using CellAssign
cellassign_fit <- cellassign::cellassign(exprs_obj = expression_matrix_filt, 
                                         marker_gene_info = data_wide_matrix_filt, 
                                         s = sum_factors, 
                                         learning_rate = 1e-2, 
                                         shrinkage = TRUE,
                                         verbose = TRUE)
