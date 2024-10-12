# Author: Anna Lyubetskaya. Date: 22-02-16

# This script:
# - Calculates pseudo-bulk expression profiles for a set of ST samples and normalize them
# - Define two groups of samples by selecting a meta data column and specifying which values from that column goes into each group
# - Calculate mean edgeR normalized gene expression values for each group
# - Calculate a lm fit and residuals and find any gene that is more than 3SD deviations away from the mean


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_pca.R")
source("code/utils/utils_stats.R")
source("code/utils/utils_ggplot.R")
source("code/utils/utils_rna_diff_expr.R")

source("code/R/Utils/utils_10X_matrix.R")


## PARAMETERS ----


# Run name to name the output folder
run_name <- "polyA_v_probes_norm"

# Seurat RDS files are tagged as follows
cohort_name <- "PDAC"
# Samples to include
cohort_regex <- "PDAC"

# Normalize data
do_norm <- TRUE

# Meta data column and values within it to use to define groups
meta_data_col <- "user.Protocol"

group1 <- c("FF-polyA", "FFPE-polyA", "FF-polyA-woCS")
group2 <- c("FFPE-probes-v1")

group1_name <- "PDAC_polyA"
group2_name <- "PDAC_probes"


## PATHS ----


# Input folder
input_paths <- c("XXXX",
                 "XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_name, "_", run_name, "/")

dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Find all files following the regex
file_list <- unname(unlist(sapply(input_paths, function(p) dir(p, pattern=paste0(cohort_regex, ".*.rds"), full.names=TRUE))))


# Ingest RDS objects with Seurat data and tranform them into a pseudo-count vector
pseudo_count_list <- list()
meta_data_list <- list()
for(rds_file in file_list){
  
  # Load the dataset
  data_seurat <- readRDS(rds_file)
  
  # Sample name
  sample_name <- data_seurat@misc[["user.Sample_Name"]]
  
  # Extract meta data information
  meta_data_list[[sample_name]] <- data_seurat@misc
  
  # Calculate pseudo bulk counts for a 10X sample
  pseudo_count_list[[sample_name]] <- pseudo_bulk_counts_my(data_seurat, log_transform=FALSE, assay="Spatial", slot="data")
  
}


## WRANGLE DATA ----


# Create a meta data tibble
meta_df <- list2tibble_my(meta_data_list)
meta_df$user.Sample_Name <- names(meta_data_list)

# Check the column for grouping
unique(meta_df[[meta_data_col]])
group1 %in% unique(meta_df[[meta_data_col]])
group2 %in% unique(meta_df[[meta_data_col]])

# Adjust meta data if 
# Combine a named list of lists into an expression wide tibble
data_wide_df <- pseudo_bulk_vectors_combine_my(pseudo_count_list, names(pseudo_count_list)) %>%
  dplyr::mutate_at(vars(dplyr::all_of(names(pseudo_count_list))), as.numeric)


## CREATE AND WRANGLE ESET ----


# Create an annotated data frame object of gene annotations
featureData <- new("AnnotatedDataFrame", data = data.frame(Symbol = data_wide_df$Symbol,
                                                           row.names = data_wide_df$Symbol))

# Create an annotated data frame object of sample phenotypes
phenoData <- new("AnnotatedDataFrame", data = meta_df %>% 
                   tibble::column_to_rownames("user.Sample_Name") %>%
                   as.data.frame())

# Create an eSet object
data_eset <- Biobase::ExpressionSet(assayData = as.matrix(data_wide_df %>% 
                                                            tibble::column_to_rownames("Symbol")),
                                    phenoData = phenoData,
                                    featureData = featureData)


# Calculate CPM values
data_norm_eset <- edgeR_normalize_my(data_eset)

# Extract the raw or CPM matrix from the eSet
if(do_norm == TRUE){
  data_norm_matrix <- Biobase::exprs(data_norm_eset)
} else{
  data_norm_matrix <- Biobase::exprs(data_eset)
}

dim(data_norm_matrix)
table(colnames(data_norm_matrix) == meta_df$user.Sample_Name)

# Find samples for each group
group1_indx <- which(meta_df[[meta_data_col]] %in% group1)
group2_indx <- which(meta_df[[meta_data_col]] %in% group2)

# Find group means
group1_mean <- rowMeans(data_norm_matrix[, group1_indx])
group2_mean <- rowMeans(data_norm_matrix[, group2_indx])

# Create a tibble of expression values
expr_df <- tibble::tibble(Symbol = names(group1_mean),
                          Group1 = unname(group1_mean),
                          Group2 = unname(group2_mean))

colnames(expr_df) <- c("Symbol", group1_name, group2_name)


## VISUALIZE AND FIND OUTLIERS ----


# Find the correlation coefficient
corr_coef <- round(cor(group1_mean, group2_mean), 3)

# Calculate linear regression and residuals, identify outliers
expr_df <- lm_residual_outliers(expr_df, group1_name, group2_name, sd_num=3)

# Write data to file
filename <- paste0(output_path, "scatter_pair_", group1_name, "_", group2_name, ".txt")
readr::write_delim(expr_df, filename, delim = "\t")

# Write data to file
filename <- paste0(output_path, "outliers_", group2_name, ".txt")
readr::write_delim(expr_df %>%
                     dplyr::filter(IsOutlier == paste0(group2_name, "_Outlier")) %>%
                     dplyr::select(Symbol) %>%
                     dplyr::arrange(Symbol), 
                   filename, delim = "\t")

# Create a scatter plot of 2 most dissimilar samples within the cohort
filename <- paste0(output_path, "scatter_pair_", group1_name, "_", group2_name)
create_scatter_plot_my(expr_df, x_label=group1_name, y_label=group2_name, 
                       fill_label="IsOutlier", filename=filename, 
                       labels=c(group1_name, group2_name, paste("R =", corr_coef)),
                       size=1, do_fit=TRUE, stroke=0)


## PCA ----


# Create a long tibble of expression CPMs
data_df <- tibble::as_tibble(data_norm_matrix, rownames="Symbol") %>%
  df_wide2long_my(key="user.Sample_Name", val="CPM")

# Identify genes present in all objects and that are not outliers
genes_present <- Reduce(intersect, sapply(pseudo_count_list, function(x) names(x)))
genes_represented <- names(which(rowSums(data_norm_matrix >= 1) >= 10))
genes_fit <- expr_df %>% 
  dplyr::filter(IsOutlier == "Fit") %>% 
  dplyr::pull(Symbol)
genes_pca <- intersect(genes_present, intersect(genes_represented, genes_fit))

# PCA parameters
params <- list(sample_value = "user.Sample_Name", 
               sample_filter = unique(data_df$user.Sample_Name),
               feature_value = "Symbol", 
               feature_filter = genes_pca,
               cell_value = "CPM", 
               color_by = c(meta_data_col))

# Perform PCA on the whole cohort
filename <- paste0(output_path, "pca_by_group")
pc_df <- pca_from_long_tibble_my(data_df, meta_df, params, filename=filename, table2file=FALSE, extra_plot=TRUE)
