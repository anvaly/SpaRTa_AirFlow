# Author: Anna Lyubetskaya. Date: 20-02-04

# This script:
# - Calculates pseudo-bulk expression profiles for a set of ST samples
# - Calculates pairwise correlations between pseudo-bulk expression vectors
# - Fits a linear regression to each pair of comparisons and identifies outliers as mean+3sd of residuals abs values
# - Output log10(UMI+1) psedo-bulk expression vectors, a correlation heatmap, and pairwise scatter plots with a linear regression fit and outlier calls


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_stats.R")
source("code/utils/utils_rna_diff_expr.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_signatures.R")
source("code/utils/utils_in_out.R")
source("code/utils/utils_ggplot.R")
source("code/utils/utils_specialized_plots.R")

source("code/R/Utils/utils_10X_matrix.R")


## PARAMETERS ----


# Seurat RDS files are tagged as follows
cohort_name <- "PDAC84_path12"
cohort_regex <- "PDAC"

# Exclude certain samples from analysis
sample_exclude <- NULL
# A list of all cohorts
cohort_col <- ""

# Make sample scatter plots
do_scatter_plots <- TRUE
# Number of pairwise sample scatter plots to generate
plot_num <- 10

## A way to abstract data from individual samples to an average of another user-defined grouping
# Make mean of groups scatter plots
do_mean_plots <- FALSE
# Define mean sample pairs to compare
cohort_pair_list <- NULL

# A list of genes to plot for pseudobulk spot checking
gene_select <- c("FZD5", "CEACAM5", "MAP3K14", "CD8A")  # NULL


## PATHS ----


# Input folder
input_paths <- c("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_name, "/")

dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Find all files following the regex
file_list <- unname(unlist(sapply(input_paths, function(p) dir(p, pattern=paste0(cohort_regex, ".*.rds"), full.names=TRUE))))

# Ingest RDS objects with Seurat data and transform them into a pseudo-count vector
pseudo_count_list <- list()
occupancy_count_list <- list()
meta_data_list <- list()
for(rds_file in file_list){
  # Load the dataset
  data_seurat <- readRDS(rds_file)
  
  if(!data_seurat@misc[["user.Sample_Name"]] %in% sample_exclude){
    # Sample name
    if(cohort_col != ""){
      sample_name <- paste(data_seurat@misc[[cohort_col]], data_seurat@misc[["user.Sample_Name"]])
    } else{
      sample_name <- data_seurat@misc[["user.Sample_Name"]]
    }
    
    # Extract meta data information
    meta_data_list[[sample_name]] <- data_seurat@misc
    
    # Calculate pseudo bulk counts for a 10X sample
    pseudo_count_list[[sample_name]] <- pseudo_bulk_counts_my(data_seurat, log_transform=TRUE, assay="Spatial", slot="data")
    occupancy_count_list[[sample_name]] <- occupancy_counts_my(data_seurat, log_transform=FALSE, assay="Spatial", slot="data")
  }
}


## WRANGLE META DATA ----


# Create a meta data tibble
meta_df <- list2tibble_my(meta_data_list)

# List of sample groups
cohort_list <- unique(meta_df[[cohort_col]])

# All combinations of cohort pairs
if(length(cohort_list) >= 2 && is.null(cohort_pair_list)){
  cohort_pair_list <- combn(cohort_list, 2)
  cohort_pair_list <- cohort_pair_list[, unique(round(runif(plot_num, 1, ncol(cohort_pair_list))))]
}


## WRANGLE COUNTS DATA ----


# Combine a named list of lists into an expression wide tibble
data_wide_df <- pseudo_bulk_vectors_combine_my(pseudo_count_list, names(pseudo_count_list)) %>%
  dplyr::mutate_at(vars(dplyr::all_of(names(pseudo_count_list))), as.numeric)

# Find gene abundances across all samples
data_df <- data_wide_df %>%
  df_wide2long_my(key="Sample", val="Value") %>%
  dplyr::group_by(Symbol) %>%
  dplyr::summarise(Sum = sum(Value))

hist(data_df$Sum)

# Find genes that have no counts
genes_remove <- data_df %>%
  dplyr::filter(Sum > 0) %>%
  dplyr::pull(Symbol)

# Filter expression matrix for genes that have counts
# Remove mitochondrial and ribosomal genes from the count matrix
data_wide_df <- data_wide_df %>%
  dplyr::filter(!Symbol %in% identify_mt_ribo_genes_my(data_wide_df$Symbol)) %>%
  dplyr::filter(Symbol %in% genes_remove)

# Write pseudo-bullk data to file
filename <- paste0(output_path, "matrix_pseudo_sct_log10.txt")
readr::write_delim(data_wide_df, filename, delim = "\t")

# Combine a named list of lists into a gene occupancy wide tibble
# Highlight gene symbols appearing in our signature DB
occupancy_wide_df <- pseudo_bulk_vectors_combine_my(occupancy_count_list, names(occupancy_count_list)) %>%
  dplyr::mutate_at(vars(dplyr::all_of(names(occupancy_count_list))), as.numeric)


## CALCULATE PAIRWISE CORRELATION ----


filename <- paste0(output_path, "corr_matrix_reps")
# Calculate pairwise correlations between every sample represented as a pseudo bulk vector
corr_wide_df <- correlation_calculate_plot_my(data_wide_df, 
                                              row_col="Symbol", output_file=filename)

# Wide to long correlation tibble
corr_df <- corr_wide_df %>%
  dplyr::rename(Sample1 = term) %>%
  df_wide2long_my(key="Sample2", val="R") %>%
  tidyr::drop_na()

# Create a unique sample pair name
#corr_df["SamplePair"] <- sapply(1:nrow(corr_df), function(x) paste(sort(c(corr_df[[x, "Sample1"]], corr_df[[x, "Sample2"]])), collapse=":"))
corr_df["SamplePair"] <- sapply(1:nrow(corr_df), function(x) paste(sort(c(corr_df[[x, "Sample1"]])), collapse=":"))

# Remove replicate correlation plots
corr_df <- corr_df %>%
  dplyr::group_by(SamplePair) %>%
  dplyr::summarise(R = unique(R)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Sample1 = gsub(":.+$", "", SamplePair),
                Sample2 = gsub("^.+:", "", SamplePair)) %>%
  dplyr::select(-SamplePair)

cat("Correlation:", min(corr_df$R), max(corr_df$R), "\n")


## PLOT PAIRWISE CORRELATION ----


# Plot a correlation heatmap
midpoint <- round(1 - (1 - min(corr_df$R))/2, 2)
correlation_plot_my(corr_wide_df, scale = c(round(min(corr_df$R), 2), 1, midpoint), cols=c("white", "darkblue", "lightblue"), 
                    filename=filename, rowname_col="term")


## SCATTER PLOT BETWEEN TWO MOST DISSIMILAR SAMPLES IN EACH SET ----


# For each cohort, find the most dissimilar samples and compare their expression via a scatter plot
if(do_scatter_plots == TRUE){
  
  for(i in 1:nrow(corr_df)){
    
    # Select data for the comparison
    data_loc_df <- data_wide_df %>%
      dplyr::select(all_of(c("Symbol", corr_df[[i, "Sample1"]], corr_df[[i, "Sample2"]])))
    
    # Find the correlation coefficient
    corr_loc_df <- corr_df %>% 
      dplyr::filter(Sample1 == corr_df[[i, "Sample1"]] & Sample2 %in% corr_df[[i, "Sample2"]]) %>%
      dplyr::arrange(R, Sample1, Sample2)
    
    # Calculate linear regression and residuals, identify outliers
    data_loc_df <- lm_residual_outliers(data_loc_df, corr_df[[i, "Sample1"]], corr_df[[i, "Sample2"]])

    sample1 <- gsub("^.* ", "", corr_df[[i, "Sample1"]])
    sample2 <- gsub("^.* ", "", corr_df[[i, "Sample2"]])
    
    # Create a scatter plot of 2 most dissimilar samples within the cohort
    filename <- paste0(output_path, "scatter_pair_", sample1, "_", sample2)
    create_scatter_plot_my(data_loc_df, x_label=corr_df[[i, "Sample1"]], y_label=corr_df[[i, "Sample2"]], 
                           fill_label="Is`Outlier`", filename=filename, 
                           labels=c(sample1, sample2, paste("R =", round(corr_df[[i, "R"]], 3))),
                           size=0.5, do_fit=TRUE, stroke=0)
    
    filename <- paste0(output_path, "scatter_pair_", sample1, "_", sample2, ".txt")
    readr::write_delim(data_loc_df, filename, delim = "\t")
    
  }
  
}


## MEAN EXPRESSION FOR SAMPLE GROUPS ----


if(!is.null(cohort_pair_list) && do_mean_plots == TRUE){
  
  
  ## CALCULATE MEAN DATASET COUNTS ----
  
  
  data_mean_list <- list()
  
  # For each cohort, find all corresponding experiments and average the expression across all samples
  for(cohort in cohort_list){
    # Find all samples that belong to the cohort
    col_list <- colnames(data_wide_df)[grep(cohort, colnames(data_wide_df))]
    
    if(length(col_list) > 1){
      # Calculate mean expression for every gene across all samples in the cohort
      data_mean_list[[cohort]] <- as.list(round(rowMeans(data_wide_df %>%
                                                           dplyr::select(all_of(c("Symbol", col_list))) %>% 
                                                           tibble::column_to_rownames("Symbol")), 3))
    } else if(length(col_list) == 1){
      data_mean_list[[cohort]] <- as.list(data_wide_df[[col_list]])
      names(data_mean_list[[cohort]]) <- data_wide_df$Symbol
    }
  }
  
  # Generate a tibble of mean expression of each dataset
  data_mean_df <- list2tibble_my(data_mean_list, rownames="Symbol", transpose=TRUE) %>%
    dplyr::mutate_at(vars(dplyr::all_of(names(data_mean_list))), as.numeric)
  
  filename <- paste0(output_path, "matrix_mean_", cohort_name, ".txt")
  readr::write_delim(data_mean_df, filename, delim = "\t")
  
  
  ## PLOT MEAN EXPRESSION FOR SAMPLE GROUPS ----
  
  
  # Create a scatter plot of mean expression values for input datasets
  for(j in 1:ncol(cohort_pair_list)){
    
    # Select a pair of experiments from a list of pairs
    value1 <- cohort_pair_list[1, j]
    value2 <- cohort_pair_list[2, j]
    
    # Select expression data for a pair of experiments
    # data_delta_df <- data_mean_df %>%
    #   dplyr::select(dplyr::all_of(c("Symbol", value1, value2, "InSignature")))
    data_delta_df <- data_mean_df %>%
      dplyr::select(dplyr::all_of(c("Symbol", value1, value2)))
    
    # Calculate the delta and other stats between conditions
    data_delta_df[["Delta"]] <- round(data_delta_df[[value1]] - data_delta_df[[value2]], 3)
    data_delta_df[["DeltaAbs"]] <- abs(data_delta_df[["Delta"]])
    data_delta_df[["Mean"]] <- round((data_delta_df[[value1]] + data_delta_df[[value2]]) / 2, 3)
    threshold <- sd(data_delta_df[["Delta"]]) * 3
    
    # Define strong outliers
    data_delta_df <- data_delta_df %>%
      dplyr::mutate(Outlier = abs(Delta) >= threshold)
    
    # Scatter plot of mean values between different protocols applied to the same cohort
    filename <- paste0(output_path, "scatter_mean_", value1, "_", value2)
    create_scatter_plot_my(data_delta_df, x_label=value1, y_label=value2, 
                           fill_label="Outlier", filename=filename, do_fit=TRUE, size=0.5, stroke=0,
                           labels=c(value1, value2, paste("R =", round(cor(data_delta_df[[value1]], data_delta_df[[value2]]), 3))))
    
    # Write information about expression differences between a pair of experiments to file
    data_delta_df <- data_delta_df %>%
      dplyr::filter(Mean >= 2) %>%
      dplyr::arrange(desc(DeltaAbs)) %>%
      dplyr::mutate(IsRiboMito = grepl("^MT-|^RPL|^RPS", Symbol))
    
    filename <- paste0(output_path, "scatter_mean_", value1, "_", value2, ".txt")
    readr::write_delim(data_delta_df, filename, delim = "\t")
    
  }
}
