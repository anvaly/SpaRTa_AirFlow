# Author: Anna Lyubetskaya. Date: 20-02-04

# This script:
# - Calculates pseudo-bulk expression profiles for a set of ST samples
# - Calculates pairwise correlations between pseudo-bulk expression vectors
# - Fits a linear regression to each pair of comparisons and identifies outliers as mean+3sd of residuals abs values
# - Output log2(UMI+1) psedo-bulk expression vectors, a correlation heatmap, and pairwise scatter plots with a linear regression fit and outlier calls


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

library(Seurat)


## PARAMETERS ----


# Seurat RDS files are tagged as follows
cohort_name <- "PDAC_FERMT"
cohort_regex <- ""

# Break data by pathology
break_path_classes <- TRUE

# Make sample scatter plots
do_scatter_plots <- FALSE
# Number of pairwise sample scatter plots to generate
plot_num <- 10

# A list of genes to plot for pseudobulk spot checking
gene_select <- c("FERMT1", "FERMT2", "FHL1")

# Meta data to sum across columns and add to sample-level meta data
meta_select <- paste0("Pathology.", c("IntestineAdj", "MuscleAdj", "NormalAdj",  "LymphNode", "TLSAggregate", "TLSImmature", "TLSMature", 
                                      "Adipose", "Vessel", "Muscle", "Nerve", "NonEpi", "BenignEpi", "Blood", "ExoEndo", "LuminalNec", "Tumor"), ".percent")


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

# Ingest RDS objects with Seurat data and tranform them into a pseudo-count vector
pseudo_count_list <- list()
meta_data_list <- list()
for(rds_file in file_list){
  
  # Load the dataset
  data_seurat <- readRDS(rds_file)
  
  # Sample name
  sample_name <- data_seurat@misc[["user.Sample_Name"]]
  
  
  # Find major pathology classes in the sample
  if(break_path_classes == TRUE){
    path_group_counts <- table(data_seurat@meta.data$Pathology.Group)
    path_groups <- names(path_group_counts[which(path_group_counts >= 100)])
  } else{
    path_groups <- c("skip")
  }
  
  
  # Calculate pseudo bulk counts for a 10X sample
  for(p in path_groups){
    
    if(p == "skip"){
      
      # Sum spot level meta data to sample level meta data
      meta_spot_list <- Matrix::colSums(data_seurat@meta.data[intersect(colnames(data_seurat@meta.data), meta_select)]) / nrow(data_seurat@meta.data)

      # Extract meta data information
      meta_data_list[[sample_name]] <- c(data_seurat@misc, as.list(meta_spot_list))
      
      pseudo_count_list[[sample_name]] <- pseudo_bulk_counts_my(data_seurat, log_transform=TRUE, assay="Spatial", slot="data")
      
    } else{
      
      # Select spots corresponding to a specific pathology compartment
      data_seurat_subset <- subset(data_seurat, cells=rownames(data_seurat@meta.data)[which(data_seurat@meta.data$Pathology.Group == p)])
      
      # Sum spot level meta data to sample level meta data
      meta_spot_list <- Matrix::colSums(data_seurat_subset@meta.data[intersect(colnames(data_seurat_subset@meta.data), meta_select)]) / nrow(data_seurat_subset@meta.data)

      # Extract meta data information
      meta_data_list[[paste(sample_name, p)]] <- c(data_seurat_subset@misc, as.list(meta_spot_list))
      
      pseudo_count_list[[paste(sample_name, p)]] <- pseudo_bulk_counts_my(data_seurat_subset, log_transform=TRUE, assay="Spatial", slot="data")
      
    }
    
  }
}


## WRANGLE META DATA ----


# Create a meta data tibble
meta_df <- list2tibble_my(meta_data_list)


## WRANGLE COUNTS DATA ----


# Combine a named list of lists into an expression wide tibble
data_wide_df <- pseudo_bulk_vectors_combine_my(pseudo_count_list, names(pseudo_count_list)) %>%
  dplyr::mutate_at(vars(dplyr::all_of(names(pseudo_count_list))), as.numeric)

# Long tibble of gene pseudo-count expression
data_df <- data_wide_df %>%
  df_wide2long_my(key="Sample", val="Value")

# Sum gene abundances across all samples
data_sum_df <- data_df %>%
  dplyr::group_by(Symbol) %>%
  dplyr::summarise(Sum = sum(Value))

hist(data_sum_df$Sum)

# Find genes that have at least some counts
genes_keep <- data_sum_df %>%
  dplyr::filter(Sum > 0) %>%
  dplyr::pull(Symbol)

# Filter expression matrix for genes that have counts
# Remove mitochondrial and ribosomal genes from the count matrix
data_wide_df <- data_wide_df %>%
  dplyr::filter(!Symbol %in% identify_mt_ribo_genes_my(data_wide_df$Symbol)) %>%
  dplyr::filter(Symbol %in% genes_keep)

# Write pseudo-bullk data to file
filename <- paste0(output_path, "matrix_pseudo_sct_log2.txt")
readr::write_delim(data_wide_df, filename, delim = "\t")


## CALCULATE PAIRWISE CORRELATION ----


# Calculate pairwise correlations between every sample represented as a pseudo bulk vector
filename <- paste0(output_path, "corr_pairwise_matrix")
corr_wide_df <- correlation_calculate_plot_my(data_wide_df, row_col="Symbol", output_file=filename)

# Subset meta data to only relevant fields
meta_loc_df <- meta_df %>%
  dplyr::select(user.Sample_Name, user.Region_ID, user.Block_ID)

# Wide to long correlation tibble
corr_df <- corr_wide_df %>%
  dplyr::rename(Sample1 = term) %>%
  df_wide2long_my(key="Sample2", val="R") %>%
  tidyr::drop_na() %>%
  dplyr::arrange(R) %>%
  dplyr::mutate(term1 = gsub(" .+", "", Sample1)) %>%
  dplyr::mutate(term2 = gsub(" .+", "", Sample2)) %>%
  dplyr::inner_join(meta_loc_df, by=c("term1" = "user.Sample_Name")) %>%
  dplyr::rename(Group1 = user.Region_ID, Block1 = user.Block_ID) %>%
  dplyr::inner_join(meta_loc_df, by=c("term2" = "user.Sample_Name")) %>%
  dplyr::rename(Group2 = user.Region_ID, Block2 = user.Block_ID) %>%
  dplyr::mutate(SameGroup = Group1 == Group2,
               SameBlock = Block1 == Block2,
               R = round(R, 3))

# Save a long tibble of correlations too
filename <- paste0(output_path, "corr_pairwise_list.txt")
readr::write_delim(corr_df, filename, delim="\t")


## PLOT PAIRWISE CORRELATION ----


# Plot a correlation heatmap
filename <- paste0(output_path, "corr_pairwise")
midpoint <- round((max(corr_df$R) - min(corr_df$R))/2 + min(corr_df$R), 2)
correlation_plot_my(corr_wide_df, scale = c(round(min(corr_df$R), 2), 1, midpoint), 
                    cols=c("yellow", "navyblue", "lightblue"), 
                    filename=filename, rowname_col="term")

# Histogram of correlation plots
filename <- paste0(output_path, "hist_corr_coef")
create_hist_plot_my(corr_df, x_label="R", fill_label="SameGroup", 
                    intercept=c(round(min(corr_df$R), 2), midpoint, 1), binwidth=0.01, 
                    labels=c("Pairwise correlation coefficient, R", "Number of comparisons", 
                             "Pairwise pseudo-bulk comparison between samples"),
                    filename=filename)


## BOX PLOT OF SELECT GENES ----


# Create a box plot of gene signature levels by experiment group, fixed axes
filename <- paste0(output_path, "violin_genes")
create_violin_plot_my(data_df %>%
                        dplyr::filter(Symbol %in% gene_select), 
                      x_label="Symbol", y_label="Value", fill_label="Symbol", filename=filename,
                      labels=c("Symbol", "UMI sum + 1, log2", ""))

if(break_path_classes == TRUE){
  filename <- paste0(output_path, "violin_genes_by_path")
  create_violin_plot_my(data_df %>%
                          dplyr::filter(Symbol %in% gene_select) %>%
                          dplyr::mutate(Pathology = gsub("^.+ ", "", Sample)), 
                        x_label="Pathology", y_label="Value", fill_label="Pathology", 
                        facet_var=c("Symbol", "free_x"), filename=filename,
                        labels=c("Symbol", "UMI sum + 1, log2", ""))

  
  data_loc_df <- data_df %>% 
    dplyr::filter(grepl("Tumor", Sample)) %>%
    dplyr::filter(Symbol %in% gene_select[c(1,2)]) %>%
    dplyr::mutate(Sample = gsub(" .+", "", Sample)) %>%
    dplyr::left_join(meta_df %>%
                       dplyr::select(user.Sample_Name, user.Block_ID),
                     by=c("Sample" = "user.Sample_Name"))
  
  filename <- paste0(output_path, "box_genes_by_samples")
  p <- create_box_plot_my(data_loc_df, 
                          x_label="user.Block_ID", y_label="Value", fill_label="Symbol", 
                          filename=NULL,
                          labels=c("Patient", "UMI sum + 1, log2", ""), reorder_x=TRUE)
  
  write_plot2file_my(p, filename, num_row=1, num_col=length(gene_select))
  
  filename <- paste0(output_path, "scatter_genes_by_samples")
  create_scatter_plot_my(data_loc_df %>%
                           dplyr::select(Sample, Symbol, Value) %>%
                           df_long2wide_my(rows="Sample", cols="Symbol", value="Value") %>%
                           dplyr::inner_join(meta_df %>%
                                               dplyr::select(user.Sample_Name, Pathology.Tumor.percent),
                                             by=c("Sample" = "user.Sample_Name")) %>%
                           dplyr::mutate(`Tumor%` = factor(round(Pathology.Tumor.percent, -1), levels=seq(0,100,10))), 
                         x_label=gene_select[1], y_label=gene_select[2], fill_label="Tumor%", 
                         filename=filename, do_fit=TRUE, size=3, stroke=0.1)
  
}


## SCATTER PLOT BETWEEN SAMPLE PAIRS ----


# For each cohort, find the most dissimilar samples and compare their expression via a scatter plot
if(do_scatter_plots == TRUE){
  
  for(i in 1:plot_num){
    
    # Select data for the comparison
    data_loc_df <- data_wide_df %>%
      dplyr::select(all_of(c("Symbol", corr_df[[i, "Sample1"]], corr_df[[i, "Sample2"]])))
    
    # Find the correlation coefficient
    corr_loc_df <- corr_df %>% 
      dplyr::filter(Sample1 == corr_df[[i, "Sample1"]] & Sample2 %in% corr_df[[i, "Sample2"]]) %>%
      dplyr::arrange(R, Sample1, Sample2)
    
    # Calculate linear regression and residuals, identify outliers
    data_loc_df <- lm_residual_outliers(data_loc_df, corr_df[[i, "Sample1"]], corr_df[[i, "Sample2"]])

    # Create a scatter plot of 2 most dissimilar samples within the cohort
    filename <- paste0(output_path, "scatter_pair_", corr_df[[i, "Sample1"]], "_", corr_df[[i, "Sample2"]])
    create_scatter_plot_my(data_loc_df, x_label=corr_df[[i, "Sample1"]], y_label=corr_df[[i, "Sample2"]], 
                           fill_label="IsOutlier", filename=filename, 
                           labels=c(corr_df[[i, "Sample1"]], corr_df[[i, "Sample2"]], paste("R =", round(corr_df[[i, "R"]], 3))),
                           size=0.5, do_fit=TRUE, stroke=0)
    
    filename <- paste0(output_path, "scatter_pair_", corr_df[[i, "Sample1"]], "_", corr_df[[i, "Sample2"]], ".txt")
    readr::write_delim(data_loc_df, filename, delim = "\t")
    
  }
  
}
