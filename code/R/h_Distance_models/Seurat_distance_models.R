# Author: Anna Lyubetskaya. Date: 20-04-22

# This script compares normalized gene expression values in each spot to a select spot distance metric (e.g., distance between a spot and tumor center)
# The comparison is performed in the following ways:
# - Correlation between distance and each gene expression profile across spots
# - Linear regression fit: GeneX expression ~ Distance + Cell counts
# - ElNet feature selection: Distance ~ GeneA expression + ... + GeneZ expression + Cell counts
# Differences between tissue sections can be accounted by including section identity as a model covariate or by normalizing the distance metric to within the section


## ENVIRONMENT ----


library(Seurat)
library(foreach)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`


source("code/utils/utils_ggplot.R")
source("code/utils/utils_signatures.R")
source("code/utils/utils_lin_regress.R")

source("code/R/Utils/utils_10X_matrix.R")


## PARAMETERS ----


# Path to processed Seurat data
# sample_name <- "Syng_B16_cca_sct"
sample_name <- "PDAC_FFPE-probes_merge_dist"
sample_exclude <- NULL  # c("HumanPanc_ROI1_FFPE_A_Apr21")  # c("HumanPanc_ROI2_FFPE_D_Dec20_T")

# Distance variable to use for correlations
dist_var <- "Distance"  # Pathology.Distance.Tissue.filled; Pathology.Distance.Epithelium.filled.invertTRUE
# Normalize pathology distances by sample
normalize_dist_by_sample <- FALSE

# Variables to include in the ElNet model as confounders
confounder_vars_num <- NULL  # c("CellCounts")
confounder_vars_cat <- "user.Sample_Name"  # e.g., NULL or c("user.Sample_Name", "user.Tissue")

# Filters to select genes
sct_threshold <- 0.5
spot_threshold <- 5

# List of genes to highlight in the scatter plot
# MC38
#gene_highlight_list <- c("CD74", "H2-AB1", "COL1A2", "SRGN", "ARG1", "H2-AA", "H2-EB1", "CRIP1", "COL1A1", "FCRLS", "COL3A1")
# B16
#gene_highlight_list <- c("VIM","CAR6","COPS9","FTL1","FOS","S100A6","TMSB4X","LGALS1","ATP5MPL","FTH1","HBB-BS","FN1")
gene_highlight_list <- NULL


## PATHS ----


# Path to a signature file
#sig_path <- "data/import/Signatures/signatures_syngeneics_t100_Aug21.txt"
sig_path <- "data/import/Signatures/signatures_pdac_t1000_May21.txt"

# Location of pre-processed data
input_path <- paste0("XXXX")
#input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, sample_name, "_", dist_var, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Seurat data
data_seurat <- readRDS(input_path)

# Find samples in the merged object that don't satisfy the gene selection criteria
sample_exclude <- c(sample_exclude, names(table(data_seurat@meta.data$user.Sample_Name))[which(table(data_seurat@meta.data$user.Sample_Name) <= spot_threshold)])

# Remove samples if necessary
if(!is.null(sample_exclude)){
  barcode_list <- rownames(data_seurat@meta.data[which(!data_seurat@meta.data[["user.Sample_Name"]] %in% sample_exclude),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# List of samples and number of spots in each
table(data_seurat@meta.data$user.Sample_Name)

# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, NULL, sig_length_min=1, sig_length_max=1000, ratio_threshold=0)

# Invert signatures to get an annotated gene list
sig_gene_df <- invert_list_my(signature_list)


## WRANGLE DATA ----


# Abundant gene list
gene_list <- unique(seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold, 
                                                    assay="SCT", slot="data", split_by="user.Sample_Name"))

# Exclude MT/Ribo genes for ease of interpretation
gene_list <- setdiff(gene_list, gene_list[grep("^MT-|^RP[SL]|^GM\\d+", gene_list)])
# Remove NA value
gene_list <- gene_list[which(!is.na(gene_list))]

# Extract a wide matrix of SCT normalized expression values
data_wide_df <- tibble::as_tibble(t(as.matrix(Seurat::GetAssayData(data_seurat, assay="SCT", slot="data"))), rownames="Coordinate")

# Subset data to only relevant genes
data_wide_df <- data_wide_df %>%
  dplyr::select(dplyr::all_of(intersect(c("Coordinate", gene_list), colnames(data_wide_df))))

# Check that metadata and expression data have the same order our of paranoia
table(rownames(data_seurat@meta.data) == data_wide_df$Coordinate)

# Extract the relevant meta data
meta_data <- data_seurat@meta.data[c(dist_var, confounder_vars_num, confounder_vars_cat)] %>%
  dplyr::filter(!is.na(!!rlang::sym(dist_var)) & !is.infinite(!!rlang::sym(dist_var))) %>%
  tibble::rownames_to_column("Coordinate")

# Re-normalize pathology distance by sample (% of max)
if(normalize_dist_by_sample == TRUE){
  dist_norm_df <- data_seurat@meta.data[c(dist_var, "user.Sample_Name")] %>%
    dplyr::group_by(user.Sample_Name) %>%
    dplyr::mutate(MaxDist = max(!!rlang::sym(dist_var)),
                  DISTANCE = !!rlang::sym(dist_var) / MaxDist * 100)
  
  meta_data[dist_var] <- dist_norm_df$DISTANCE
}

# Add meta data to the expression data
data_wide_df <- data_wide_df %>%
  dplyr::inner_join(meta_data, by="Coordinate")


## FIND DISTANCE CORRELATES ----


# Calculate correlation between the distance and each of the gene expression vectors
cor_vec <- sapply(gene_list, function(x) cor(data_wide_df[[x]], meta_data[[dist_var]]))
names(cor_vec) <- gene_list

# Create a correlation tibble
cor_df <- tibble::tibble(Symbol = gene_list,
                         R2 = cor_vec) %>%
  dplyr::mutate(R2 = round(R2, 3)) %>%
  dplyr::arrange(R2) %>%
  dplyr::left_join(sig_gene_df, by="Symbol")

# Write correlations to file
filename <- paste0(output_path, "correlation_result.txt")
readr::write_delim(cor_df, filename, delim="\t")

# Create a correlation histogram
filename <- paste0(output_path, "correlation_hist")
create_hist_plot_my(cor_df, x_label="R2", fill_label="InSignature", intercept=c(-0.25, 0, 0.25), 
                    binwidth=0.01, filename=filename, 
                    labels=c("R2", "Gene number", paste0("Correlation with", dist_var)))

# Cleanup before running more analysis
rm(data_seurat)
gc()


## ELNET OF EXPRESSION ON DISTANCE ----


## Prep data for ElNet

# Define the outcome and explicitly define row names
data_wide_df <- data_wide_df %>% 
  dplyr::select(-Coordinate)

# Factorize factor fields
if(!is.null(confounder_vars_cat) && length(confounder_vars_cat) > 0){
  data_wide_df[confounder_vars_cat] <- lapply(data_wide_df[confounder_vars_cat], as.factor) 
}


# Save gene ElNets to a file as a list of models
filename <- paste0(output_path, "gene_elnet_result.rds")
gene_model_list <- list()

# Go through all genes and perform individual fits: Gene Expression ~ Distance to Center
if(!file.exists(filename)){
  
  # Parallelize the process for speed
  cl <- parallel::makeCluster(parallel::detectCores() / 2)
  doParallel::registerDoParallel(cl)
  
  gene_model_list <- foreach(gene = gene_list, .combine='c') %dopar% {

    # Define the outcome and explicitly define row names
    data_wide_loc_df <- data_wide_df %>%
      dplyr::select(dplyr::all_of(c(dist_var, gene, confounder_vars_cat, confounder_vars_num))) %>% 
      dplyr::rename("outcome" = !!rlang::sym(gene))
    
    # Run ElNet
    model_res <- run_all_models_my(data_wide_loc_df, c("norm"), output_loc=NULL, formula=NULL)
    list(gene = summary(model_res$norm)$coefficients)

  }
  
  parallel::stopCluster(cl)
  
  # Add gene names to list names
  names(gene_model_list) <- gene_list
  
  # Save a list of models to file
  saveRDS(gene_model_list, filename)
} else{
  gene_model_list <- readRDS(filename)
}

# Analyze the result of the linear model
coef_gene_df <- t(sapply(names(gene_model_list), function(g) gene_model_list[[g]][dist_var, c(1,4)])) %>%
  tibble::as_tibble(rownames="Symbol") %>%
  dplyr::mutate(Score = -log10(`Pr(>|t|)`),
                Labels = ifelse(Symbol %in% gene_highlight_list,
                                Symbol, ""))

# Establish the estimate threshold as mean +/- 3SD
estimate_threshold <- mean(abs(coef_gene_df$Estimate)) + sd(abs(coef_gene_df$Estimate)) * 3

# Call significant coefficients
coef_gene_df <- coef_gene_df %>%
  dplyr::mutate(IsSignificant = `Pr(>|t|)` <= 0.05 / length(gene_list) & abs(Estimate) >= estimate_threshold) %>%
  dplyr::arrange(desc(IsSignificant), desc(Score))

# Visualize the result of linear modeling of gene expression ~ distance to center
filename <- paste0(output_path, "gene_elnet_result_scatter")
create_scatter_plot_my(coef_gene_df, x_label="Estimate", y_label="Score", 
                       fill_label="IsSignificant", shape=21, size=2, dot_labels="Labels", 
                       filename=filename, labels=NULL, do_fit=NULL, stroke=0)

filename <- paste0(output_path, "gene_elnet_result_scatter.txt")
readr::write_delim(coef_gene_df, filename, delim="\t")


## ELNET OF DISTANCE ON EXPRESSION ----


# Define the outcome and explicitly define row names
data_wide_loc_df <- data_wide_df %>% 
  dplyr::rename("outcome" = !!rlang::sym(dist_var))

# Run ElNet
filename <- paste0(output_path, "dist_elnet_result.rds")
if(!file.exists(filename)){
  #model_res <- run_all_models_my(data_wide_loc_df, c("elnet"), output_loc=NULL, formula=NULL)
  #saveRDS(model_res, filename)
} else{
  model_res <- readRDS(filename)
}

model_res$elnet

# Analyze feature importance
# https://topepo.github.io/caret/variable-importance.html
elnet_imp <- caret::varImp(model_res$elnet)

# Analyze the result of ElNet
coef_df <- analyze_model_my(model_res$elnet, analysis_type="elnet", lambda_type="lambda.1se") %>%
  dplyr::left_join(tibble::as_tibble(elnet_imp$importance, rownames="variable"), by="variable") %>%
  dplyr::left_join(sig_gene_df, by=c("variable" = "Symbol")) %>%
  dplyr::mutate(Overall = round(Overall, 2))

# Number of top features to visualize
n <- 30

# Visualize the result of ElNet
filename <- paste0(output_path, "dist_elnet_result")
coef_df <- visualize_model_my(coef_df, filename, feature_list=gene_list, n=n)

# Write correlations to file
filename <- paste0(output_path, "dist_elnet_result_coefficients_filt.txt")
readr::write_delim(coef_df %>%
                     dplyr::arrange(-Overall) %>%
                     dplyr::filter(Overall >= 5), filename, delim="\t")

# Plot importance
filename <- paste0(output_path, "dist_elnet_result_imp_bar")
p <- create_bar_plot_my(coef_df %>%
                          dplyr::filter(Overall >= 25), x_label="variable", y_label="Overall", 
                        fill_label="IsLineage", filename=filename, reorder_x=TRUE, 
                        labels=c("Feature", "Importance", "ElNet, coefficient importance"))

# Default feature imporance visualization
filename <- paste0(output_path, "dist_elnet_result_imp_line.png")
png(filename, width = 10, height = 5, units = "in", res = 300)
plot(elnet_imp, top = n)
dev.off()
