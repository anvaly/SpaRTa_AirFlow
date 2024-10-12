# Author: Anna Lyubetskaya. Date: 20-04-22
# Compare signature scores between clusters to characterize clusters better
# Useful link: https://drsimonj.svbtle.com/exploring-correlations-in-r-with-corrr


## ENVIRONMENT ----

setwd("/repos/P02567_TBIO-3021_10X_pilot/")
library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

# source("code/utils/utils_stats.R")
source("code/R/Utils/utils_10X_signatures.R")
source("code/R/c_Signatures/ST_analysis_signatures_corr_utils.R")


## PARAMETERS ----


# Path to processed Seurat data
sample_name <- "pdac_cohort_spotclean_harmony_integrated"

# SCT and spot number thresholds for gene filtering
sct_threshold <- 1
spot_threshold <- 5

# Slot to use to calculate signature score: data or scale.data
slot <- "data"

# Select clustering resolution to work with; if NULL, look in the misc slot
resolution <- "SCT_snn_res.0.5"  # "integrated_snn_res.0.2"


## PATHS ----


# Path to a signature file
#sig_path <- "data/import/Signatures/signatures_colon_t100_Mar21_short.txt"
sig_path <-  "/repos/P02567_TBIO-3021_10X_pilot/data/import/Signatures/signatures_220210.txt"

# Location of pre-processed data
input_path <- paste0("XXXX")
#input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, sample_name, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Seurat data
data_seurat <- readRDS(input_path)

if(is.null(resolution)){
  resolution <- data_seurat@misc$user.Clustering
}


## WRANGLE DATA ----


# Calculate signature scores and add them to Seurat meta data
data_seurat <- annotate_seurat_with_signatures(data_seurat, sig_path, 
                                               sct_threshold=sct_threshold, spot_threshold=spot_threshold)

# Find column names in the Seurat object
sig_names <- colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))])


## EXTRACT SIGNATURE DATA ----


# Wide tibble of signature scores
sig_wide_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate") %>%
  dplyr::select(dplyr::all_of(c("Coordinate", sig_names))) %>%
  dplyr::mutate_if(is.numeric, round, 2)

# Create a long tibble of signature scores
sig_df <- sig_wide_df %>%
  df_wide2long_my(key="Signature_name", val="Score")

# Extract clustering meta data
res_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate") %>%
  dplyr::select(dplyr::all_of(c("Coordinate", resolution)))

# Add cluster data to the long tibble
sig_clust_df <- sig_df %>%
  dplyr::inner_join(res_df, by="Coordinate")


## CALCULATE SIGNATURE STATS BY CLUSTER ----


## USE ANOVA TO SELECT SIGNATURES ENRICHED IN CLUSTERS ----
# https://www.datanovia.com/en/lessons/anova-in-r/


# Perform one-sided ANOVA to compare signature scores by cluster
# Couldn't make the code work without the loop
anova_pval_list <- list()
for(sig in sig_names){
  anova_df <- sig_clust_df %>%
    dplyr::filter(Signature_name == sig) %>%
    dplyr::rename(resolution = !!rlang::sym(resolution)) %>%
    dplyr::select(resolution, Score)
  
  anova_pval_list[[sig]] <- rstatix::anova_test(anova_df, Score ~ resolution)$p
}
