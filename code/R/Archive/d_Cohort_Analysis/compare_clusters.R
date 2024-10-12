# Author: Anna Lyubetskaya. Date: 20-08-03


##_ SETUP ENVIRONMENT _##


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_in_out.R")
source("code/utils/utils_heatmap.R")
source("code/utils/utils_signatures.R")
source("code/utils/utils_specialized_plots.R")

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_matrix.R")

source("code/R/Cohort_Analysis/compare_clusters_utils.R")


## PARAMETERS ----


# The cohort of interest regex ID
cohort_name <- "FF_HumanPanc"
# Substitute original cluster for user defined or not
do_sub_clusters <- FALSE

cluster_order_list <- rev(c(2, 6, 9, 1, 5, 3, 4, 7, 8, 0))  # HumanPanc FF
#cluster_order_list <- rev(c(6, 2, 1, 4, 7, 8, 5, 0, 3))  # HumanPanc FFPE


## PATHS ----


# Input folder
input_folder <- "C:/Users/lyubetsa/Documents/Projects_NewTech/10X_ST_Reports/"

# DEA files are tagged as follows
dea_regex_files <- "_dea.csv"
# Seurat RDS files are tagged as follows
data_regex_files <- "_all.Sobj.rds"

# Output folder
output_folder <- "C:/Users/lyubetsa/Documents/Projects_NewTech/10X_ST_Cluster/"

# DEA thresholds
fc_threshold <- 1
fdr_threshold <- 2  # in -log10(FDR) space


## INGEST DATA ----


# Ingest clean Seurat data

seurat_list <- read_seurat_rds_my(input_folder, file_tag=paste0(cohort_name, ".*", data_regex_files), suffix=data_regex_files)
names(seurat_list) <- gsub(data_regex_files, "", names(seurat_list))

# Read DEA data
dea_df <- read_dir2file_my(input_folder, in_regex=paste0(cohort_name, ".*", dea_regex_files), do_search=TRUE) %>%
  dplyr::mutate(Sample = gsub("^.+/|_dea.csv", "", File)) %>%
  dplyr::select(-File)


## WRANGLE DATA ----


# Select genes that pass DEA thresholds
dea_filt_df <- dea_df %>%
  dplyr::filter(avg_logFC >= fc_threshold & p_val_adj_neg_log10 >= fdr_threshold)

# List of DEA genes
dea_gene_list <- dea_filt_df %>% 
  dplyr::pull(Symbol) %>% 
  unique

# Extract gene expression values for select genes from each sample in the cohort
expression_list <- list()
for(n in names(seurat_list)){
  # Seurat clusters information
  clusters_df <- tibble::as_tibble(seurat_list[[n]]@meta.data, rownames="id") %>%
    dplyr::rename(Coordinate = id) %>%
    dplyr::select(Coordinate, seurat_clusters)
  
  # Long tibble of expression values
  sample_expression_df <- seurat_expression_to_long_tibble_my(seurat_list[[n]], assay="SCT")
  
  # Select genes of interest
  # Add seurat cluster annotation
  expression_list[[n]] <- sample_expression_df %>%
    dplyr::filter(Symbol %in% dea_gene_list) %>%
    dplyr::mutate(Sample = n) %>%
    dplyr::inner_join(clusters_df, by="Coordinate")
}

# Bind all expression data into a single long tibble
expression_df <- dplyr::bind_rows(expression_list) %>%
  dplyr::mutate(Sample = gsub("P-\\d+-\\d+_ST_", "", Sample)) %>%
  dplyr::mutate(Sample_cluster = paste(Sample, seurat_clusters),
                Cohort = gsub("_.+$", "", Sample_cluster))

# Calculate mean expression zscore for DE genes for each cluster in each sample
cluster_profiles_df <- expression_df %>%
  dplyr::group_by(Sample_cluster, Symbol, Cohort, Sample, seurat_clusters) %>%
  dplyr::summarise(SCT_zscore_mean = mean(SCT_zscore)) %>%
  dplyr::ungroup()


## LOAD AND FILTER SIGNATURE DATA ----


# Path to signatures
sig_path <- "C:/Users/lyubetsa/Documents/Data/Signatures/processed_integrated_panc_Oct2020.txt"

# Load and filter signatures to only those genes that are well represented 
# Note: These signatures are filtered down to contain only DE genes
signature_list <- read_filter_signatures_my(sig_path, unique(dea_df$Symbol),
                                            sig_length_max=1000, ratio_threshold=0)  # dea_gene_list

# Annotate DE genes by their signatures
gene_sig_df <- invert_list_my(signature_list) %>%
  dplyr::mutate(Sig_Name = gsub("drokhle_|brownea_|kumarn_|panc_|prostate_|siemersn_|", "", Sig_Name)) %>%
  dplyr::right_join(dea_filt_df %>% 
                     dplyr::select(Symbol, pct_1) %>% 
                     dplyr::group_by(Symbol) %>% 
                     dplyr::summarise(pct_mean = mean(pct_1)), by="Symbol") %>%
  tidyr::replace_na(list(Num_Sigs = 0)) %>%
  dplyr::mutate(InSignature = Num_Sigs > 0)

# Write signature stats to file
readr::write_delim(gene_sig_df, path=paste0(output_folder, "/gene_sig_annotations.txt"), 
                   delim = "\t", append=FALSE, col_names = TRUE)


## CORRELATE SAMPLES TO CONSOLIDATE CLUSTERS ----


# Correlate sample clusters using DE genes
filename <- paste0(output_folder, cohort_name, "_corr_markers_zscore_mean")
correlate_samples_my(cluster_profiles_df, filename, cluster_col="Sample_cluster")


## HEATMAP OF SAMPLES USING DEA GENES ----


# Manually-curated dictionary for cluster name consolidation
consolidation_dict_file <- paste0(output_folder, cohort_name, "_cluster_consolidation.txt")

if(!is.null(consolidation_dict_file)){
  
  # Ingest manually curated cluster dictionary
  cluster_rename_df <- read_file2df_my(consolidation_dict_file, delim="\t")

  # Calculate mean expression zscore for DE genes for each cluster in each sample
  # Merge cluster data with the manually curated labels
  if(do_sub_clusters == TRUE){
    cluster_profiles_plot_df <- expression_df %>%
      dplyr::inner_join(cluster_rename_df, by="Sample_cluster") %>%
      dplyr::mutate(my_clusters = gsub("^.+ ", "", Sample_cluster_my)) %>%
      dplyr::group_by(Sample_cluster_my, Symbol, Cohort, Sample, my_clusters) %>%
      dplyr::summarise(SCT_zscore_mean = mean(SCT_zscore)) %>%
      dplyr::ungroup() %>%
      dplyr::rename(Sample_cluster = Sample_cluster_my)
  } else{
    cluster_profiles_plot_df <- cluster_profiles_df %>%
      dplyr::inner_join(cluster_rename_df, by="Sample_cluster") %>%
      dplyr::mutate(my_clusters = gsub("^.+ ", "", Sample_cluster_my))
  }
  
  # Perform clustering with new labels
  filename <- paste0(output_folder, cohort_name, "_corr_markers_zscore_mean_upd")
  correlate_samples_my(cluster_profiles_plot_df, filename, cluster_col="Sample_cluster")
  
  vec_new <- c(sapply(cluster_order_list, function(x) paste(x, unique(cluster_profiles_plot_df$Sample))))
  vec_new <- unname(sapply(vec_new, function(x) paste(unlist(strsplit(x, " "))[2], unlist(strsplit(x, " "))[1])))
  
  # Heatmap parameters
  params <- list(cell_value = "SCT_zscore_mean",
                 row_label = "Symbol", 
                 col_label = "Sample_cluster", 
                 distance = "pearson",
                 row_annotation = c("InSignature"),
                 col_annotation = c("Sample", "my_clusters", "Cohort"),
                 range = c(-2, 0, 2),
                 colors = c("blue4", "white", "red4"),
                 row_expected_groups = 10,
                 col_expected_groups = 10,
                 column_order = intersect(vec_new, cluster_profiles_plot_df$Sample_cluster)
                 )
  
  # Column meta data
  col_meta_df <- cluster_profiles_plot_df %>% 
    dplyr::select(dplyr::all_of(c(params$col_label, params$col_annotation))) %>% 
    unique
  
  # Create heatmap of sample clusters using DE genes
  filename <- paste0(output_folder, cohort_name, "_hm_markers_zscore_mean_upd.png")
  hm <- create_heatmap_my(cluster_profiles_plot_df, params, row_meta_df = gene_sig_df, col_meta_df = col_meta_df, filename=filename)

  hm_gene_list <- rownames(hm@ht_list$SCT_zscore_mean@matrix)
  hm_row_index <- ComplexHeatmap::row_order(hm)
  names(hm_row_index)
  hm_markers_df <- invert_list_my(sapply(hm_row_index, function(x) hm_gene_list[x])) %>%
    dplyr::inner_join(gene_sig_df, by="Symbol")
  b = hm_markers_df %>% dplyr::ungroup() %>% dplyr::group_by(Sig_Name.x) %>% dplyr::summarise(group = paste(unique(Sig_Name.y)))
}
