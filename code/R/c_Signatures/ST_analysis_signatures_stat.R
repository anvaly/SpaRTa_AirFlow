# Author: Anna Lyubetskaya. Date: 20-04-22

# Calculate various statistics for a set of signatures in a set of Seurat objects
# Stats include: min, mean, sd, max


## SETUP ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`


source("code/utils/utils_signatures.R")
source("code/utils/utils_ggplot.R")
source("code/utils/utils_heatmap.R")

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_signatures.R")


## PARAMETERS ----


# The cohort of interest regex ID
# Seurat RDS files are tagged as follows
cohort_name <- "PDAC108_path14"
cohort_regex <- "PDAC"

# Exclude specific samples from the analysis
samples_exclude <- NULL

# Gene abundance filters
sct_threshold <- 0.5
spot_threshold <- 5

# Signatures to select
sig_select <- c("PDAC.collisson.classical", "PDAC.moffitt.basal", "PDAC.moffitt.activatedstroma", "PDAC.moffitt.normalstroma",
                "PDAC.P19.Ductal_1", "PDAC.P19.Acinar", "PDAC.P19.Endocrine", "PDAC.P19.Endothelial", "PDAC.P19.Stellate", 
                "PDAC.P19.Bcell", "PDAC.P19.Fibroblast", "PDAC.P19.Macrophage", "PDAC.P19.Tcell",
                "BMS.Pathway.TNFa", "BMS.Pathway.TGFB", "BMS.Pathway.IFNg", "BMS.Pathway.IFNa",
                "BMS.PDAC.Hypoxia.CL", "PDAC.U.Immune.Tcell.exhausted"
)


## PATHS ----


# Input folder
input_paths <- c("XXXX")

# Output information
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_name, "/")

# Create output folder
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)

# Path to a signature file
sig_path <- "data/import/Signatures/signatures_220210.txt"


## INGEST DATA ----


# Find all files following the regex
file_list <- unname(unlist(sapply(input_paths, function(p) dir(p, pattern=paste0(cohort_regex, ".*.rds"), full.names=TRUE))))

# Identify samples user wants to exclude and remove them from the file list
if(!is.null(samples_exclude)){
  samples_exclude_indx <- c(unname(sapply(samples_exclude, function(x) grep(x, file_list))))
  file_list <- setdiff(file_list, file_list[samples_exclude_indx])
}

# Ingest a set of RDS Seurat objects
seurat_list <- read_rds_list_simple_my(file_list)


## WRANGLE AND VISUALIZE DATA ----


# Gather signature scores and stats into a list
stat_list <- list()

# Cycle through Seurat datasets
for(dataset in names(seurat_list)){
  
  # Select a dataset
  data_seurat <- seurat_list[[dataset]]
  spot_num <- ncol(data_seurat)
  
  
  ## Calculate signature scores and add them to Seurat meta data ----
  
  
  # Find genes abundant in this sample
  gene_list <- seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold,
                                               assay="SCT", slot="data")
  
  # Load signatures and filter them down to only well represented genes
  signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select, 
                                              sig_length_min=1, sig_length_max=1000, ratio_threshold=0)
  
  # Signature names in the object
  sig_names <- names(signature_list)
  
  # Invert signature list
  sig_invert_list <- invert_list_my(signature_list)
    
  # Add signature scores to a seurat object
  data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.", assay="SCT")
  
  # Seurat renames column names, this step finds new signature names
  sig_names_upd <- colnames(data_seurat@meta.data)[grep("sig.", colnames(data_seurat@meta.data))]
  
  
  ## Calculate signature stats ----
  
  
  # Wide tibble of signature scores
  sig_wide_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate") %>%
    dplyr::select(dplyr::all_of(c("Coordinate", sig_names_upd)))
  
  # Create a long tibble of signature scores
  sig_df <- sig_wide_df %>%
    df_wide2long_my(key="Signature_name", val="Score")
  
  # Calculate various signature stats
  sig_stat_df <- sig_df %>%
    dplyr::group_by(Signature_name) %>%
    dplyr::summarise(score_min = min(Score),
                     score_mean = mean(Score),
                     score_median = median(Score),
                     score_sd = sd(Score),
                     score_max = max(Score),
                     spot_above0.5_num = round(sum(Score >= score_median + 0.5) / spot_num * 100),
                     spot_above1sd_num = round(sum(Score >= score_median + score_sd * 1) / spot_num * 100)) %>%
    dplyr::mutate_if(is.numeric, round, 2) %>%
    dplyr::arrange(desc(score_mean), desc(score_max)) %>%
    dplyr::mutate(HighMax = score_max >= 2,
                  Sample_Name = data_seurat@misc$user.Sample_Name)
  

  ## Extract expression data ----
  
  
  # Extract expression matrix, transform it into a long tibble, and calculate z-score for each gene across all barcodes
  data_df <- seurat_expression_to_long_tibble_my(data_seurat, assay="SCT", slot="data")
  
  # Top 20 highest expressed genes in the signature
  gene_stat_df <- data_df %>%
    dplyr::filter(Symbol %in% sig_invert_list$Symbol) %>%
    dplyr::group_by(Symbol) %>%
    dplyr::summarise(score_mean = mean(SCT_data)) %>%
    dplyr::filter(score_mean >= 0) %>%
    dplyr::arrange(desc(score_mean))

  # Sort and select signature gene lists
  signature_select_list <- lapply(names(signature_list), function(x) gene_stat_df %>% 
           dplyr::filter(Symbol %in% signature_list[[x]]) %>% 
           dplyr::pull(Symbol) %>%
           paste(collapse=";"))
  
  # Signature lists sorted by individual gene abundance
  signature_upd_list <- tibble::tibble(Signature_name = paste0("sig.", names(signature_list)),
                                       Gene_list = unlist(signature_select_list))
  
  # Add gene names to the tibble with stats
  sig_stat_df <- sig_stat_df %>%
    dplyr::left_join(signature_upd_list, by="Signature_name")
  
  # Gather signature stats into a list
  stat_list[[data_seurat@misc$user.Sample_Name]] <- sig_stat_df
  
  
  ## Outputs and visualizations ----
  
  
  filename <- paste0(output_path, "table_sig_stats_", dataset, ".txt")
  readr::write_delim(sig_stat_df, filename, delim="\t")
  
  # Select top 50 signatures
  #sig_stat_df <- sig_stat_df %>%
  #  dplyr::top_n(n=50, wt=spot_above1_num)
  
  # Mean scores with stdev
  filename <- paste0(output_path, "sig_hist_", dataset)
  create_bar_plot_my(sig_stat_df, x_label="Signature_name", y_label="score_mean", 
                     fill_label="HighMax", filename=filename, 
                     labels=c("Signature name", "Signature score mean and standard deviation", 
                              "Signature stats distribution"), 
                     error_label="score_sd", reorder_x=TRUE)

  # Signature score box plots
  filename <- paste0(output_path, "sig_box_", dataset)
  create_box_plot_my(sig_df, x_label="Signature_name", y_label="Score", 
                     fill_label="Signature_name", filename=filename, 
                     labels=c("Signature name", "Signature score", "Signature score distribution"),
                     reorder_x=TRUE)
}


# Join signature stats across samples
stat_df <- dplyr::bind_rows(stat_list) %>%
  dplyr::arrange(score_mean)

filename <- paste0(output_path, "joint_table_sig_stats.txt")
readr::write_delim(stat_df, filename, delim="\t")


# Create a box plot of gene signature levels by sample
filename <- paste0(output_path, "joint_signature_boxplot")

# Signature means
p1 <- create_box_plot_my(stat_df, x_label="Signature_name", y_label="score_mean", fill_label="Signature_name", 
                         labels=c("Signature", "Mean signature score", "Mean score across spots in a sample"), 
                         filename=NULL, with_dots=FALSE, outlier_shape=NULL, reorder_x=FALSE)

# Signature medians
p2 <- create_box_plot_my(stat_df, x_label="Signature_name", y_label="score_median", fill_label="Signature_name", 
                         labels=c("Signature", "Median signature score", "Median score across spots in a sample"), 
                         filename=NULL, with_dots=FALSE, outlier_shape=NULL, reorder_x=FALSE)

# Number of spots with a signature score higher than median + 0.5
p3 <- create_box_plot_my(stat_df, x_label="Signature_name", y_label="spot_above0.5_num", fill_label="Signature_name", 
                         labels=c("Signature", "Spot # above threshold", "# spots with signature score higher than median + 0.5"),
                         filename=NULL, with_dots=FALSE, outlier_shape=NULL, reorder_x=FALSE)

# Number of spots with a signature score higher than median + 1SD
p4 <- create_box_plot_my(stat_df, x_label="Signature_name", y_label="spot_above1sd_num", fill_label="Signature_name", 
                         labels=c("Signature", "Spot # above threshold", "# spots with signature score higher than median + 1SD"), 
                         filename=NULL, with_dots=FALSE, outlier_shape=NULL, reorder_x=FALSE)

write_plot2file_my(patchwork::wrap_plots(list(p1, p2, p3, p4), ncol=2, nrow=2), filename, num_row=2, num_col=2)


# Create a heatmap of signature scores by sample
params <- list(cell_value = "spot_above0.5_num",
               row_label = "Sample_Name", 
               col_label = "Signature_name", 
               distance = "pearson",
               range = c(0, median(stat_df$spot_above0.5_num), round(median(stat_df$spot_above0.5_num) + sd(stat_df$spot_above0.5_num) * 3)),
               colors = c("white", "yellow", "royalblue4"),
               show_column_dend = TRUE,
               show_row_dend = TRUE)

filename <- paste0(output_path, "joint_signature_heatmap.png")
create_heatmap_my(stat_df %>%
                    dplyr::filter(spot_above0.5_num > 0),
                  params, row_list=NULL, col_list=NULL, col_meta_df=NULL, row_meta_df=NULL, filename=filename)
