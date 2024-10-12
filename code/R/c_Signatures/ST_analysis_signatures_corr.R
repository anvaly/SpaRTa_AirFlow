# Author: Anna Lyubetskaya. Date: 20-04-22
# Perform correlative analysis between a set of signatures
# Useful link: https://drsimonj.svbtle.com/exploring-correlations-in-r-with-corrr


## ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`


source("code/utils/utils_ggplot.R")
source("code/utils/utils_heatmap.R")

source("code/R/Utils/utils_10X_signatures.R")

source("code/R/c_Signatures/ST_analysis_signatures_corr_utils.R")


## PARAMETERS ----


# Output name ID
cohort_name <- "PDAC108_path14_P19ChijiCosMx_ADM"

# Regex to gather objects
cohort_regex <- "PDAC"

# Exclude certain samples from analysis
samples_exclude <- NULL

# SCT and spot number thresholds for gene filtering
sct_threshold <- 0.5
# Number of spots threshold
spot_threshold <- 5

# Signature filter: filter out signature that don't ever reach this intensity
sig_threshold <- 0.25

# Correlation threshold
corr_threshold <- 0
# Cluster size
cluster_size <- 10

# Assay and slot to use to calculate signature score
assay <- "SCT"
slot <- "data"

# Signatures to select
sig_select <- c("PDAC.collisson.classical", "PDAC.moffitt.classical", "PDAC.P19.Ductal_2",
                "PDAC.moffitt.basal", "PDAC.collisson.quasimesenchymal", 
                "PDAC.moffitt.activatedstroma", "PDAC.P19.Fibroblast", 
                "PDAC.moffitt.normalstroma", "PDAC.P19.Stellate", 
                "BMS.PDAC.Hypoxia.CL", "Syng.U.State.Hypoxia.metabolism", "PDAC.P19.Endothelial", 
                "PDAC.P19.Acinar", "PDAC.P19.Endocrine", "PDAC.P19.Ductal_1", 
                "PDAC.P19.Bcell", "PDAC.P19.Macrophage", "PDAC.P19.Tcell",
                "BMS.Pathway.TNFa", "BMS.Pathway.TGFB", "BMS.Pathway.IFNg", "BMS.Pathway.IFNa", 
                "Oxford.Tcell_exhausted",
                "PDAC.CosMx.Bcell", "PDAC.CosMx.TumorBasal",
                "PDAC.CosMx.TumorClassical", "PDAC.CosMx.TumorClassicalCycling",
                "PDAC.CosMx.Ductal", "PDAC.CosMx.Endocrine",
                "PDAC.CosMx.Endothelial", "PDAC.CosMx.Exocrine",
                "PDAC.CosMx.Fibroblast", "PDAC.CosMx.Lowquality",
                "PDAC.CosMx.Macrophage", "PDAC.CosMx.Mast",
                "PDAC.CosMx.MDSC", "PDAC.CosMx.TumorMixed",
                "PDAC.CosMx.TumorOther", "PDAC.CosMx.Plasma", "PDAC.CosMx.Tcell",
                "PDAC.Chiji.Acinar", "PDAC.Chiji.Bcell", "PDAC.Chiji.Ductal_1",
                "PDAC.Chiji.Ductal_2", "PDAC.Chiji.Endocrine", "PDAC.Chiji.Endothelial",
                "PDAC.Chiji.Fibroblast", "PDAC.Chiji.Macrophage", "PDAC.Chiji.Stellate",
                "PDAC.Chiji.Tcell", "PDAC.Chiji.TumorClassical", "PDAC.Chiji.TumorBasal"
)

# Rename signatures for visualizations
sig_rename <- c("collisson.classical", "moffitt.classical", "P19.Ductal_2",
                "moffitt.basal", "collisson.quasimesenchymal", 
                "moffitt.activatedstroma", "P19.Fibroblast", 
                "moffitt.normalstroma", "P19.Stellate", 
                "Hypoxia.CL", "Hypoxia.metabolism", "P19.Endothelial", 
                "P19.Acinar", "P19.Endocrine", "P19.Ductal_1", 
                "P19.Bcell", "P19.Macrophage", "P19.Tcell",
                "Pathway.TNFa", "Pathway.TGFB", "Pathway.IFNg", "Pathway.IFNa", 
                "Tcell_exhausted",
                "CosMx.Bcell", "CosMx.TumorBasal",
                "CosMx.TumorClassical", "CosMx.TumorClassicalCycling",
                "CosMx.Ductal", "CosMx.Endocrine",
                "CosMx.Endothelial", "CosMx.Exocrine",
                "CosMx.Fibroblast", "CosMx.Lowquality",
                "CosMx.Macrophage", "CosMx.Mast",
                "CosMx.MDSC", "CosMx.TumorMixed",
                "CosMx.TumorOther", "CosMx.Plasma", "CosMx.Tcell",
                "Chiji.Acinar", "Chiji.Bcell", "Chiji.Ductal_1",
                "Chiji.Ductal_2", "Chiji.Endocrine", "Chiji.Endothelial",
                "Chiji.Fibroblast", "Chiji.Macrophage", "Chiji.Stellate",
                "Chiji.Tcell", "Chiji.TumorClassical", "Chiji.TumorBasal"
)

# Add one-off genes to the plots
# Adding a "Gene." prefix so that AddModuleScore doesn't protest
gene_select_symbol <- NULL
gene_select <- NULL

# Specific pairs of signatures to investigate
# NULL or tibble::as_tibble(t(combn(sig_select, 2)))
sig_pairs <- tibble::as_tibble(t(combn(c(sig_select, gene_select), 2)))  
if(!is.null(sig_pairs)){
  colnames(sig_pairs) <- c("Signature_name", "Signature_other")
}

# A way to select signature pairs of particular importance
sig_pair_keyword <- "Chiji."


## ADJUST PARAMETERS ----


# Set absent signature names
if(is.null(sig_rename)){
  sig_rename <- sig_select
}

length(sig_select) == length(sig_rename)

# Create a dictionary to rename signatures
sig_dict <- c(sig_rename, gene_select)
names(sig_dict) <- c(sig_select, gene_select)


## PATHS ----


# Path to a signature file
sig_path <- "data/import/Signatures/signatures_pdac_20230720.txt"

# Location of pre-processed data
input_paths <- c("XXXX")

# Output folder
output_path_init1 <- "XXXX"
output_path_init2 <- paste0(output_path_init1, cohort_name, "/")

# Create output folders
dir.create(output_path_init1, showWarnings = FALSE)
dir.create(output_path_init2, showWarnings = FALSE)


## FIND DATA ----


# Find all files following the regex
file_list <- unname(unlist(sapply(input_paths, function(p) dir(p, pattern=paste0(cohort_regex, ".*.rds"), full.names=TRUE))))

# Identify samples user wants to exclude and remove them from the file list
if(!is.null(samples_exclude)){
  samples_exclude_indx <- c(unname(sapply(samples_exclude, function(x) grep(x, file_list))))
  file_list <- setdiff(file_list, file_list[samples_exclude_indx])
}

# Gather all correlation coefficient tibbles into a list
corr_list <- list()


## PROCESS INDIVIDUAL SAMPLES ----


for(f in file_list){
  
  
  ## INGEST DATA ----
  
  
  # Seurat data
  data_seurat <- readRDS(f)
  
  # Sample name
  sample_name <- data_seurat@misc$user.Sample_Name
  if(is.null(sample_name)){
    sample_name <- cohort_name
  }
  
  # Sample specific output path
  output_path <- paste0(output_path_init2, sample_name, "/")
  dir.create(output_path, showWarnings = FALSE)
  
  # File to write the correlation data
  filename_cor <- paste0(output_path, "cor_matrix_", sample_name, ".txt")
  
  if(!file.exists(filename_cor)){
    
    
    ## WRANGLE DATA ----
    
    
    ## Calculate signature scores and add them to Seurat meta data
    
    # Find genes abundant in this sample
    gene_list <- seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold, 
                                                 assay=assay, slot=slot)
    
    # Load signatures and filter them down to only well represented genes
    signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select,
                                                sig_length_min=1, sig_length_max=1000, ratio_threshold=0)
    
    # Add user defined one-off genes to the signature list for processing
    if(!is.null(gene_select)){
      for(i in 1:length(gene_select)){
        if(gene_select_symbol[i] %in% gene_list){
          signature_list[[gene_select[i]]] <- c(gene_select_symbol[i])
        }
      }
    }
    
    # Add signature scores to a seurat object
    data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.", assay=assay,
                                           method="ADM")
    
    # Find column names in the Seurat object - Seurat can rename original signatures
    sig_names <- colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))])
    
    
    ## EXTRACT SIGNATURE DATA ----
    
    
    # Wide tibble of signature scores
    sig_wide_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate") %>%
      dplyr::select(dplyr::all_of(c("Coordinate", sig_names))) %>%
      dplyr::mutate_if(is.numeric, round, 2)
    
    # Write signature matrix to file
    filename <- paste0(output_path, "sig_matrix_", sample_name, ".txt")
    readr::write_delim(sig_wide_df, filename, delim = "\t", append=FALSE, col_names=TRUE)
    
    
    # Create a long tibble of signature scores
    sig_df <- sig_wide_df %>%
      df_wide2long_my(key="Signature_name", val="Score")
    
    # Find signatures that reach a maximum threshold
    signatures_max <- sig_df %>%
      dplyr::group_by(Signature_name) %>%
      dplyr::summarise(ScoreMax = max(Score))
    
    # Find signature names left after thresholding by user-definied maximum intensity
    sig_max_list <- signatures_max %>%
      dplyr::filter(ScoreMax >= sig_threshold) %>%
      dplyr::pull(Signature_name)
    
    # Filter signature data to only include strong signature profiles
    sig_wide_df <- sig_wide_df[, c("Coordinate", sig_max_list)]
    
    
    ## ANALYZE SIGNATURE CORRELATIONS ----
    
    
    # Calculate gene pairwise correlations between all distinct coding sequences
    corrr_wide_df <- corrr::correlate(sig_wide_df %>% 
                                        dplyr::select(-Coordinate), quiet = FALSE) %>%
      dplyr::mutate_if(is.numeric, round, 3)
    
    # Create a long tibble of correlations
    corr_df <- corrr_wide_df %>% 
      df_wide2long_my(key="Signature_other", val="R2") %>%
      dplyr::rename(Signature_name = term) %>%
      tidyr::replace_na(list(R2 = 0)) %>%
      dplyr::mutate(Sample_Name = sample_name) %>%
      dplyr::arrange(Signature_name, Signature_other)
    
    # Gather all correlation coefficient tibbles into a list
    corr_list[[sample_name]] <- corr_df
    
    # Write signature matrix to file
    readr::write_delim(corr_df, filename_cor, delim = "\t", append=FALSE, col_names=TRUE, ".txt")
    
    
    ## FIND CORRELATED SIGNATURE GROUPS ----
    
    
    # Join signature pair properties with signature pair correlation
    # Select correlated groups of signatures of certain size
    sig_groups_df <- define_correlated_signature_groups_my(corr_df, corr_threshold=corr_threshold) %>%
      dplyr::mutate(Cluster_ID = paste0("c_", Cluster_ID)) %>%
      dplyr::group_by(Cluster_ID) %>%
      dplyr::mutate(Cluster_size = dplyr::n_distinct(Signature_name)) %>%
      dplyr::filter(Cluster_size >= cluster_size) %>%
      dplyr::ungroup()
    
    
    ## SIGNATURE CORRELATED PLOT ----
    
    
    # For each signature group, create a correlation heatmap and spatial visualization
    for(group_name in unique(sig_groups_df$Cluster_ID)){
      
      # Information on group of signatures
      # List of signatures in the group
      sig_select_list <- sig_groups_df %>%
        dplyr::filter(Cluster_ID == group_name) %>%
        dplyr::pull(Signature_name) %>%
        unique
      
      # Create a correlation heatmap of all signatures in the group
      filename <- paste0(output_path, "corr_matrix_", sample_name, "_", group_name)
      cor_subset_heatmap_my(corrr_wide_df, sig_select_list, filename)
      
    }
  } else{
    
    # Read the sample-level correlation data
    corr_list[[sample_name]] <- readr::read_delim(filename_cor, delim="\t")
    
  }
  
  gc()
} 


## PROCESS SAMPLE COHORT ----


# Create additional plotting for signature pair correlations across multiple samples
if(!is.null(sig_pairs)){
  
  # Gather all generated correlation coefficients
  corr_join_df <- dplyr::bind_rows(corr_list) %>%
    dplyr::mutate(Signature_name = gsub("sig.", "", Signature_name),
                  Signature_other = gsub("sig.", "", Signature_other)) %>%
    dplyr::inner_join(sig_pairs, by=c("Signature_name", "Signature_other"))
  
  # Rename signatures if user provided alternative names
  corr_join_df[["Signature_name"]] <- unname(sig_dict[corr_join_df[["Signature_name"]]])
  corr_join_df[["Signature_other"]] <- unname(sig_dict[corr_join_df[["Signature_other"]]])
  
  # Create signature pair names
  corr_join_df <- corr_join_df %>%
    dplyr::mutate(Pair_name = paste(Signature_name, Signature_other, sep=" - "),
                  Keyword = grepl(sig_pair_keyword, Pair_name))
  
  # Write signature matrix to file
  filename <- paste0(output_path_init2, "cor_cohort_", cohort_name, ".txt")
  readr::write_delim(corr_join_df, filename, delim = "\t", append=FALSE, col_names=TRUE, ".txt")
  
  
  # Calculate stats for every signature pair
  corr_pair_stat_df <- corr_join_df %>%
    dplyr::group_by(Signature_name, Signature_other, Pair_name) %>%
    dplyr::summarise(RMean = round(mean(R2), 3),
                     RMedian = round(median(R2), 3),
                     RSD = round(sd(R2), 3),
                     RAbsMax = max(abs(R2)),
                     RAbove0.25 = sum(abs(R2) >= 0.25)) %>%
    dplyr::arrange(desc(abs(RMean)))
  
  # Write signature matrix to file
  filename <- paste0(output_path_init2, "cor_cohort_stat_", cohort_name, ".txt")
  readr::write_delim(corr_pair_stat_df, filename, delim = "\t", append=FALSE, col_names=TRUE, ".txt")
  
  # Select specific pairs by their correlation properties
  corr_pair_stat_df <- corr_pair_stat_df  %>%
    dplyr::filter(abs(RMean) >= 0.25 | RAbsMax >= 0.5)
  
  # Boxplot of pair correlations across samples
  for(s in unique(corr_join_df$Signature_name)){
    filename <- paste0(output_path_init2, "box_corr_cohort_", cohort_name, "_", s)
    p <- create_box_plot_my(corr_join_df %>%
                              dplyr::filter((Signature_name == s | Signature_other == s) & 
                                              Pair_name %in% corr_pair_stat_df$Pair_name), 
                            x_label="Pair_name", y_label="R2", fill_label="Keyword", 
                            labels=c("Signature pair name", "Correlation coefficient", s),
                            filename=filename, reorder_x=TRUE, outlier_shape=19, with_dots=TRUE)
    
    write_plot2file_my(p, filename, num_row=2, num_col=2)
  }
  
}
