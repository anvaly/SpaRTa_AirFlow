# Author: Anna Lyubetskaya. Date: 21-02-01


# To enable %dopar%
library(foreach)

source("code/utils/utils_signatures.R")
source("code/R/Utils/utils_10X_matrix.R")

if(!"UCell" %in% rownames(installed.packages())){
  BiocManager::install("UCell")
}

signature_empirical_pvalue_my <- function(data_seurat, signature_list, gene_list, output_path, sample_name, sig_name, sig_random_filename=NULL,
                                          n_simulations=1000, assay="SCT", output_col_name="Sig_Emp_Pvalue", save2file=FALSE){
  ## Create a random signature background distribution and score the target signature against it
  
  
  # Write bootstrapped gene expression scores to this file
  if(is.null(sig_random_filename)){
    sig_random_filename <- paste0(output_path, "/sig_random_", sample_name, "_", sig_name, ".txt")
  }
  
  # Make sure the object has the user.Sample_Name column in the meta.data
  if(!"user.Sample_Name" %in% colnames(data_seurat@meta.data)){
    data_seurat@meta.data[["user.Sample_Name"]] <- data_seurat@misc$user.Sample_Name
  }
  
  # Read from file or generate random gene expression scores matching target signature
  if(!file.exists(sig_random_filename)){
    
    # Generate random signature statistics for each sample in the cohort
    sig_random_stat_list <- list()
    for(s in unique(data_seurat@meta.data[["user.Sample_Name"]])){
      # Subset Seurat object by sample
      barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data[["user.Sample_Name"]] == s),])
      data_seurat_subset <- subset(data_seurat, cells=barcode_list)
      
      # Calculate random gene signature
      sig_random_stat_list[[s]] <- calculate_random_signature_my(data_seurat_subset, gene_list, sig_length=length(signature_list), n_simulations=n_simulations, assay=assay)
    }
    
    # Stats from random signatures of length matching the selected target signature
    sig_random_stat_df <- dplyr::bind_rows(sig_random_stat_list) %>%
      tibble::rownames_to_column("Coordinate")
    
    if(save2file == TRUE){
      # Write to file random gene expression scores
      readr::write_delim(sig_random_stat_df, sig_random_filename, delim="\t")
    }
    
  } else{
    sig_random_stat_df <- readr::read_delim(sig_random_filename, delim="\t")
  }
  
  # Move Coordinate to rownames for the next comparison
  sig_random_stat_df <- sig_random_stat_df %>%
    tibble::column_to_rownames("Coordinate")
  
  # Target signature scores tibble
  meta_loc_df <- data_seurat@meta.data %>%
    dplyr::select(dplyr::all_of(sig_name))
  
  # Calculate empirical p-value of the target signature being significant in each spot compared to the bootstrapped background distribution of random gene expression scores
  sig_random_num <- ncol(sig_random_stat_df)
  
  # Parallelize the empirical p-value calculation for efficiency
  cl <- parallel::makeCluster(parallel::detectCores())
  doParallel::registerDoParallel(cl)
  sig_pv_emp_list <- foreach::foreach(x = rownames(meta_loc_df), .combine=c) %dopar% {
    sum(sig_random_stat_df[x,] >= meta_loc_df[x, sig_name]) / sig_random_num
  };
  parallel::stopCluster(cl)
  
  # Add target signature significance empirical p-values to the meta data and select spots with high signature using the thresholds
  data_seurat@meta.data[output_col_name] <- sig_pv_emp_list
  
  return(data_seurat)
}


calculate_random_signature_my <- function(data_seurat, gene_list, sig_length=10, n_simulations=100, assay="SCT"){
  ## Create a random set of genes of pre-defined length and calculate the corresponding signature score
  
  # Generate X random gene lists of list Y; X and Y - function inputs
  signature_random_list <- lapply(1:n_simulations, function(x) sample(gene_list, sig_length))
  names(signature_random_list) <- paste0("random.", 1:n_simulations)
  
  # Add signature scores to a seurat object
  data_seurat <- add_signature_scores_my(data_seurat, signature_random_list, prefix="sig.", assay=assay)
  
  # Extract the randomly generated signature scores
  signature_random_df <- round(data_seurat@meta.data[paste0("sig.random.", 1:n_simulations)], 3)
  
  # Calculate stats
  # sig_mean_list <- sapply(1:ncol(signature_random_df), function(x) mean(signature_random_df[[x]]))
  # sig_sd_list <- sapply(1:ncol(signature_random_df), function(x) sd(signature_random_df[[x]]))
  # stats_list <- c(mean(sig_mean_list), sd(sig_sd_list))
  
  return(signature_random_df)
}


add_signature_scores_my <- function(data_seurat, signature_list, prefix="sig.", assay="SCT", cleanup=TRUE, method="UCell"){
  ## Add a list of signature scores to a Seurat object
  ## Methods: AddModuleScore (AMS) and UCell
  
  # Remove any old signature information
  if(cleanup == TRUE){
    data_seurat@meta.data <- data_seurat@meta.data[which(!grepl("sig.", colnames(data_seurat@meta.data)))]
  }
  
  # Column names originally
  col_names <- colnames(data_seurat@meta.data)
  
  # Remove signatures that are already part of the object
  signature_list <- signature_list[setdiff(names(signature_list), gsub(prefix, "", col_names))]
  
  if(length(signature_list) > 0){
    # Calculate signature score and add it to the seurat object  
    if(method == "AMS"){
      data_seurat <- Seurat::AddModuleScore(data_seurat, features=signature_list, name="Signature", assay=assay, nbin=10)
    } else{
      data_seurat <- UCell::AddModuleScore_UCell(data_seurat, features=signature_list, name="Signature", assay=assay)
    }
    
    # Adjust file names
    colnames(data_seurat@meta.data) <- c(col_names, paste0(prefix, names(signature_list)))
  }
  
  return(data_seurat)
}


annotate_seurat_with_signatures <- function(data_seurat, sig_path, sct_threshold=0.5, spot_threshold=10, 
                                            sig_select=NULL, assay="SCT"){
  ## Calculate signature scores and add them to Seurat meta data
  
  # Find genes abundant in this sample
  gene_list <- seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold)
  
  # Load signatures and filter them down to only well represented genes
  signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select,
                                              sig_length_min=1, sig_length_max=1000, ratio_threshold=0)
  
  # Add signature scores to a seurat object
  data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.",
                                         assay=assay)
  
  return(data_seurat)
}


calculate_signature_stats_my <- function(data_seurat, sig_names){
  ## Calculate signature stats
  
  # Wide tibble of signature scores + nCount_Spatial and seurat_clusters
  sig_wide_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate") %>%
    dplyr::select(c("Coordinate", dplyr::all_of(sig_names)))
  
  # Create a long tibble of signature scores
  sig_df <- sig_wide_df %>%
    df_wide2long_my(key="Signature_name", val="Score")
  
  # Calculate various signature stats
  sig_stat_df <- sig_df %>%
    dplyr::group_by(Signature_name) %>%
    dplyr::summarise(score_min = min(Score),
                     score_max = max(Score),
                     score_delta = abs(max(Score)) - abs(min(Score)),
                     score_mean = mean(Score),
                     score_sd = sd(Score),
                     num_above_threshold = sum(Score >= score_mean + score_sd)) %>%
    dplyr::mutate_if(is.numeric, round, 2)
  
  return(sig_stat_df)
}
