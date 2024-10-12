# Author: Anna Lyubetskaya. Date: 20-02-10


source("code/utils/utils_biomart.R")
source("code/utils/utils_tibble.R")


calculate_tpm_my <- function(counts, lengths){
  ## Calculate TPM from a wide matrix of gene counts
  ## https://support.bioconductor.org/p/91218/
  
  x <- counts / lengths * 1e6
  mat <- t(t(x) * 1e6 / colSums(x))
  
  return(x)
}


tibble_to_matrix <- function(wide_df, field="Symbol"){
  ## Tranform wide expression tibble to a matrix for CPM and TPM calculations
  
  data_matrix <- wide_df %>%
    dplyr::arrange(!!sym(field)) %>%
    tibble::column_to_rownames(field) %>%
    data.matrix()
  
  return(data_matrix)
}


cpm_from_matrix_to_tibble_my <- function(data_matrix, field="Symbol"){
  ## Calculate CPM or TPM from wide expression matrix, make resulting tibble long, and add log value

  # Calculate CPMs from raw counts and report a long tibble
  cpm_df <- data_matrix %>%
    edgeR::cpm() %>%
    tibble::as_tibble(rownames="id") %>%
    dplyr::rename(field = "id") %>%
    df_wide2long_my(key=sample_field, val="CPM") %>%
    dplyr::mutate(log2CPM = round(log2(CPM+0.001), 2)) %>%
    dplyr::mutate(CPM = round(CPM, 2))
  
  return(cpm_df)
}


tpm_from_matrix_to_tibble <- function(data_matrix, gene_lengths_df, field="Symbol"){
  ## Calculate TPMs from raw counts and biomaRt gene lengths and report a long tibble
  ## It's important that expression data and gene length data is sorted by Symbol and contains the same number of rows
  ## https://support.bioconductor.org/p/91218/
  
  # Check that measurement data and gene length data vectors are for the same genes in the same order
  gene_lengths_df <- gene_lengths_df %>%
    dplyr::arrange(!!sym(field))
  
  # Compare gene symbols in the measurement matrix and in the gene length tibble
  gene_match_vector <- rownames(data_matrix) == gene_lengths_df %>% 
    dplyr::pull(!!sym(field))
  
  if(!"FALSE" %in% names(table(gene_match_vector))){
    # Calculate TPMs from raw counts and biomaRt gene lengths and report a long tibble
    tpm_df <- data_matrix %>%
      calculate_tpm_my(gene_lengths_df$transcript_length_max) %>%
      tibble::as_tibble(rownames="id") %>%
      dplyr::rename(field = "id") %>%
      df_wide2long_my(key=sample_field, val="TPM") %>%
      dplyr::mutate(log2TPM = round(log2(TPM+0.001), 3)) %>%
      dplyr::mutate(TPM = round(TPM, 3))
  } else{
    print("ERROR: Gene length data doesn't have the same symbols as the measurement matrix")
    tpm_df <- NULL
  }
  
  return(tpm_df)
}


find_gene_lengths_from_ensembl_my <- function(genes_meta_df){
  ## Given a list of ARCHS4 Ensembl IDs, find gene length information
  ## Minimal input fields: Symbol and gene_enseblid
  
  # Extract Ensembl IDs and gene symbols
  gene_ids_df <- genes_meta_df %>%
    dplyr::select(gene_ensemblid, Symbol) %>%
    dplyr::rename("ensembl_gene_id" = "gene_ensemblid")
  
  # Extract all Ensembl IDs
  gene_list <- gene_ids_df$ensembl_gene_id %>% unique
  
  # Translate ARCHS4 Ensembl IDs to RefSeq transcript IDs for compatibility with the CRISPR library
  transcript_lengths <- ensembl2length_my(gene_list, dataset="mmusculus_gene_ensembl")
  
  # Report genes that were not annotated by biomaRt
  genes_not_found <- setdiff(gene_list, transcript_lengths$ensembl_gene_id %>% unique)
  cat("biomaRt failed to find gene lengths for ", length(genes_not_found), " genes of ", length(gene_list), "\n")
  # cat(genes_not_found)
  
  # For each gene Symbol find the longest transcript
  gene_lengths <- transcript_lengths %>%
    dplyr::group_by(ensembl_gene_id) %>%
    dplyr::summarise(transcript_length_max = max(transcript_length)) %>%
    dplyr::inner_join(gene_ids_df, by="ensembl_gene_id") %>%
    dplyr::arrange(Symbol)
  
  return(genes_meta_df %>%
           dplyr::left_join(gene_lengths, by="Symbol"))
}


filter_no_annotation_my <- function(data_df){
  ## Remove genes with no annotations
  ## Input is a wide matrix with the first column called "Symbol"
  
  # Remove symbols like "Gm\d+" (predicted genes with no annotations)
  genes_gm <- stringr::str_detect(data_df$Symbol %>% unique, stringr::regex("Gm\\d+"))
  
  # Remove weak predicted genes from original data
  if(!is.null(genes_gm)){
    data_df <- data_df %>%
      dplyr::filter(!Symbol %in% data_df$Symbol[genes_gm])
  }
  
  return(data_df)  
}


filter_mean_reads_my <- function(data_df, field="Counts", threshold=0){
  ## Remove genes with no reads
  ## Input is a wide matrix with the first column called "Symbol"
  
  # Detect genes with no reads at all
  genes_filtered <- data_df %>%
    dplyr::group_by(Symbol) %>%
    dplyr::summarise(mean_reads = mean(!!rlang::sym(field))) %>%
    dplyr::filter(mean_reads <= threshold)
  
  # Remove filtered genes from original data
  if(!is.null(genes_filtered)){
    data_df <- data_df %>%
      dplyr::filter(!Symbol %in% genes_filtered$Symbol)
  }
  
  return(data_df)
}


counts_cpm_tpm_merge_filter_my <- function(expression_df, sample_field="GEO_ID", sample_meta_df=NULL, gene_lengths_df=NULL, gene_list=NULL, 
                                           mean_threshold=-100, output_folder=NULL){
  ## From a wide matrix of counts, calculate CPM, log2CPM, and log2TPM values
  ## Merge all data into a long matrix
  
  ## Column for gene symbols = Symbol
  ## Column for samples = GEO_ID
  ## Measurements included in the final table: Counts, CPM, log2CPM, log2TPM
  
  ## Inputs:
  ## sample_meta_df = tibble containing meta data (joined to measurements using sample ID = GEO_ID)
  ## gene_list = a list of genes to filter measurements by (using gene ID = Symbols)
  ## gene_lengths_df = tibble containing gene lengths (columns = [Symbol, transcript_length_max], e.g. from function find_gene_lengths_from_ensembl_my())
  
  if(is.null(gene_list)){
    gene_list <- expression_df$Symbol %>% 
      unique
  }
  
  # Tranform wide expression tibble to a matrix for CPM and TPM calculations
  expression_matrix <- tibble_to_matrix(expression_df, field="Symbol")

  # Long tibble of raw counts
  counts_df <- expression_df %>% 
    df_wide2long_my(key=sample_field, val="Counts")
  
  # Calculate CPMs from raw counts and report a long tibble
  cpm_df <- cpm_from_matrix_to_tibble(expression_matrix)
  
  # Create a long table of counts and CPMs with cell line annotation
  expression_long_df <- cpm_df %>%
    dplyr::inner_join(counts_df, by=c("Symbol", sample_field))
  
  # Calculate TPMs from raw counts and biomaRt gene lengths and report a long tibble
  tpm_df <- tpm_from_matrix_to_tibble(data_matrix, gene_lengths_df, field="Symbol")

  if(!is.null(tpm_df)){
    # Add TPMs to the long tibble
    expression_long_df <- expression_long_df %>%
      dplyr::inner_join(tpm_df, by=c("Symbol", sample_field))
  }
  
  if(!is.null(sample_meta_df)){
    expression_long_df <- expression_long_df %>%
      dplyr::inner_join(sample_meta_df, by=sample_field)
  }
  
  # Write full data tables to file
  if(!is.null(output_folder)){
    readr::write_delim(expression_long_df, path = paste0(output_folder, "data.txt"), delim = "\t")
  }
  
  # Remove genes not present in the reference
  # Remove genes with no annotations and no reads
  expression_long_filt_df <- expression_long_df %>%
    dplyr::filter(tolower(Symbol) %in% tolower(gene_list)) %>%
    filter_no_annotation_my() %>%
    filter_mean_reads_my(threshold = mean_threshold)
  
  # Write pre-filtered data to file
  if(!is.null(output_folder)){
    readr::write_delim(expression_long_filt_df, path = paste0(output_folder, "data_filt.txt"), delim = "\t")
  }
  
  # Total number of genes
  gene_total_number <- expression_long_filt_df$Symbol %>% unique %>% length()
  cat("Total number of genes = ", gene_total_number, "\n")
  
  # Total number of samples
  sample_total_number <- expression_long_filt_df$GEO_ID %>% unique %>% length()
  cat("Total number of samples = ", sample_total_number, "\n")
  
  return(expression_long_filt_df)
}


summarize_expression_metric_my <- function(expression_df, dimension="Symbol", metric="Counts", threshold=1, output_file=NULL){
  ## From a long matrix containing Counts, log2CPM, and TPM, calculate various statistics for a given dimension
  # Suggested thresholds for various metrics: log2TPM = 1, log2CPM = 0.5
  
  col_list <- colnames(expression_df)
  
  data <- list()
  if(metric %in% col_list){
    # Calculate basic stats for each sample
    stat_df <- expression_df %>%
      dplyr::group_by_at(dplyr::vars(dimension)) %>%
      dplyr::summarise(metric_total = round(sum(!!sym(metric)), 3),
                       metric_mean = round(mean(!!sym(metric)), 3),
                       metric_sd = round(sd(!!sym(metric)), 3),
                       metric_over_threshold = sum(!!sym(metric) > threshold),
                       metric_under_threshold = sum(!!sym(metric) <= threshold))
    
    # Make column names specific to the metric
    colnames(stat_df) <- gsub("metric", metric, colnames(stat_df))
  }
  
  # Print sample stats to file
  if(!is.null(output_file)){
    readr::write_delim(stat_df, path = output_file, delim = "\t")
  }
  
  return(stat_df)
}
