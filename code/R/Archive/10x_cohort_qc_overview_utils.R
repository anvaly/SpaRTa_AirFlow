# Author: Anna Lyubetskaya. Date: 20-03-05


measurement_qc_analysis_my <- function(data_long_df, genes_interest){
  ## Gather 10X measurements and plot QCs
  ## For each gene and for each barcode, calculate standard stats: mean, stdev, total, and =0. Plot as a histogram
  
  # Calculate stats for barcodes (mean, stdev, total)
  barcode_stat_df <- summarize_expression_metric_my(data_long_df, dimension="Barcode", metric="Counts", threshold=0,
                                                    output_file=paste0(output_folders$Matrix_QC, "stat_Barcode.txt"))
  
  # Number of barcodes
  barcode_num <- length(unique(data_long_df$Barcode))
  
  # Calculate stats for genes (mean, stdev, total)
  gene_stat_withzero_df <- summarize_expression_metric_my(data_long_df, dimension="Symbol", metric="Counts", threshold=0) %>%
    dplyr::mutate(Counts_mean_log2 = round(log2(Counts_mean + 0.5), 3),
                  Counts_sd_sqrt = round(sqrt(Counts_sd), 3),
                  IsMarker = Symbol %in% genes_interest,
                  Percent_zero = round(Counts_under_threshold / barcode_num * 100))
  
  # Cut off UMI long tail and zero values
  data_long_nozero_df <- data_long_df %>% 
    dplyr::filter(Counts > 0)
  
  # Calculate gene stats (mean, stdev, total) with zeros removed
  gene_stat_nozero_df <- summarize_expression_metric_my(data_long_nozero_df, dimension="Symbol", metric="Counts", threshold=0) %>%
    dplyr::mutate(Counts_mean_log2 = round(log2(Counts_mean + 0.5), 3),
                  Counts_sd_sqrt = round(sqrt(Counts_sd), 3))
  
  # Make columns unique for the merge
  colnames(gene_stat_nozero_df) <- paste0(colnames(gene_stat_nozero_df), "_nozero")
  
  # Merge gene statistics with and without zeros
  gene_stat_df <- gene_stat_withzero_df %>%
    dplyr::inner_join(gene_stat_nozero_df, by=c("Symbol" = "Symbol_nozero"))
  
  for(stat_type in c("", "_nozero")){
    # Create a voom plot of log2 mean of counts v sqrt standard deviation of counts for each gene
    create_scatter_plot_my(gene_stat_df, 
                           x_label=paste0("Counts_mean_log2", stat_type), y_label=paste0("Counts_sd_sqrt", stat_type), 
                           fill_label="IsMarker", shape="IsMarker", size=1, 
                           filename=paste0(output_folders$Matrix_QC, "scatter_counts_voom", stat_type), 
                           labels=c("Counts, mean, log2", "Counts, standard deviation, sqrt", 
                                    "Voom plot: log2 mean v sqrt SD of counts"))
    
    # Create a scatter plot of number of zero samples v mean UMI value for the gene
    create_scatter_plot_my(gene_stat_df, 
                           x_label=paste0("Counts_mean_log2", stat_type), y_label="Percent_zero", 
                           fill_label="IsMarker", shape="IsMarker", size=1, 
                           filename=paste0(output_folders$Matrix_QC, paste0("scatter_mean_v_dropout", stat_type)), 
                           labels=c("Counts, mean, log2", "Number of zero samples, %", "Log2 mean v number of zero measurements per gene"))
  }
  
  # Write gene stats to file
  readr::write_delim(gene_stat_df, path = paste0(output_folders$Matrix_QC, "stat_genes.txt"), delim = "\t")

  # Write target gene stats to file
  readr::write_delim(gene_stat_df %>% 
                       dplyr::filter(Symbol %in% genes_interest), 
                     path = paste0(output_folders$Matrix_QC, "stat_markers.txt"), delim = "\t")

  return(gene_stat_df)
}
