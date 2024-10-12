# Author: Anna Lyubetskaya. Date: 20-08-06


correlate_samples_my <- function(cluster_profiles_df, filename, cluster_col="Sample_cluster"){
  
  # Calculate gene pairwise correlations between samples using DEA genes
  corr_wide_df <- corrr::correlate(cluster_profiles_df %>% 
                                     df_long2wide_my(rows="Symbol", cols=cluster_col, value="SCT_zscore_mean") %>%
                                     dplyr::select(-Symbol)) %>%
    dplyr::mutate_if(is.numeric, round, 3)
  
  # Write correlation matrix to file
  readr::write_delim(corr_wide_df, path=paste0(filename, ".txt"), delim = "\t", append=FALSE, col_names=TRUE)
  
  # Plot the correlation matrix
  correlation_plot_my(corr_wide_df, filename=filename)
  
}
