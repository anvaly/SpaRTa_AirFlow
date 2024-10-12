# Author: Anna Lyubetskaya. Date: 20-05-11


source("code/utils/utils_specialized_plots.R")


find_sig_groups <- function(sig_corr_df, corr_threshold, sig_grep){
  ## Select correlated signature pairs using thresholds and find signature clusters
  
  # Find all positively correlated signatures with non intersecting gene lists
  sig_corr_select_df <- sig_corr_df %>%
    dplyr::filter(R2 >= corr_threshold & 
                    grepl(sig_grep, Signature_name) & 
                    grepl(sig_grep, Signature_other))
  
  # Find all signatures that belong to the same cluster
  sig_graph <- igraph::graph_from_data_frame(data.frame(sig_corr_select_df[, c("Signature_name", "Signature_other")]),
                                             directed=FALSE, vertices=NULL)
  
  # Signature group lists
  sig_groups_list <- igraph::clusters(sig_graph)$membership
  
  return(sig_groups_list)
}


define_correlated_signature_groups_my <- function(sig_corr_df, corr_threshold=0.5, sig_grep="sig."){
  ## Using a matrix of pairwise correlations, define cliques of (anti)correlated signatures
  
  # Select positively correlated signature pairs using thresholds and find signature clusters
  sig_groups_pos_list <- find_sig_groups(sig_corr_df, corr_threshold, sig_grep)

  sig_groups_pos_df <- NULL
  if(length(sig_groups_pos_list) > 0){
    # Process signature clusters as a tibble
    sig_groups_pos_df <- tibble::tibble(Signature_name = names(sig_groups_pos_list),
                                        Cluster_ID = paste0("pos", unname(sig_groups_pos_list)))
  }

  # # Select negatively correlated signature pairs using thresholds and find signature clusters
  # sig_groups_neg_list <- find_sig_groups(sig_corr_df %>%
  #                                          dplyr::mutate(R2 = R2*-1), corr_threshold, sig_grep)
  # 
  # sig_groups_neg_df <- NULL
  # if(length(sig_groups_neg_list) > 0){
  #   # Process signature clusters as a tibble
  #   sig_groups_neg_df <- tibble::tibble(Signature_name = names(sig_groups_neg_list),
  #                                       Cluster_ID = paste0("neg", unname(sig_groups_neg_list)))
  # }

  # return(rbind(sig_groups_pos_df, sig_groups_neg_df))
  return(sig_groups_pos_df)
}


cor_subset_heatmap_my <- function(corrr_wide_df, sig_list, filename, scale=c(-1,1,0)){
  ## Subset correlations and visualize
  ## Link: http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2
  ## Input is a wide square matrix of correlations with the first columns named "term"
  
  # Select appropriate data
  corrr_wide_select_df <- corrr_wide_df %>%
    dplyr::filter(term %in% sig_list) %>%
    dplyr::select(term, sort(sig_list)) %>%
    dplyr::arrange(term)
  
  # Write correlation matrix to file
  readr::write_delim(corrr_wide_select_df, 
                     paste0(filename, ".txt"), 
                     delim = "\t", append=FALSE, col_names=TRUE)

  # Create a clustregram of correlations
  p <- correlation_plot_my(corrr_wide_select_df, scale=scale, cols=NULL, filename=filename)
}
