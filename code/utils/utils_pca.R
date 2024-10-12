# Author: Anna Lyubetskaya. Date: 19-10-01
# Perform PCA

# todo: impute missing values before PCA
# UMAP: https://cran.r-project.org/web/packages/umap/vignettes/umap.html


source("code/utils/utils_ggplot.R")


pca_from_long_tibble_my <- function(data_df, meta_df, params, filename="~/pca", table2file=FALSE, 
                                    extra_plot=FALSE, dot_labels="Sample_ID"){
  ## Perform PCA analysis on RNA-seq data for selected genes and cell lines
  ## Row = samples, columns = features
  
  ## params <- list(sample_value = "DepMap_ID", sample_filter = c(),
  ##                feature_value = "Hugo_Symbol", feature_filter = c(),
  ##                cell_value = "Expression", color_by = c())
  
  # Select data for PCA
  expr_filt_df <- data_df %>% 
    dplyr::filter(!!rlang::sym(params$sample_value) %in% params$sample_filter & 
                    !!rlang::sym(params$feature_value) %in% params$feature_filter)
  
  # Go from long to wide format
  expr_filt_wide_df <- df_long2wide_my(expr_filt_df, rows=params$sample_value, cols=params$feature_value, value=params$cell_value)
  
  # Write the wide matrix to file
  if(table2file == TRUE){
    readr::write_delim(expr_filt_wide_df, paste0(filename, "_expr.txt"), delim = "\t")
  }
  
  # Run PCA analysis
  pca_analysis_my(expr_filt_wide_df, meta_df, sample=params$sample_value, feature=params$feature_value, 
                  col_groups=params$color_by, filename=filename, extra_plot=extra_plot, dot_labels=dot_labels)
}


pca_analysis_my <- function(df_wide, meta_df, sample="DepMap_ID", feature="Symbol", col_groups=c(), filename="pca", 
                            extra_plot=FALSE, dot_labels="Sample_ID"){
  ## Perform PCA analysis
  ## Don't forget to prefilter your df before calling this function
  ## https://tbradley1013.github.io/2018/02/01/pca-in-a-tidy-verse-framework/
  ## https://cmdlinetips.com/2019/05/how-to-do-pca-in-tidyverse-framework/
  
  # Remove input columns that have NA: PCA requires all values to be present. Consider imputing instead
  df_wide <- df_wide[, colSums(is.na(df_wide)) == 0]
  
  # Perform PCA analysis. Make sure to remove non-numerical columns
  pca_df <- df_wide %>% 
    dplyr::select(-dplyr::all_of(sample)) %>% 
    tidyr::nest(data = everything()) %>% 
    dplyr::mutate(pca = purrr::map(data, ~ prcomp(.x, center = TRUE, scale = TRUE))) %>% 
    dplyr::mutate(pca_aug = purrr::map2(pca, data, ~broom::augment(.x, data = .y)))
  
  # Find and plot percent variance explained by component
  pca2plot_var_df <- pca_df$pca %>% 
    purrr::map(~broom::tidy(.x, data = .y, "pcs")) %>% 
    as.data.frame() %>%
    dplyr::mutate(percent = round(percent * 100, 2),
                  High = percent >= 5)
  
  # Plot variance explained by PC
  if(extra_plot == TRUE){
    create_bar_plot_my(pca2plot_var_df %>% 
                         dplyr::filter(percent >= 1), 
                       x_label="PC", y_label="percent", fill_label="High", 
                       filename=paste0(filename, "_perc_var"), 
                       labels=c("PC", "% variance", "Percent variance explained by components"))
  }
  
  # Extract variance captured by first 2 components
  pca2plot_var <- pca2plot_var_df %>% 
    dplyr::filter(PC <= 2) %>% 
    dplyr::pull(percent) %>% 
    round()
  
  # Extract full PCA information relative to original variables
  pca2plot_df <- pca_df %>% 
    tidyr::unnest(pca_aug)
  
  # Backfill sensible rownames - they get removed during PCA execution
  pca2plot_df[sample] <- df_wide[sample]
  
  # Join PCA data with meta data
  pca2plot_df <- pca2plot_df %>% 
    dplyr::left_join(meta_df, by=sample)
  
  # First two PCs with meta data
  pc2file <- pca2plot_df[c(colnames(meta_df), ".fittedPC1", ".fittedPC2")]
  readr::write_delim(pc2file, paste0(filename, ".txt"), delim = "\t")
  
  for(gr in col_groups){
    x_label <- paste0("PC 1=", pca2plot_var[1], "%")
    y_label <- paste0("PC 2=", pca2plot_var[2], "%")
    
    title <- paste0("PCA of ", sample, " using ", feature, " colored by ", gr)
    
    # Plot first 2 components
    p <- create_scatter_plot_my(pca2plot_df, x_label=".fittedPC1", y_label=".fittedPC2", 
                                fill_label=gr, size=2, filename=NULL, 
                                labels=c(x_label, y_label, title), do_fit=NULL, 
                                dot_labels=dot_labels, cols=define_cols_my(n=length(unique(pca2plot_df[[gr]]))))
    
    write_plot2file_my(p, paste0(filename, "_", gr), num_row=1, num_col=2)
    
  }
  
  return(pc2file)
}
