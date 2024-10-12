# Author: Anna Lyubetskaya. Date: 19-11-22
# todo: FIX FORMULA


# Installation for Domino instances
if(!"spatialEco" %in% rownames(installed.packages())){
  install.packages("spatialEco")
}
if(!"ggcorrplot" %in% rownames(installed.packages())){
  install.packages("ggcorrplot")
}

library(ggplot2)
source("code/utils/utils_ggplot.R")


find_density_modes_my <- function(x){
  ## Find modes of a density distribution
  ## https://stackoverflow.com/questions/27418461/calculate-the-modes-in-a-multimodal-distribution-in-r
  
  modes <- NULL
  for (i in 2:(length(x)-1)){
    if ((x[i] > x[i-1]) & (x[i] > x[i+1])) {
      modes <- c(modes,i)
    }
  }
  if (length(modes) == 0){
    modes <- 'This is a monotonic distribution'
  }
  
  return(modes)
}


density_analysis_my <- function(density_vec){
  ## Analyze a vector density distribution: mode and local minima
  
  # Calculate the density distribution for color values of the low resolution image
  value_density <- density(density_vec)
  
  # Find density modes
  mode_indx <- find_density_modes_my(value_density$y)
  mode_x <- value_density$x[mode_indx]  
  mode_y <- value_density$y[mode_indx]
  
  # Find the local minima of the color value density distribution
  min_indx <- which(value_density$y %in% spatialEco::local.min.max(value_density$y, plot=FALSE)$minima)
  min_x <- value_density$x[min_indx]
  min_y <- value_density$x[min_indx]
  
  return(list(min_x = min_x, 
              min_y = min_y,
              mode_x = mode_x, 
              mode_y = mode_y))
}


correlation_calculate_plot_my <- function(data_wide_df, row_col, output_file=NULL){
  ## Calculate gene pairwise correlations between all distinct coding sequences
  ## Write correlations to file and plot a correlation heatmap
  
  ## Input is a wide square matrix of correlations with the first columns named "rowname"
  ## Link: http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2
  
  # Prepare data for pairwise correlation 
  data_wide_df <- data_wide_df %>% 
    dplyr::select(-!!rlang::sym(row_col))
  
  # Calculate gene pairwise correlations between all distinct coding sequences
  corrr_wide_df <- corrr::correlate(data_wide_df, quiet = FALSE)
  
  if(!is.null(output_file)){
    # Write correlation matrix to file
    readr::write_delim(corrr_wide_df %>% 
                         dplyr::mutate_if(is.numeric, round, 3), 
                       paste0(output_file, ".txt"), delim = "\t", append=FALSE, col_names=TRUE)
  }
  
  return(corrr_wide_df)
}


wilcoxon_from_tibble_my <- function(df, feature, group, measurement){
  ## Perform Wilcoxon test on a long tibble
  ## Individual Wilcoxon test is performed on every feature, between two groups of the input measurement
  ## For example: feature = Hugo_Symbol, group = user-defined groups of cell lines, measurement = CERES
  ## https://www.rdocumentation.org/packages/stats/versions/3.6.1/topics/wilcox.test
  
  # Formula to perform Wilcoxon
  wilcox_formula <- as.formula(paste0(measurement, " ~ ", group))
  significance_threshold <- 0.05
  fdr_method <- "bonferroni"
  
  # Perform Wilcoxon test - formula variable not working!
  wilcoxon_df <- df %>%
    dplyr::select(dplyr::all_of(c(feature, group, measurement))) %>%  # Select only relevant data: feature, group, and measurement
    dplyr::group_by_at(dplyr::all_of(c(feature, group))) %>%  # Group by feature and group
    dplyr::rename(Variable = !!rlang::sym(group), Measurement = !!sym(measurement)) %>%
    dplyr::mutate(mean = mean(Measurement)) %>%  # Calculate mean of the metric
    dplyr::group_by(!!sym(feature)) %>%  # Group by gene
    dplyr::summarize(wilcox_pval = wilcox.test(Measurement ~ Variable)$p.value,
                     stat = wilcox.test(Measurement ~ Variable)$statistic,
                     mean_delta = round(abs(max(mean) - min(mean)), 2),
                     mean_min = round(min(mean), 2)) %>%  # Calculate Wilcoxon p-value, delta between group means, and the smallest mean
    dplyr::mutate(wilcox_padj_neglog10 = round(-log10(p.adjust(wilcox_pval, method=fdr_method)), 2)) %>%
    dplyr::mutate(IsSignificant = wilcox_pval <= significance_threshold)
  
  return(wilcoxon_df)
}


lm_residual_outliers <- function(data_df, var1, var2, name="", sd_num=3){
  ## Calculate linear regression and residuals, identify outliers
  
  
  # Calculate and display the regression line
  regression <- lm(data_df[[var1]] ~ data_df[[var2]])
  
  # Calculate residuals
  data_df[[paste0("Residuals", name)]] <- round(residuals(regression), 3)
  
  # Choose a threshhold
  outlier_threshold <- mean(residuals(regression)) + sd(residuals(regression)) * sd_num
  
  cat(name, mean(residuals(regression)), sd(residuals(regression)), c(outlier_threshold), "\n")
  
  # Print only names of outliers
  data_df[[paste0("IsOutlier", name)]] <- sapply(data_df[[paste0("Residuals", name)]], function(x) 
    if(x >= outlier_threshold){ paste0(var1, "_Outlier"); }
    else if(x <= -outlier_threshold){ paste0(var2, "_Outlier"); }
    else{"Fit"}
    )

  return(data_df)
}
