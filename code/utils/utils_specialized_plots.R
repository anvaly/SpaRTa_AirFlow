# Author: Anna Lyubetskaya. Date: 20-05-01


if(!"ggcorrplot" %in% rownames(installed.packages())){
  install.packages("ggcorrplot")
}

source("code/utils/utils_ggplot.R")


scatter2boxplot_my <- function(data_df, x_label, y_label, filename, fill_label=NULL, facet_var=NULL, round_up=FALSE){
  ## For two numerical vectors, bin X-axis values and create a series of boxplots
  
  col_names <- c("Category", x_label, y_label)
  
  if(!is.null(facet_var)){
    col_names <- c(col_names, facet_var[1])
  }
  
  if(!is.null(fill_label)){
    col_names <- c(col_names, fill_label)
  }
  
  # Split X-axis into interval categories  
  if(round_up == TRUE){
    data_loc_df <- data_df %>%
      dplyr::mutate(Category = as.character(ceiling(!!sym(x_label) * 2) / 2))
  } else{
    data_loc_df <- data_df %>%
      dplyr::mutate(Category = as.character(round(!!sym(x_label))))
  }
  
  if(is.null(fill_label)){
    fill_label <- "Category"
  }
  
  # Split X-axis into interval categories
  data_loc_df <- data_loc_df %>%
    dplyr::arrange(!!sym(x_label)) %>%
    dplyr::select(dplyr::all_of(col_names)) %>%
    tidyr::drop_na()
  
  # Enforce the order the X axis
  data_loc_df$Category <- factor(data_loc_df$Category, levels=unique(data_loc_df$Category))
  
  # Box plot to file
  title <- paste0(x_label, " v ", y_label)
  
  p <- create_box_plot_my(data_loc_df, x_label="Category", y_label=y_label, 
                          fill_label=fill_label, facet_var = facet_var,
                          filename=filename, labels=c(x_label, y_label, title), with_dots=TRUE)
  
  return(p)
}


survival_plot_my <- function(data_df, time="OS", event="OS.CNSR", factors=c("Category"), 
                             covariates=c("Trial_ID", "TumorType_Organ", "TreatmentStandard"), label="", output_folder=output_folder){
  ## Perform survival analysis
  
  # Stop analysis if one of the inputs is not present in the dataframe
  stopifnot(all(c(time, event, factors, covariates) %in% colnames(data_df)))
  
  # Create a fit formula for the survival plot
  formula_my <- as.formula(paste0("survival::Surv(", time, ", ", event, ")", " ~ ", paste0(factors, collapse=" + "), " - ", paste0(covariates, collapse=" - ")))
  
  # Compute the predicted survivor function for a Cox proportional hazards model
  # Analysis types survfit v coxph
  fit <- eval(substitute(survival::coxph(formula_my, data = data_df), list(formula_my = formula_my)))
  
  # Don't draw confidence intervals if one of the groups has only 1 sample
  if(length(fit$n) >= 1){
    if(min(fit$n) == 1){
      conf_int <- FALSE
    } else{
      conf_int <- TRUE
    }
    
    # Plot a survival curve
    survival_plot <- survminer::ggsurvplot(fit, data = data_df,
                                           conf.int = conf_int, conf.int.alpha = 0.05, pval = TRUE,
                                           risk.table = TRUE, risk.table.height=0.3, risk.table.y.text=FALSE,
                                           legend = "right", title = label)  # palette="Paired"
    
    # Save a survival curve to file
    filename <- paste0(output_folder, paste("/survival", label, time, 
                                            paste0(gsub("-|_", "", factors), collapse="-"), 
                                            paste0(gsub("-|_", "", covariates), collapse="-"), sep="_"))
    write_plot2file_my(print(survival_plot), filename, num_row=2, num_col=2)
  }
}


correlation_plot_my <- function(corr_wide_df, scale=NULL, cols=NULL, filename=NULL, rowname_col="term", hc.order=TRUE){
  ## Create a clustregram of correlations
  ## scale = c(min, max, midpoint)
  
  # Set a colormap if not provided
  if(is.null(cols)){
    cols <- c("darkblue", "darkred", "white")
  }
  
  # Create a correlation heatmap
  p <- ggcorrplot::ggcorrplot(corr_wide_df %>%
                                tibble::column_to_rownames(rowname_col),
                              hc.order = hc.order, outline.col = "white") +
    theme(axis.text.y = element_text(angle=0, hjust=1, size=8),
          axis.text.x = element_text(angle=90, hjust=1, size=8))
  
  # Adjust the scale
  if(!is.null(scale)){
    p <- p + 
      scale_fill_gradient2(midpoint = scale[3], limit = c(scale[1], scale[2]), 
                           low = cols[1], mid = cols[3], high = cols[2])
  }
  
  if(!is.null(filename)){
    num_row <- ceiling(nrow(corr_wide_df) / 30)
    write_plot2file_my(p, filename, num_row=num_row, num_col=num_row)
  }
  
  return(p)
}
