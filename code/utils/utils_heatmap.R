# Author: Anna Lyubetskaya. Date: 19-10-16
# https://jokergoo.github.io/ComplexHeatmap-reference
# todo: Establish a better default palette


library(ComplexHeatmap)
source("code/utils/utils_tibble.R")


define_cols_my <- function(n=NULL, col_type="jet", cols=NULL){
  ## Define color palette for the plot
  ## Color types = c("random", "spectral", "rich", "rainbow", "viridis", "rgb", "jet")
  
  if(is.null(cols)){
    if(col_type == "random"){
      cols <- randomcoloR::distinctColorPalette(n, runTsne = FALSE)
    } else if(col_type == "spectral"){
      cols <- rev(RColorBrewer::brewer.pal(n, "Spectral"))
    } else if(col_type == "rich"){
      cols <- gplots::rich.colors(n, palette="temperature")
    } else if(col_type == "rainbow"){
      cols <- rainbow(n)
    } else if(col_type == "viridis"){
      cols <- viridis::plasma(n)
    } else if(col_type == "rgb"){
      cols <- colorRamps::rgb.tables(n)
    } else if(col_type == "jet"){
      col_func <- colorRampPalette(c("navy", "blue", "cyan", "darkgreen", "green", 
                                     "yellow", "orange", "brown", "red", "magenta", "purple", "gray"))
      cols <- col_func(n)
    }
  }
  
  # To visualize:
  # barplot(1:n, col=cols)
  
  return(cols)
}


pull_named_list_my <- function(df, key, value, key_list){
  ## Pull a named list where keys are one parameters and values are another
  
  # Pull a dataframe with two relevant columns
  select_df <- df %>% 
    dplyr::filter(!!rlang::sym(key) %in% key_list) %>% 
    dplyr::select(dplyr::all_of(c(key, value)))
  
  # Create a named list (dictionary), with keys as col names and values as col values
  named_list <- stats::setNames(select_df %>% dplyr::pull(!!rlang::sym(value)), 
                                select_df %>% dplyr::pull(!!rlang::sym(key)))
  
  return(named_list)
}


create_annotation_my <- function(meta_df, label, values, val_list, names, annotation_type="column", annotation_name_side="left"){
  ## Gather annotation from a tibble, assign it a color scheme, and create a heatmap annotation
  
  labels_data <- list()
  palettes <- list()
  
  for(val in values){
    # Add annotation values to a named list
    data_list <- pull_named_list_my(meta_df, label, val, val_list)
    # Sort labels explicitly to match heatmap column names
    data_list <- data_list[names]
    
    # Find unique values of the named list
    labels_data_unique <- sort(unique(data_list))
    
    if((length(labels_data_unique) > 1 && !is.numeric(labels_data_unique)) || (length(labels_data_unique) > 2 && is.numeric(labels_data_unique))){
      
      print(val)
      print(labels_data_unique)
      
      labels_data[[val]] <- data_list
      
      # Add colors to a named list of palettes
      if(length(labels_data_unique) <= 20 && !is.numeric(labels_data_unique)){
        palettes[[val]] <- stats::setNames(define_cols_my(length(labels_data_unique)), labels_data_unique)
      } else{
        palettes[[val]] <- circlize::colorRamp2(c(min(labels_data[[val]]), max(labels_data[[val]])), c("white", "red4"))
        labels_data[[val]] <- as.numeric(labels_data[[val]])
      }
      
    }
  }
  
  annotation <- ComplexHeatmap::HeatmapAnnotation(df=data.frame(labels_data), col=palettes,
                                                  which=annotation_type, annotation_name_side=annotation_name_side)
  
  return(annotation)
}


setup_parameters <- function(params, row_names, col_names){
  ## Add missing parameters
  
  if(!"row_labels" %in% names(params)){
    if(length(row_names) <= 100){
      params[["row_labels"]] <- row_names      
    } else{
      params[["row_labels"]] <- rep("", length(row_names))
    }
  }
  
  if(!"column_labels" %in% names(params)){
    if(length(col_names) <= 100){
      params[["column_labels"]] <- col_names      
    } else{
      params[["column_labels"]] <- rep("", length(col_names))
    }
  }
  
  if(!"column_km" %in% names(params)){
    params[["column_km"]] <- 1
  }
  
  if(!"row_km" %in% names(params)){
    params[["row_km"]] <- 1
  }
  
  if(!"row_split" %in% names(params)){
    params[["row_split"]] <- NULL
  }
  
  if(!"show_column_dend" %in% names(params)){
    params[["show_column_dend"]] <- FALSE
  }
  
  if(!"show_row_dend" %in% names(params)){
    params[["show_row_dend"]] <- FALSE
  }
  
  if(!"row_order" %in% names(params)){
    params[["row_order"]] <- NULL
    params[["cluster_rows"]] <- TRUE
  } else{
    params[["cluster_rows"]] <- FALSE
  }
  
  if(!"column_order" %in% names(params)){
    params[["column_order"]] <- NULL
    params[["cluster_columns"]] <- TRUE
  } else{
    params[["cluster_columns"]] <- FALSE
  }
  
  return(params)
}


create_heatmap_my <- function(df, params, row_list=NULL, col_list=NULL, col_meta_df=NULL, row_meta_df=NULL, filename="hm",
                              width=10, height=10){
  ## Create an annotated clustergram
  ## Parameter named list example
  ## params <- list(cell_value = "Expression",
  ##               row_label = "Hugo_Symbol", 
  ##               col_label = "DepMap_ID", 
  ##               distance = "pearson",
  ##               row_annotation = c("Node"),
  ##               col_annotation = c("lineage_Yuan18", "KRAS_subtype"),
  ##               range = c(-3, 0, 3),
  ##               colors = c("royalblue4", "white", "red3"))
  
  # install.packages("remotes")
  # remotes::install_github("jmw86069/jamba")
  
  # If no row and col values provided, use all values in the data frame
  if(is.null(row_list)){
    row_list <- df[[rlang::sym(params$row_label)]] %>%
      unique
  }
  if(is.null(col_list)){
    col_list <- df[[rlang::sym(params$col_label)]] %>%
      unique
  }
  if(is.null(params$range)){
    params$range <- c(min(df[[params$cell_value]]), max(df[[params$cell_value]]))
    cat(params$range)
  }
  
  # Extract data from a tibble and prepare it for a heatmap visualization
  data_matrix <- df %>%
    dplyr::filter(!!rlang::sym(params$row_label) %in% row_list & !!rlang::sym(params$col_label) %in% col_list) %>% # Select rows and columns to plot
    dplyr::select(dplyr::all_of(c(params$row_label, params$col_label, params$cell_value))) %>% # Select the 3 columns relevant to the heatmap
    dplyr::group_by_at(c(params$row_label, params$col_label)) %>% # Group by row and col labels
    dplyr::mutate(mean_val = mean(!!rlang::sym(params$cell_value))) %>% # Calculate mean value to make sure there are no duplicate rows
    dplyr::select(dplyr::all_of(c(params$row_label, params$col_label, "mean_val"))) %>% # Select the 3 columns relevant to the heatmap including mean value
    df_long2wide_my(rows=params$row_label, cols=params$col_label, value="mean_val") %>% # Go from long to wide matrix
    dplyr::mutate_if(is.numeric, list(~tidyr::replace_na(., 0))) %>%
    tibble::column_to_rownames(params$row_label) %>% # Explicitly name rows
    data.matrix() # Switch to matrix
  
  cat("Creating heatmap:", filename, "Matrix size=", nrow(data_matrix), ncol(data_matrix), "\n")
  
  # Select column and row features
  col_names <- colnames(data_matrix)
  row_names <- rownames(data_matrix)
  
  # Create heatmap color scheme
  col_fun <- circlize::colorRamp2(params$range, params$colors)
  
  # Row and column annotation place holders
  row_a <- NULL
  col_a <- NULL
  
  # Create row annotation and add to heatmap
  if(!is.null(params$row_annotation) && !is.null(row_meta_df)){
    row_a <- create_annotation_my(row_meta_df, params$row_label, params$row_annotation, row_list, row_names, 
                                  annotation_type="row", annotation_name_side="top")
  }
  
  # Create column annotation and add to heatmap
  if(!is.null(params$col_annotation) && !is.null(col_meta_df)){
    col_a <- create_annotation_my(col_meta_df, params$col_label, params$col_annotation, col_list, col_names, 
                                  annotation_type="column", annotation_name_side="left")
  }
  
  # Add missing parameters
  params <- setup_parameters(params, row_names, col_names)
  
  # Initiate a heatmp
  hm <- ComplexHeatmap::Heatmap(data_matrix, 
                                name = params$cell_value, 
                                col = col_fun, 
                                na_col = "grey", 
                                column_title = paste0(params$col_label, " ", params$cell_value), 
                                row_title = params$row_label,
                                
                                row_names_gp = gpar(fontsize = 6),
                                row_names_side = "left",
                                
                                column_labels = params$column_labels,
                                row_labels = params$row_labels,
                                show_column_dend = params$show_column_dend, 
                                show_row_dend = params$show_row_dend,
                                cluster_columns = params$cluster_columns,
                                clustering_distance_columns = params$distance,                                 
                                cluster_rows = params$cluster_rows,
                                clustering_distance_rows = params$distance,
                                column_order = params$column_order,
                                row_order = params$row_order, 
                                column_km = params$column_km,
                                row_km = params$row_km,
                                row_split = params$row_split,
                                right_annotation = row_a
                                
                                # cell_fun = jamba::cell_fun_label(params$cell_text)
                                
                                #heatmap_legend_param = list(
                                #  legend_width = unit(4, "in")
                                #)
  )
  
  if(!is.null(params$col_annotation) && !is.null(col_meta_df)){
    hm <- hm %v% col_a
  }
  
  # Definte heatmap size  
  width = width
  height = height
  
  # Write heatmap to a file
  png(filename, width = width, height = height, units = "in", res = 300)
  dev.control(displaylist="enable")
  ComplexHeatmap::draw(hm)
  dev.off()
  
  cat("Heatmap created with the following parameters\n")
  print(params)
  cat("Heatmap output is here:\n", filename, "\n")
  
  return(hm)
}
