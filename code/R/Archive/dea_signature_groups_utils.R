# Author: Anna Lyubetskaya. Date: 20-05-11


vector_categories_hi_lo_my <- function(df, col, coef=1, dn_threshold=NULL, up_threshold=NULL){
  ## Partition a vector in a tibble into hi-lo quadrants based on mean and SD of the vector
  
  df <- df %>%
    tidyr::drop_na()
  
  if(is.null(up_threshold)){
    up_threshold <- mean(df[[col]]) + sd(df[[col]]) * coef
  }
  if(is.null(dn_threshold)){
    dn_threshold <- mean(df[[col]]) - sd(df[[col]]) * coef
  }
  
  cat("Thresholds:", dn_threshold, up_threshold, "\n")
  
  category_df <- df %>% 
    dplyr::mutate(Hi = ifelse(!!sym(col) >= up_threshold, "Hi", ""),
                  Lo = ifelse(!!sym(col) <= dn_threshold, "Lo", ""),
                  Category = paste0(Hi, "_", Lo)) %>%
    dplyr::select(-Hi, -Lo) %>%
    dplyr::mutate(Category = gsub("^_|_$", "", Category),
                  Category = gsub("^$", "Med", Category))
  
  return(category_df)
}


find_neighborhood_mean_signal_my <- function(data_seurat, signature_wide_df, sig, sig_select_df, radius=1, coef=0.5){
  ## For each coordinate, find its proximal coordinates and calculate a mean of a signature in that selection of spots
  ## Radius = Size of the area
  
  
  # All physical coordinates for this sample
  physical_coordinates <- tibble::as_tibble(data_seurat@images$slice1@coordinates, rownames="rowname") %>%
    dplyr::rename(Coordinate = rowname)
  
  # Create a named vector for the signature of interest
  sig_area_list <- signature_wide_df[[sig]]
  names(sig_area_list) <- signature_wide_df$Coordinate
  
  # For each spot of interest, calculate a neighborhood average score of signature of interest
  for(i in 1:nrow(sig_select_df)){
    neighbors_list <- find_spatial_neighbors_my(physical_coordinates, sig_select_df[[i, "Coordinate"]], radius)
    sig_select_df[[i, "Sig_area"]] <- mean(sig_area_list[neighbors_list])
  }
  
  return(sig_select_df)
}


find_spatial_neighbors_my <- function(physical_coordinates, barcode, radius=1, units_image=FALSE){
  ## For a spot barcode find all barcodes of neighboring spots within a certain radius
  ## Radius is defined as the number of spots away from the target spot
  
  ## Coordinate information is located within ST Seurat structure: data_seurat@images$slice1@coordinates
  ## Zero-centric coordinates: "row" and "col"; [0; Inf]
  ## Image coordinates: "imagerow" and "imagecol"
  
  
  # Find coordinates of the target barcode
  center_coordinate <- physical_coordinates %>%
    dplyr::filter(Coordinate == barcode)
  
  # Define X, Y coordinate columns
  row <- "row"
  col <- "col"
  if(units_image == TRUE){
    row <- "imagerow"
    col <- "imagecol"
  }
  
  # Define neighborhood coordinate intervals
  row_start <- center_coordinate[[row]] - radius
  row_stop <- center_coordinate[[row]] + radius
  col_start <- center_coordinate[[col]] - radius
  col_stop <- center_coordinate[[col]] + radius
  
  # Find all barcodes in the neighborhood
  neighbors_list <- physical_coordinates %>%
    dplyr::filter(!!sym(row) <= row_stop & 
                    !!sym(row) >= row_start & 
                    !!sym(col) >= col_start & 
                    !!sym(col) <= col_stop) %>%
    dplyr::pull(Coordinate)
  
  return(setdiff(neighbors_list, barcode))
}


calculate_slide_spatial_distances_my <- function(data_seurat, filter_coord_list, filename, run_min=TRUE){
  ## Calculate pairwise spatial distance between every point in the spatial image
  
  # Extract spot positions
  physical_coordinates <- tibble::as_tibble(data_seurat@images$slice1@coordinates[c("row", "col")], rownames="rowname") %>%
    dplyr::rename(Coordinate = rowname)

  # Define all pairs of spots
  coord_combos <- combn(physical_coordinates$Coordinate, 2)
  
  # Create a tibble to contain all pairwise distances on a slide
  coord_df <- tibble::tibble("Coordinate1" = coord_combos[1,],
                             "Coordinate2" = coord_combos[2,],
                             "Distance" = 0) %>%
    dplyr::filter(Coordinate1 %in% filter_coord_list &
                  Coordinate2 %in% filter_coord_list)
  
  for(i in 1:nrow(coord_df)){
    coord1 <- physical_coordinates %>%
      dplyr::filter(Coordinate == coord_df[[i, "Coordinate1"]])
    coord2 <- physical_coordinates %>%
      dplyr::filter(Coordinate == coord_df[[i, "Coordinate2"]])
    
    coord_df[[i, "Distance"]] <- raster::pointDistance(c(coord1$row, coord1$col), 
                                                       c(coord2$row, coord2$col), lonlat=FALSE)
  }
  
  if(run_min == FALSE){
    # Write distances to file
    filename <- paste0(global_path, "Input_data_processed/", sample_name, "_dist_pairs.txt")
    readr::write_delim(coord_df, path=filename, delim = "\t", append=FALSE, col_names = TRUE)
  }
}


quadrant_group_spots_on_two_signatures_my <- function(signature_wide_df, sig_stat_df, sig_center, sig_area){
  ## Partition spots into hi-lo quadrants based on scores of two signatures
  
  # Stat values for the 1st signature
  sig_stat1 <- sig_stat_df %>% 
    dplyr::filter(Signature_name == sig_center)
  
  # Stat values for the 2nd signature
  sig_stat2 <- sig_stat_df %>% 
    dplyr::filter(Signature_name == sig_area)
  
  # Partition spots into hi-lo quadrants based on scores of two signatures
  sig_select_df <- signature_wide_df %>% 
    dplyr::select(c("Coordinate", sig_center, sig_area)) %>%
    dplyr::mutate(Sig1_Hi = ifelse(!!sym(sig_center) >= sig_stat1$score_mean + sig_stat1$score_sd, "Hi_", ""),
                  Sig1_Lo = ifelse(!!sym(sig_center) <= sig_stat1$score_mean - sig_stat1$score_sd, "Lo_", ""),
                  Sig2_Hi = ifelse(!!sym(sig_area) >= sig_stat2$score_mean + sig_stat2$score_sd, "Hi", ""),
                  Sig2_Lo = ifelse(!!sym(sig_area) <= sig_stat2$score_mean - sig_stat2$score_sd, "Lo", ""),
                  Category = paste0(Sig1_Hi, Sig1_Lo, Sig2_Hi, Sig2_Lo))
  
  # Number of spots in each quadrant
  group_sizes <- sig_select_df %>%
    dplyr::group_by(Category) %>%
    dplyr::summarise(Count = dplyr::n_distinct(Coordinate))
  
  head(group_sizes)
  
  return(sig_select_df %>%
           dplyr::select(c("Coordinate", "Category", sig_center, sig_area)))
}


reduce_visualize_distances_my <- function(coord_df, output_folder){
  ## Create histograms of various distance reductions

  # Create a histogram of all signature scores to get a sense of the distribution
  filename <- paste0(output_folder, "/hist_group_distance")
  create_hist_plot_my(coord_df, x_label="Distance", fill_label="Category_pair", 
                      facet_var=c("Category_pair", "fixed"), intercept=0, binwidth=1, 
                      filename=filename, labels=c("Distance", "Number of spots", ""))
  
  # Distance to the closest spot with a different category
  coord_reduce_df <- coord_df %>%
    dplyr::filter(Category_center != Category_area) %>%
    dplyr::group_by(Coordinate1) %>%
    dplyr::summarise(Distance_min = min(Distance),
                     closest_category = paste0(unique(ifelse(Distance == Distance_min, Category_pair, "")), collapse=""),)
  
  # Create a histogram of all signature scores to get a sense of the distribution
  filename <- paste0(output_folder, "/hist_group_distance_min_opposite")
  create_hist_plot_my(coord_reduce_df %>%
                        dplyr::mutate(Distance_min = ifelse(Distance_min <= 10, Distance_min, 10)), 
                      x_label="Distance_min", fill_label="closest_category", 
                      facet_var=c("closest_category", "fixed"), intercept=0, binwidth=0.5, 
                      filename=filename, labels=c("Distance", "Number of spots", ""))
  
  # Distance to the closest spot with the same category
  coord_reduce_df <- coord_df %>%
    dplyr::filter(Category_center == Category_area) %>%
    dplyr::group_by(Coordinate1) %>%
    dplyr::summarise(Distance_min = min(Distance),
                     closest_category = paste0(unique(ifelse(Distance == Distance_min, Category_pair, "")), collapse=""),)
  
  # Create a histogram of all signature scores to get a sense of the distribution
  filename <- paste0(output_folder, "/hist_group_distance_min_same")
  create_hist_plot_my(coord_reduce_df %>%
                        dplyr::mutate(Distance_min = ifelse(Distance_min <= 10, Distance_min, 10)), 
                      x_label="Distance_min", fill_label="closest_category", 
                      facet_var=c("closest_category", "fixed"), intercept=0, binwidth=0.5, 
                      filename=filename, labels=c("Distance", "Number of spots", ""))
}


spatial_subgroups_threshold_my <- function(){
  ## THIS FUNCTION IS PROBABLY DEPRECIATED
  
  # Calculate pairwise spatial distance between every point in the spatial image
  filename <- paste0(global_path, "Input_data_processed/", sample_name, "_dist_pairs.txt")
  if(!file.exists(filename)){
    coord_init_df <- calculate_slide_spatial_distances_my(data_seurat, sig_select_df$Coordinate, filename, run_min)
  } else{
    coord_init_df <- read_file2df_my(filename)
  }
  
  # Trim coordinates only to those of interest
  coord_df <- coord_init_df %>%
    dplyr::inner_join(sig_select_df %>%
                        dplyr::rename(Category_center = Category), 
                      by=c("Coordinate1" = "Coordinate")) %>%
    dplyr::inner_join(sig_select_df %>%
                        dplyr::rename(Category_area = Category), 
                      by=c("Coordinate2" = "Coordinate")) %>%
    dplyr::mutate(Category_pair = paste0(Category_center, "-", Category_area))
  
  # Create histograms of various distance reductions
  if(run_min == FALSE){
    reduce_visualize_distances_my(coord_df, output_folders$group_dea)
  }
  
  # Define neighborhood radius
  # raster::pointDistance(c(0,0), c(1,1), lonlat=FALSE)
  # 1 cell away = 1.42; 2 cells away = 2.83; 3 cells away = 4.25
  neighborhood_radius <- 5
  lib_size <- 15000
  
  # Define groups of interest
  coord_reduce_df <- coord_df %>%
    dplyr::filter(Distance <= neighborhood_radius) %>%
    dplyr::mutate(Coordinate_ID = paste0(Coordinate1, "_", Category_center)) %>%
    dplyr::group_by(Coordinate_ID, Category_area) %>%
    dplyr::summarise(Count = dplyr::n_distinct(Coordinate2)) %>%
    dplyr::ungroup() %>%
    df_long2wide_my(rows="Coordinate_ID", cols="Category_area", value="Count") %>%
    tidyr::replace_na(list(Hi_Lo = 0, Lo_Hi = 0)) %>%
    dplyr::mutate(Ratio = log2(round((Hi_Lo+1) / (Lo_Hi+1), 2)),
                  Category_neighbor = ifelse(abs(Ratio) > neighborhood_radius, "Outlier", "Else")) %>%
    dplyr::mutate(Coordinate = gsub("_.+$", "", Coordinate_ID),
                  Category_spot = gsub("^.+-1_", "", Coordinate_ID)) %>%
    dplyr::inner_join(signature_wide_df %>%
                        dplyr::select(c(Coordinate, nCount_RNA)), by="Coordinate") %>%
    dplyr::filter(nCount_RNA <= lib_size)
  
  # Create a histogram of the ratio between cells with one signature over the other
  if(run_min == FALSE){
    filename <- paste0(output_folders$group_dea, "/hist_group_distance_ratio")
    create_hist_plot_my(coord_reduce_df, 
                        x_label="Ratio", fill_label="Category_neighbor", 
                        facet_var=NULL, intercept=0, binwidth=0.5, 
                        filename=filename, labels=c("Ratio, log2", "Number of spots", ""), 
                        add_density=FALSE, log_scale=FALSE)
  }
  
  # Not a single high stroma spot within the radius; high cancer content
  cancer_no_stroma <- coord_reduce_df %>% 
    dplyr::filter(Category_spot == "Hi_Lo" & Lo_Hi == 0)
  
  # At least one high stroma spot with the radius; low high cancer content
  cancer_w_stroma <- coord_reduce_df %>% 
    dplyr::filter(Category_spot == "Hi_Lo" & Lo_Hi > 0)
  
}
