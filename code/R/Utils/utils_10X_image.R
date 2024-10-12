# Author: Anna Lyubetskaya. Date: 20-12-23

## 10X and Seurat deal with 3 different image resolution for an ST capture area:
## - "original 10X resolution" - barcode image coordinates are provided in this coordinate system
## -- Barcodes are associated with circles of a certain radius mapped to this image
## - "high 10X resolution" - an image compressed down to ~2000x2000 pixels - generally speaking, this resolution is not really used
## - "low 10X resolution" - an image compressed down to ~600x600 pixels - this resolution is used by Seurat for all visualizations


# source("code/utils/utils_stats.R")

if(packageVersion("spatstat") < 2){
  update.packages("spatstat")
}


update_seurat_image_properties_from_json_my <- function(data_seurat, folder_10X){
  ## Load image JSON, calculate spot diameter in low 10X resolution, and add that to the Seurat object
  
  image_json <- rjson::fromJSON(file = paste0(folder_10X, "/spatial/scalefactors_json.json"))
  
  spot_diameter_fullres <- image_json$spot_diameter_fullres
  spot_radius_lowres <- floor(image_json$spot_diameter_fullres * image_json$tissue_lowres_scalef / 2)
  
  data_seurat@images$slice1@scale.factors[["spot_diameter_fullres"]] <- spot_diameter_fullres
  data_seurat@images$slice1@scale.factors[["spot_radius_lowres"]] <- spot_radius_lowres
  
  return(data_seurat)
}


seurat_spot_neigbors_my <- function(data_seurat=NULL, image=NULL, spot_start=0, spot_distance=2, n=1){
  ## For each Seurat spot find all neighbors within a given distance
  
  # Extract Seurat coordinates
  if(!is.null(data_seurat)){
    xy <- data_seurat@images[[n]]@coordinates[, c("row", "col")]
  } else{
    xy <- image@coordinates[, c("row", "col")]
  }
  
  barcode_list <- rownames(xy)
  xy <- mapply(xy, FUN=as.numeric)
  
  # Find the matrix with the indices of points belonging to the set of the k nearest neighbours of each other
  spot_nb <- spdep::dnearneigh(xy, spot_start, spot_distance, row.names=barcode_list)
  
  return(spot_nb) 
}


seurat_image_to_clusters_my <- function(data_seurat, spot_distance=2, n=1){
  ## Identify contiguous clusters from an ST Seurat image using a fixed distance between spots
  
  # For each Seurat spot find all neighbors within a given distance
  spot_nb <- seurat_spot_neigbors_my(data_seurat, spot_distance=spot_distance, n=n)
  
  # Transform neighborhoods into a weighted list
  spot_nb_list <- as(spdep::nb2mat(spot_nb, style="B", zero.policy=TRUE), "CsparseMatrix")
  
  # Calculate graph adjacency
  ## set add.rownames = NULL to pass the cell ID to the output dgcMatrix
  spot_graph <- igraph::graph.adjacency(spot_nb_list, mode="undirected", add.rownames = NULL)
  
  # Identify graph subclusters
  spot_clusters <- igraph::clusters(spot_graph)
  
  return(spot_clusters)
}


tissue_contiguity_filter_my <- function(data_seurat, spot_num=200){
  ## Subset Seurat data to only contiguious clusters above certain size
  ## Working on the first image only
  ## Assuming that the barcode order is preserved
  
  ## https://www.rdocumentation.org/packages/spdep/versions/1.1-5
  ## https://r-spatial.github.io/spdep/reference/dnearneigh.html
  ## https://www.rdocumentation.org/packages/spdep/versions/1.1-3/topics/Graph%20Components
  
  # Identify contiguious clusters from an ST Seurat image
  spot_clusters <- seurat_image_to_clusters_my(data_seurat)
  
  if(length(spot_clusters$csize)==1){
    write.csv(spot_clusters, paste0(output_path, "sample_name", "_spot_clusters.csv"))
  }
  
  # Select clusters of spots above a certain size
  cluster_index <- which(spot_clusters$csize > spot_num)
  
  # Select barcodes within big clusters
  barcode_list <- rownames(data_seurat@images[[1]]@coordinates)[which(spot_clusters$membership %in% cluster_index)]
  
  # Subset the Seurat object
  data_seurat_select <- subset(data_seurat, cells=barcode_list)
  
  return(data_seurat_select)
}


rgb_matrix_to_color_value_df_my <- function(rgb_matrix){
  ## Turn an RGB matrix into a long matrix of color values
  
  # Find image dimensions in low resolution
  image_dim <- dim(rgb_matrix[,,1])
  
  # Find color values from Seurat low resolution RGB rendering of the image
  # Reminder 1: data_seurat@images$slice1@image is a matrix dim = [600, Y, 3] corresponding to 3 channel rendering, RGB
  # Reminder 2: when converting a 2D matrix to 1D vector, the vector contains the elements in the following order: column1: row 1 to n; column2: row1 to n, etc
  color_value_vec <- rgb2hsv(as.vector(rgb_matrix[,,1]), as.vector(rgb_matrix[,,2]), as.vector(rgb_matrix[,,3]))[3, ]
  
  # Create a tibble of color values from Seurat low resolution image with explicit row and column indices
  color_value_df <- tibble::tibble(Y = rep(1:image_dim[1], image_dim[2]),
                                   X = as.vector(sapply(1:image_dim[2], function(x) rep(x, image_dim[1]))),
                                   Value = color_value_vec)
  
  return(color_value_df)
}


nearest_neighbor_df_my <- function(df, ref_df, value="Barcode"){
  ## Find the nearest neighbors of a set of query points from a tibble in a reference tibble
  
  # Calculate the distance between every low resolution image coordinate and its nearest spot center
  nearest_dist <- Biobase::matchpt(as.matrix(df[c("X", "Y")]), as.matrix(ref_df[c("X", "Y")]))
  
  # Update image data with the nearest barcode information and distance
  df[[value]] <- ref_df[[value]][nearest_dist$index]
  df[["Distance"]] <- nearest_dist$distance
  
  return(df)
}


image_mean_color_value_my <- function(data_seurat, num_sd=1){
  ## Analyze color values in Seuart image matrix and identify spots that seem too light to be under tissue
  
  # Seurat image object
  image_structure <- data_seurat@images$slice1
  
  # Low resolution coefficient - allows to translate original 10X coordinates to low-resolution 10X coordinates
  coef <- image_structure@scale.factors$lowres
  # Spot radius in low resolution
  spot_radius <- image_structure@scale.factors$spot_radius_lowres
  
  # Seurat spot center coordinates in low resolution
  spot_df <- tibble::as_tibble(image_structure@coordinates, rownames="Barcode") %>%
    dplyr::mutate(Y = coef * imagerow,
                  X = coef * imagecol) %>%
    dplyr::select(-tissue, -row, -col, -imagerow, -imagecol)
  
  # Turn an RGB matrix into a long matrix of color values
  color_value_df <- rgb_matrix_to_color_value_df_my(image_structure@image)
  
  # Find the nearest neighbors of a set of query points from a tibble in a reference tibble
  color_value_df <- nearest_neighbor_df_my(color_value_df, spot_df, value="Barcode")
  
  # Map color data to spot center data
  integrated_df <- color_value_df %>%
    dplyr::filter(Distance <= spot_radius) %>%
    dplyr::inner_join(spot_df, by="Barcode")
  
  # Find pixels not under the tissue
  other_df <- color_value_df %>%
    dplyr::filter(Distance >= spot_radius*2,
                  Value >= mean(integrated_df$Value))
  
  # Using spots not under the tissue, find a color value threshold
  color_value_cutoff <- mean(other_df$Value) - sd(other_df$Value) * num_sd
  
  # Analyze a vector density distribution: mode and local minima; identify an appropriate color value cutoff
  # density_stats <- density_analysis_my(integrated_df$Value)
  
  # Check that the color cutoff is higher than the last mode identified in the density distribution
  # x <- length(density_stats$mode_x)
  # color_value_cutoff <- max(color_value_cutoff, density_stats$mode_x[x])
  
  # Select the part of the image with signal and calculate mean color value for each barcode
  spot_color_df <- integrated_df %>%
    dplyr::group_by(Barcode) %>%
    dplyr::summarise(Value_mean = mean(Value)) %>%
    dplyr::mutate(Exclude = Value_mean >= color_value_cutoff)
  
  return(spot_color_df)
}
