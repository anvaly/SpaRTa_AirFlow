# Author: Andrew Fisher. Date: 23-08-01
# Calculate distance between every spot and a pathology identified area / event / edge and find expression correlates


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

library(Seurat)
library(raster)

source("code/R/Utils/utils_pathology_halo.R")
source("code/R/Utils/utils_10X_vis.R")


## PARAMETERS ----


# Run name
run_name <- "PDAC_TxNaive_PathDistances"

# Sample / Cohort name
cohort_name <- "PDAC"

# All classes to process
expected_classes <- c("Adipose_v4",
                      "Adjacent Intestine",
                      "Lymph_Node",
                      "Stroma_v4 - Muscle",
                      "Adjacent Muscle",
                      "Stroma_v4 - Nerve",
                      "Adjacent NonTumor Tissue",
                      "Adjacent Serosa",
                      "TLS-Aggregate",
                      "TLS-Immature",
                      "TLS-Mature",
                      "Stroma_v4 - Vessel",
                      "TME_v14 - Benign Epithelium",
                      "TME_v14 - Blood",
                      "TME_v14 - Exocrine and Endocrine",
                      "TME_v14 - Luminal Content",
                      "TME_v14 - Non-Epithelium",
                      "TME_v14 - Tumor",
                      "TME_v14 - Background",
                      "TME_v14 - White Space",
                      "Tissue",
                      "Blur_Tissue_v8",
                      "NonBlur_Tissue_v8")

# Rename classes to new labels if necessary
expected_classes_renamed <- c("Adipose",
                              "IntestineAdj",
                              "LymphNode",
                              "Muscle",
                              "MuscleAdj",
                              "Nerve",
                              "NormalAdj",
                              "NormalAdj",
                              "TLSAggregate",
                              "TLSImmature",
                              "TLSMature",
                              "Vessel",
                              "BenignEpi",
                              "Blood",
                              "ExoEndo",
                              "LuminalNec",
                              "NonEpi",
                              "Tumor",
                              "Background",
                              "Background",
                              "Tissue",
                              "TissueBlur",
                              "TissueNoBlur")

if(is.null(expected_classes_renamed)){
  expected_classes_renamed <- expected_classes
}
name_dict <- expected_classes_renamed
names(name_dict) <- expected_classes

# set classes for computation (based on simpler names in dict above)
desired_class_list <- c("Adipose",
                        "IntestineAdj",
                        "LymphNode",
                        "Muscle",
                        "MuscleAdj",
                        "Nerve",
                        "NormalAdj",
                        "TLSAggregate",
                        "TLSImmature",
                        "TLSMature",
                        "Vessel",
                        "BenignEpi",
                        "Blood",
                        "ExoEndo",
                        "LuminalNec",
                        "NonEpi",
                        "Tumor")

# Remove holes from the pathology layer
do_fill <- "unfilled" # 'filled' to do, 'unfilled' (or anything else) to not. becomes part of filename

# choose distance map computation method
# note at this time inverted distances (depths) not coded in raster mode
compute_mode <- "spatstat" # choose "spatstat" or "raster"

# If TRUE, calculate internal feature distance, if FALSE, calculate external feature distance
invert_mask <- FALSE

# Write the final RDS object
do_rds_output <- FALSE
# Write to limited csv object
do_csv_output <- TRUE

# List of specific samples to select (or NULL)
sample_list <- NULL  # c("PDAC_E1265_ROI1", "PDAC_E1265_ROI2", "PDAC_E1265_ROI3", "PDAC_E1265_ROI4")


## PATHS ----


# Input file with sample meta data
meta_file <- "XXXX"

# Input folder
input_path <- "XXXX"

# Output folder for RDS objects
output_base_path <- "XXXX"
output_path <- file.path(output_base_path, run_name)

# Output folder for images
output_figs <- file.path(output_path, 'figs')

# Create the output folder for all figures
dir.create(output_path, showWarnings = FALSE)
dir.create(output_figs, showWarnings = FALSE)


## INGEST DATA ----


# Read sample file; select only cohort-relevant samples
meta_df <- readr::read_delim(file=meta_file, delim="\t") %>%
  dplyr::filter(grepl(cohort_name, Sample_Name) &
                  Best_Looking == TRUE)

# Filter meta data for samples of interest if requested
if(!is.null(sample_list)){
  meta_df <- meta_df %>%
    dplyr::filter(Sample_Name %in% sample_list)
}


## WRANGLE DATA ----


for(row_i in c(9,40,75,76)){       #77:nrow(meta_df)){
  
  # Sample name
  sample_name <- meta_df[[row_i, "Sample_Name"]]
  # Pathology file to integrate
  path_file <- meta_df[[row_i, "Pathology_Compartment_File"]]
  
  # Find the RDS object
  file_list <- dir(input_path, pattern=paste0(sample_name), full.names=TRUE, recursive=TRUE)
  
  if(length(file_list) == 1){
    
    print(row_i)
    print(file_list[1])
    
    
    ## INGEST SEURAT DATA ----
    
    
    # Open a connection to the RDS object
    con <- gzfile(file_list[1])
    
    # Ingest the Seurat object
    data_seurat <- readRDS(con)
    
    # Close the connection to be able to overwrite
    close(con)
    
    
    ## EXTRACT SEURAT SPOT DATA ----
    
    
    # Seurat image object
    image_structure <- data_seurat@images[[1]]
    
    # Seurat spot center coordinates in full resolution: Y = imagerow, X = imagecol
    X <- image_structure@coordinates[["imagecol"]]
    Y <- image_structure@coordinates[["imagerow"]]
    
    # Define spatstat window
    W_ss <- c(min(X), max(X), min(Y), max(Y))
    
    # Spot diameter at full resolution (in pixels)
    spot_radius <- image_structure@scale.factors$spot_diameter_fullres / 2
    
    
    ## LOAD PATHOLOGY DATA ----
    
    
    # Load the annotations as a list of owin's
    # only keep expected_classes, extraneous layers are unpredictable
    polys_init <- halo_annot_to_poly_my(path_file, attr_dict = name_dict, region_list = expected_classes)
    pathology_classes_init <- sort(names(polys_init))
    
    # Drop empty annotation layers
    polys_filt <- polys_init[sapply(polys_init, function(x) length(x$bdry)>0)]
    pathology_classes_filt <- sort(names(polys_filt))
    
    # Now intersect the desired computation classes with those available
    compute_class_list <- intersect(desired_class_list, pathology_classes_filt)
    print(compute_class_list)
    
    # Looping through the available pathology classes for computation
    for (reference_class in compute_class_list){
      # User-defined name of the distance variable to add to the Seurat meta data
      pathology_var_name <- paste0("Pathology.Distance.", reference_class, ".", do_fill, ".invert", invert_mask)
      
      
      # Fill holes in the tissue
      if(do_fill == "filled"){
        polys_temp <- polys_init[[reference_class]]
        
        # Modify the bdry list (list of polygon shapes as xy coordinates) by logically filtering out the holes
        polys_temp$bdry <- polys_temp$bdry[!as.logical(lapply(polys_temp$bdry, spatstat.utils::is.hole.xypolygon))]
        
        # Split the owin object into separate shapes
        shape_list <- purrr::map(polys_temp$bdry, ~ spatstat.geom::owin(poly = .x))
        # Find the shape union
        polys_temp <- spatstat.geom::union.owin(spatstat.geom::as.solist(shape_list))
        
        polys_init[[reference_class]] <- polys_temp
      }
      
      # Select the pathology object to which calculate the distance
      polys <- polys_init[[reference_class]]
      # Ensure the domain for the distance function encompasses all spots
      polys$xrange <- polys_init[["Tissue"]]$xrange + c(-500, 500)
      polys$yrange <- polys_init[["Tissue"]]$yrange + c(-500, 500)
      
      ## CALCULATE DISTANCE ----
      
      if(compute_mode=="spatstat"){
        
        # Distance object
        dist_fun <- spatstat.geom::distfun.owin(polys, invert=invert_mask)
        
        # Report a summary for the computed distance function
        summary(dist_fun)
        
        # Calculate distance between every spot center and the pathology feature
        # Parallelize to save time (reserve at least one core so things don't go sideways)
        cl <- parallel::makeCluster(max(1, parallel::detectCores()-1))
        
        doParallel::registerDoParallel(cl)
        
        distance_list <- foreach::foreach(i = 1:length(X), .combine=c) %dopar% {
          
          # Identify X and Y ranges for this spot
          x_list <- (X[i] - spot_radius) : (X[i] + spot_radius)
          y_list <- (Y[i] - spot_radius) : (Y[i] + spot_radius)
          
          # Create a grid of all coordinate values as a square around the spot center
          square_coord <- expand.grid(x_list, y_list)
          
          # Filter the square down to circle
          square_coord <- square_coord[which(sqrt((square_coord$Var1 - X[i])^2 + (square_coord$Var2 - Y[i])^2) <= spot_radius),]
          
          # Calculate mean distance across the spot (returns that value)
          spot_centers_ppp <- spatstat.geom::as.ppp(list(x=square_coord$Var1, y=square_coord$Var2), W_ss)
          
          mean(dist_fun(spot_centers_ppp))
          
        } # end foreach instructions
        
        parallel::stopCluster(cl)
        
      } else if(compute_mode=="raster"){
        # useful: https://gis.stackexchange.com/questions/210506/how-to-calculate-distance-on-large-raster-in-r
        # convert owin polys to sp polygons
        # reverse the xy list because sp sees 'holes' order differently than owin
        poly_list <- lapply(polys$bdry, function(owin_poly) {Polygon(cbind(rev(owin_poly$x), rev(owin_poly$y)))})
        sp_polys <- SpatialPolygons(list(Polygons(poly_list, ID = 1)))
        
        if(FALSE){
          ggplot()+
            geom_polygon(data = fortify(sp_polys), aes(x = long, y = lat, group = group))
        }
        
        # Create a raster layer with the same extent and resolution as the polygons
        raster_layer <- raster(extent(c(polys$xrange, polys$yrange)), res = 1)
        # default raster sets positions as half steps (1.5, 2.5, etc.), shift back to whole numbers
        raster_layer[] <- floor(raster_layer[])
        
        # Compute the distance function from the polygons to each pixel in the raster layer
        # visualizing resulting raster suggests something is wrong
        distance_function <- distanceFromPoints(raster_layer, coordinates(sp_polys))
        
        # try rasterized polygon approaches using focal and distance
        rast_poly <- raster::rasterize(sp_polys, raster_layer, field = 1, exact = TRUE)
        distances_focal <- focal(rast_poly, w = matrix(1, nrow = 3, ncol = 3), fun = function(x) min(x, na.rm=TRUE), pad = TRUE, padValue = NA)
        distance_function <- distance(rast_poly, units = "cell")
        
        
        
        # try mask suggestion from link (failed/too long on XL instance)
        raster_layer <- raster::setValues(raster_layer, 0) # per link
        mask_layer <- mask(raster_layer, sp_polys)
        distance_layer <- distance(mask_layer)
        
        if(FALSE){
          image(distance_function)
        }
        
        # Calculate distance between every spot center and the pathology feature
        # Parallelize to save time (reserve at least one core so things don't go sideways)
        cl <- parallel::makeCluster(max(1, parallel::detectCores()-1))
        
        doParallel::registerDoParallel(cl)
        
        distance_list <- foreach::foreach(i = 1:length(X), .combine=c) %dopar% {
          
          # Identify X and Y ranges for this spot
          x_list <- (X[i] - spot_radius) : (X[i] + spot_radius)
          y_list <- (Y[i] - spot_radius) : (Y[i] + spot_radius)
          
          # Create a grid of all coordinate values as a square around the spot center
          square_coord <- expand.grid(x_list, y_list)
          
          # Filter the square down to circle
          circle_coord <- square_coord[which(sqrt((square_coord$Var1 - X[i])^2 + (square_coord$Var2 - Y[i])^2) <= spot_radius),]
          
          # round the coords to nearest pixel
          circle_px_coord <- round(circle_coord)
          
          # use the raster distance layer to get values and average them
          extracted_distances <- raster::extract(distance_function, circle_px_coord)
          
          # return the average of all px distance values
          mean(extracted_distances)
          
        } # end foreach instructions
        
        parallel::stopCluster(cl)
        
      } else {
        # Don't compute anything, mode not selected
        print("Compute mode improperly defined. Neither SpatStat or Raster.")
        distance_list <- NULL
      }
      
      # Add distance to the pathology feature as a meta data factor of the Seurat object
      data_seurat@meta.data[[pathology_var_name]] <- distance_list
      
      
      ## VISUALIZE DATA ----
      
      
      # Plot pathology feature to which we are calculating the distance
      filename <- file.path(output_figs, paste0("pathology_feature_", sample_name, "_", pathology_var_name, ".png"))
      p2 <- plot_pathology_ingestion_check_my(polys_init[reference_class], X, Y, spot_radius, filename, main_class=reference_class, title="",
                                              class_list=c(reference_class))
      
      # Spatial visualization of the distance
      p1 <- spatial_feature_plot_my(data_seurat, pathology_var_name, min.cutoff="q0", max.cutoff="q100", crop=TRUE)
      
      
      # Save distance visualization to file
      filename <- file.path(output_figs, paste0("pathology_dist_", sample_name, "_", pathology_var_name))
      write_plot2file_my(p1, filename, num_row=1, num_col=1)
      
      
      
    } # end looping per pathology class
    
    ## OVERWRITE SEURAT FILE* ----
    # * or write to csv
    
    if(do_rds_output == TRUE){
      # Write the updated Seurat object
      filename <- file.path(output_path, paste0(sample_name, "_dist_pathology.rds"))
      saveRDS(data_seurat, file=filename)
    }
    if(do_csv_output == TRUE){
      # Write the distance data to a csv (only distance data)
      filename <- file.path(output_path, paste0(sample_name, "_dist_pathology.csv"))
      data_seurat@meta.data %>%
        dplyr::select(contains("Pathology.Distance")) %>%
        write.csv(file = filename, row.names = TRUE)
    }
  } else{
    warning("More than one file found!\n", file_list)
  }
}
