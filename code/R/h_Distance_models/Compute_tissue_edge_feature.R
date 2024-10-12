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
run_name <- "PDAC_TxNaive_TissueDistances"

# Sample / Cohort name
cohort_name <- "PDAC"

# All classes to process
expected_classes <- c("TME_v14 - Background",
                      "TME_v14 - White Space",
                      "Tissue")

# Rename classes to new labels if necessary
expected_classes_renamed <- c("Background",
                              "Background",
                              "Tissue")

if(is.null(expected_classes_renamed)){
  expected_classes_renamed <- expected_classes
}
name_dict <- expected_classes_renamed
names(name_dict) <- expected_classes

# set classes for computation (based on simpler names in dict above)
desired_class_list <- c("Tissue")

# Remove holes from the pathology layer
do_fill <- "filled" # 'filled' to do, 'unfilled' (or anything else) to not. becomes part of filename

# If TRUE, calculate internal feature distance, if FALSE, calculate external feature distance
invert_mask <- TRUE

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


for(row_i in 1:nrow(meta_df)){
  
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
    
    
    # User-defined name of the distance variable to add to the Seurat meta data
    pathology_var_name <- paste0("Pathology.Distance.Tissue.", do_fill, ".invert", invert_mask)
    
    # Trim "background" from the tissue class
    tissue_owin <- spatstat.geom::setminus.owin(polys_filt[["Tissue"]], polys_filt[["Background"]])
    
    # plotting for debugging/dev
    if(FALSE){
      par(mfrow=c(1,3))
      plot(polys_filt[["Tissue"]], col = 'blue')
      plot(polys_filt[["Background"]], col = 'green')
      plot(tissue_owin, col = 'red')
    }
    
    # Fill holes in the tissue
    if(do_fill == "filled"){
      polys_temp <- tissue_owin
      
      # Modify the bdry list (list of polygon shapes as xy coordinates) by logically filtering out the holes
      polys_temp$bdry <- polys_temp$bdry[!as.logical(lapply(polys_temp$bdry, spatstat.utils::is.hole.xypolygon))]
      
      # Split the owin object into separate shapes
      shape_list <- purrr::map(polys_temp$bdry, ~ spatstat.geom::owin(poly = .x))
      # Find the shape union
      polys_temp <- spatstat.geom::union.owin(spatstat.geom::as.solist(shape_list))
      
      tissue_owin <- polys_temp
    }
    
    # Select the pathology object to which calculate the distance
    polys <- spatstat.geom::as.owin(tissue_owin)
    # Ensure the domain for the distance function encompasses all spots
    polys$xrange <- polys_filt[["Tissue"]]$xrange + c(-500, 500)
    polys$yrange <- polys_filt[["Tissue"]]$yrange + c(-500, 500)
    
    ## CALCULATE DISTANCE ----
    
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
    
    # Add distance to the pathology feature as a meta data factor of the Seurat object
    data_seurat@meta.data[[pathology_var_name]] <- distance_list
    
    
    ## VISUALIZE DATA ----
    
    
    # Plot pathology feature to which we are calculating the distance
    filename <- file.path(output_figs, paste0("pathology_feature_", sample_name, "_", pathology_var_name, ".png"))
    p2 <- plot_pathology_ingestion_check_my(polys_init['Tissue'], X, Y, spot_radius, filename, main_class="Tissue", title="",
                                            class_list="Tissue")
    
    # Spatial visualization of the distance
    p1 <- spatial_feature_plot_my(data_seurat, pathology_var_name, min.cutoff="q0", max.cutoff="q100", crop=TRUE)
    
    
    # Save distance visualization to file
    filename <- file.path(output_figs, paste0("pathology_dist_", sample_name, "_", pathology_var_name))
    write_plot2file_my(p1, filename, num_row=1, num_col=1)
    
    
    
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
