# Author: Anna Lyubetskaya. Date: 21-09-03
# Calculate distance between every spot and a pathology identified area / event / edge and find expression correlates


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

library(Seurat)
if (!("geodiv" %in% installed.packages())){
  install.packages("geodiv")
}
library(geodiv)
library(raster)

source("code/R/Utils/utils_pathology_halo.R")
source("code/R/Utils/utils_10X_vis.R")


## PARAMETERS ----


# Run name
run_name <- "PDAC_TLS"

# Sample / Cohort name
cohort_name <- "PDAC"

# Baseline class
reference_class <- "TLS"
reference_class_list <- c("TLS-Mature", "TLS-Immature")

# Remove holes from the pathology layer
do_fill <- "filled"

# If TRUE, calculate internal feature distance, if FALSE, calculate external feature distance
invert_mask <- FALSE

# User-defined name of the distance variable to add to the Seurat meta data
pathology_var_name <- paste0("Pathology.Distance.", reference_class, ".", do_fill, ".invert", invert_mask)

# Don't write the final RDS object
no_rds_output <- TRUE

# List of specific samples to select
sample_list <- c("PDAC_E5058_ROI1", "PDAC_E5058_ROI2", "PDAC_E5058_ROI3", "PDAC_E5058_ROI4",
                 "PDAC_Pt11_ROI3", "PDAC_Pt2_ROI2", "PDAC_Pt3_ROI1_s1", "PDAC_Pt3_ROI2_s1",
                 "PDAC_Pt4_ROI1", "PDAC_Pt6_ROI2", "PDAC_Pt6_ROI4", "PDAC_Pt8_ROI3")

# Merge and rename any of the pathology annotation classes
attr_dict <- list(reference_class = reference_class_list)
names(attr_dict) <- reference_class


## PATHS ----


# Input file with sample meta data
meta_file <- "XXXX"

# Input folder
input_path <- "XXXX"

# Output folder for RDS objects
output_path <- "XXXX"

# Output folder for images
output_init_figs <- "XXXX"
output_figs <- paste0(output_init_figs, run_name, "/")

# Create the output folder for all figures
dir.create(output_path, showWarnings = FALSE)
dir.create(output_init_figs, showWarnings = FALSE)
dir.create(output_figs, showWarnings = FALSE)


## INGEST DATA ----


# Read sample file; select only cohort-relevant samples
meta_df <- readr::read_delim(file=meta_file, delim="\t") %>%
  tidyr::drop_na() %>%
  dplyr::filter(grepl(cohort_name, Sample_Name) &
                  Pathology_Annotation_Name == "Compartments")

# Filter meta data for samples of interest if requested
if(!is.null(sample_list)){
  meta_df <- meta_df %>%
    dplyr::filter(Sample_Name %in% sample_list)
}


## WRANGLE DATA ----


for(row_i in 1){#:nrow(meta_df)){
  
  # Sample name
  sample_name <- meta_df[[row_i, "Sample_Name"]]
  # Pathology file to integrate
  path_file <- meta_df[[row_i, "Pathology_Full_Path"]]
  
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
    image_structure <- data_seurat@images$slice1
    
    # Seurat spot center coordinates in full resolution: Y = imagerow, X = imagecol
    X <- image_structure@coordinates[["imagecol"]]
    Y <- image_structure@coordinates[["imagerow"]]
    
    # Define spatstat window
    W_ss <- c(min(X), max(X), min(Y), max(Y))
    
    # Spot diameter at full resolution
    spot_radius <- image_structure@scale.factors$spot_diameter_fullres / 2
    
    
    ## LOAD PATHOLOGY DATA ----
    
    
    # Load the annotations as a list of owin-s
    polys_init <- halo_annot_to_poly_my(path_file, attr_dict=attr_dict)
    pathology_classes_init <- sort(names(polys_init))
    
    print(pathology_classes_init)
    
    ## CONVERT TO RASTER
    
    Mode <- function(x, na.rm=TRUE){
      if (na.rm){
        x <- x[!is.na(x)]
      }
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    
    seurat_dataframe <- data.frame(x=data_seurat@images$slice1@coordinates[,'col'], y=data_seurat@images$slice1@coordinates[,'row'], data_seurat@meta.data[-1])
    seurat_dataframe$Pathology.Group <- as.factor(seurat_dataframe$Pathology.Group)
    pathology_group_key <- seurat_dataframe$Pathology.Group
    
    meta_raster <- rasterFromXYZ(sapply(seurat_dataframe, as.numeric))
    
    # mean interpolation
    for(i in which(names(meta_raster) != "Pathology.Group")){
      meta_raster <- raster::setValues(meta_raster, raster::getValues(raster::focal(meta_raster[[i]], w=matrix(1,nrow=3,ncol=3), fun=mean, NAonly=TRUE, na.rm=TRUE)), layer=i)
    }
    # mode interpolation
    for(i in which(names(meta_raster) == "Pathology.Group")){
      meta_raster <- raster::setValues(meta_raster, raster::getValues(raster::focal(meta_raster[[i]], w=matrix(1,nrow=3,ncol=3), fun=Mode, NAonly=TRUE, na.rm=TRUE)), layer=i)
    }
    
    
    
    
    
    
    
    seurat_dataframe <- data.frame(x=data_seurat@images$slice1@coordinates[,'col'], y=data_seurat@images$slice1@coordinates[,'row'], data_seurat@meta.data[-1])
    
    gene_raster <- rasterFromXYZ()
    
    raster::focal(my_raster[["nFeature_Spatial"]], w=matrix(1, nrow = 3, ncol = 3), fun=mean, NAonly=TRUE, na.rm=TRUE)
    
    
    
    
    
    
    
    
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
    
    
    # Distance object
    dist_fun <- spatstat.geom::distfun.owin(polys, invert=invert_mask)
    
    # Report a summary for the computed distance function
    summary(dist_fun)
    
    # Calculate distance between every spot center and the pathology feature
    # Parallelize to save time
    cl <- parallel::makeCluster(parallel::detectCores()*2)
    
    doParallel::registerDoParallel(cl)
    
    distance_list <- foreach::foreach(i = 1:length(X), .combine=c) %dopar% {
      
      # Identify X and Y ranges for this spot
      x_list <- (X[i] - spot_radius) : (X[i] + spot_radius)
      y_list <- (Y[i] - spot_radius) : (Y[i] + spot_radius)
      
      # Create a grid of all coordinate values as a square around the spot center
      square_coord <- expand.grid(x_list, y_list)
      
      # Filter the square down to circle
      square_coord <- square_coord[which(sqrt((square_coord$Var1 - X[i])^2 + (square_coord$Var2 - Y[i])^2) <= spot_radius),]
      
      # Calculate mean distance across the spot
      spot_centers_ppp <- spatstat.geom::as.ppp(list(x=square_coord$Var1, y=square_coord$Var2), W_ss)
      
      mean(dist_fun(spot_centers_ppp))
      
    };
    
    parallel::stopCluster(cl)
    
    
    # Add distance to the pathology feature as a meta data factor of the Seurat object
    data_seurat@meta.data[[pathology_var_name]] <- distance_list
    
    
    ## VISUALIZE DATA ----
    
    
    # Plot pathology feature to which we are calculating the distance
    filename <- paste0(output_figs, "pathology_feature_", sample_name, "_", pathology_var_name, ".png")
    p2 <- plot_pathology_ingestion_check_my(polys_init[reference_class], X, Y, spot_radius, filename, main_class=reference_class, title="",
                                            class_list=c(reference_class))
    
    # Spatial visualization of the distance
    p1 <- spatial_feature_plot_my(data_seurat, pathology_var_name, min.cutoff="q0", max.cutoff="q100", name=pathology_var_name, crop=TRUE)
    
    # Save distance visualization to file
    filename <- paste0(output_figs, "pathology_dist_", sample_name, "_", pathology_var_name)
    write_plot2file_my(p1, filename, num_row=1, num_col=1)
    
    
    ## OVERWRITE SEURAT FILE ----
    
    
    if(no_rds_output == FALSE){
      # Write the updated Seurat object
      filename <- paste0(output_path, sample_name, "_annotated_pathology.rds")
      saveRDS(data_seurat, file=filename)
    }
  } else{
    warning("More than one file found!\n", file_list)
  }
}
