# Author: Andrew Fisher. Date: 23-08-17
# Integrates coregistered HALO exports to Visium spot observations
# Assumes object files containing coreg coordinates
# Code built specifically for trichrome assay integration
# Flow is to iteratively load through samples, match trichrome data to spots, then add appropriate metrics to observations


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

#library(foreach)
library(Seurat)
library(raster)

source("code/R/Utils/utils_10X_vis.R")


## PARAMETERS ----

# Run name
run_name <- "PDAC_TxNaive_Trichrome"

# Sample / Cohort name
cohort_name <- "PDAC"

# Write the final RDS object
do_rds_output <- FALSE
# Write to limited csv object
do_csv_output <- TRUE

# List of specific samples to select (or NULL)
sample_list <- c("PDAC_E1265_ROI1", "PDAC_E1265_ROI2", "PDAC_E1265_ROI3", "PDAC_E1265_ROI4")


## PATHS ----

# Input file with sample meta data
meta_file <- "XXXX"

# Input folder
input_path <- "XXXX"

# Input folder for Trichrome HALO exported directory
trichrome_path <- "XXXX"

# Output folder for RDS objects
output_base_path <- "XXXX"
output_path <- file.path(output_base_path, run_name)

# Output folder for images
output_figs <- file.path(output_path, 'figs')

# Create the output folder for all figures
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figs, recursive = TRUE, showWarnings = FALSE)


## INGEST DATA ----

# Read sample file; select only cohort-relevant samples
meta_df <- readr::read_delim(file=meta_file, delim="\t") %>%
  dplyr::filter(grepl(cohort_name, Sample_Name) &
                  Best_Looking == TRUE)

# Read HALO summary file (helps direct file loading per sample)
# contain both tissue summary and spot annotation summary results
halo_summary_df <-read.csv(file.path(trichrome_path, "Total_Summary_Results.csv"))

# Filter meta data for samples of interest if requested
if(!is.null(sample_list)){
  meta_df <- meta_df %>%
    dplyr::filter(Sample_Name %in% sample_list)
  # sample name was used to define spot annotations in trichrome coreg
  halo_summary_df <- halo_summary_df %>%
    dplyr::filter(Analysis.Region %in% meta_df$Sample_Name)
}


## PROCESS SAMPLES ----

for (row_i in 1:nrow(meta_df)){
  print(row_i)
  
  # Sample name
  sample_name <- meta_df[[row_i, "Sample_Name"]]
  # Trichrome file to integrate
  trichrome_job_id <- dplyr::filter(halo_summary_df, Analysis.Region==sample_name)[["Job.Id"]]
  trichrome_file <- dir(file.path(trichrome_path, "ObjectData"), pattern = paste0('job',trichrome_job_id), full.names = TRUE, recursive = TRUE)
  
  # Find the RDS object
  file_list <- dir(input_path, pattern=paste0(sample_name), full.names=TRUE, recursive=TRUE)
  
  # Ensure singular matching of found files
  if(length(file_list) == 1 & length(trichrome_file) == 1){
    
    print(file_list[1])
    print(trichrome_file[1])
    
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
    
    
    ## LOAD TRICHROME DATA ----
    
    
    # precalc the center of detected analysis objects
    halo_obj_df <- read.csv(trichrome_file) %>%
      dplyr::mutate(Visium.Registered.X = (Visium.Registered.XMin + Visium.Registered.XMax)/2) %>%
      dplyr::mutate(Visium.Registered.Y = (Visium.Registered.YMin + Visium.Registered.YMax)/2) %>%
      dplyr::mutate(X = mean(c(XMin, XMax))) %>%
      dplyr::mutate(Y = mean(c(YMin, YMax)))
      
    
    ## MAP TRICHROME OBJECTS TO SPOTS ----
    
    
    d <- raster::pointDistance(halo_obj_df[, c("Visium.Registered.X","Visium.Registered.Y")], cbind(X,Y), lonlat = FALSE)
    r <- apply(d, 1, which.min)
    # add spot call to the halo table
    halo_obj_df$Visium.Spot <- row.names(data_seurat@meta.data)[r]
    
    # now that we have spot labels, summarize halo table per spot
    # note collagen (blue) is stain 1 in results
    # note Tissue Area seems to be filtered for background (Region Area is not)
    var_list <- c("Tissue.Area.Analyzed..um.2.",
                  "Stain.1.Area..um.2.",
                  "Stain.1.Weak.Area..um.2.",
                  "Stain.1.Moderate.Area..um.2.",
                  "Stain.1.Strong.Area..um.2.")
    
    halo_spot_df <- halo_obj_df %>%
      dplyr::group_by(Visium.Spot, Classifier.Label) %>%
      dplyr::summarize_at(dplyr::all_of(var_list), sum)
    
    # Now how to treat Visium spots with no trichrome data?
    # true zeros need to be separate from spots lacking coreg tissue
    # HALO analysis modified to provide non-stromal objects (glass, tumor)
    # will need to zero out measures in tumor and glass (not in alg scope)
    # however keep Tumor area for collagen % over tissue measures
    halo_spot_df[halo_spot_df$Classifier.Label=='Tumor', 4:ncol(halo_spot_df)] <- 0
    halo_spot_df[halo_spot_df$Classifier.Label=='Glass', 3:ncol(halo_spot_df)] <- 0
    halo_spot_df$not_empty <- ifelse(halo_spot_df$Classifier.Label=="Glass", 0, 1)
    
    # compute prop stromal value (note tumor and glass will remain zero)
    halo_spot_df <- halo_spot_df %>%
      dplyr::mutate(Proportion.Stroma.Collagen.Positive = Stain.1.Area..um.2./Tissue.Area.Analyzed..um.2.) %>%
      dplyr::mutate(Proportion.Stroma.Collagen.Weak.Positive = Stain.1.Weak.Area..um.2./Tissue.Area.Analyzed..um.2.) %>%
      dplyr::mutate(Proportion.Stroma.Collagen.Moderate.Positive = Stain.1.Moderate.Area..um.2./Tissue.Area.Analyzed..um.2.) %>%
      dplyr::mutate(Proportion.Stroma.Collagen.Strong.Positive = Stain.1.Strong.Area..um.2./Tissue.Area.Analyzed..um.2.)
    
    # now Tumor and Glass have zero measures, drop classifier label and condense again
    # ignore NaN so glass values don't cancel out stromal prop values
    # use the not_empty column to drop spots that were entirely glass (no coreg tissue underneath)
    halo_spot_df <- halo_spot_df %>%
      dplyr::select(-Classifier.Label) %>%
      dplyr::group_by(Visium.Spot) %>%
      dplyr::summarize_if(is.numeric, ~sum(., na.rm = TRUE)) %>%
      dplyr::filter(not_empty > 0) %>%
      dplyr::select(-not_empty)
    
    # compute prop tissue (tumor+stroma) value
    halo_spot_df <- halo_spot_df %>%
      dplyr::mutate(Proportion.Tissue.Collagen.Positive = Stain.1.Area..um.2./Tissue.Area.Analyzed..um.2.) %>%
      dplyr::mutate(Proportion.Tissue.Collagen.Weak.Positive = Stain.1.Weak.Area..um.2./Tissue.Area.Analyzed..um.2.) %>%
      dplyr::mutate(Proportion.Tissue.Collagen.Moderate.Positive = Stain.1.Moderate.Area..um.2./Tissue.Area.Analyzed..um.2.) %>%
      dplyr::mutate(Proportion.Tissue.Collagen.Strong.Positive = Stain.1.Strong.Area..um.2./Tissue.Area.Analyzed..um.2.)
    
    # Add distance to the pathology feature as a meta data factor of the Seurat object
    trichrome_df <- halo_spot_df %>%
      dplyr::select(-Tissue.Area.Analyzed..um.2.) %>%
      tibble::column_to_rownames("Visium.Spot") %>%
      dplyr::rename_all(function(x) paste0("Trichrome.",x))
    
    # missing spots means no coreg tissue existed, should be NA in seurat
    data_seurat@meta.data <- merge(data_seurat@meta.data, trichrome_df, by="row.names", all.x = TRUE) %>%
      tibble::column_to_rownames("Row.names")
    
    
    ## VISUALIZE ----
    pathology_var_name <- "Trichrome.Proportion.Tissue.Collagen.Positive"
    p1 <- spatial_feature_plot_my(data_seurat, pathology_var_name, min.cutoff=0, max.cutoff=1, crop=TRUE)
    
    # Save distance visualization to file
    filename <- file.path(output_figs, paste0("Trichrome_", sample_name, "_", pathology_var_name))
    write_plot2file_my(p1, filename, num_row=1, num_col=1)
    
    
    ## OVERWRITE SEURAT FILE* ----
    # * or write to csv
    
    if(do_rds_output == TRUE){
      # Write the updated Seurat object
      filename <- file.path(output_path, paste0(sample_name, "_trichrome_pathology.rds"))
      saveRDS(data_seurat, file=filename)
    }
    if(do_csv_output == TRUE){
      # Write the distance data to a csv (only distance data)
      filename <- file.path(output_path, paste0(sample_name, "_trichrome_pathology.csv"))
      data_seurat@meta.data %>%
        dplyr::select(contains("Trichrome.")) %>%
        write.csv(file = filename, row.names = TRUE)
    }
    
  } else{
    warning("More than one file found!\n", file_list)
  } # end conditional check for singular file per sample
  
  
} # end looping through samples





