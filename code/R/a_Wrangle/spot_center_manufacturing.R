# Author: Hannah Pliner, Anna Lyubetskaya. Date: 23-06-27

# Compare spot center identification with and without a manufacturing file from 10x
# Pick the result with higher empirical center accuracy

# Note: If you get a permission denied from the system() line, chmod u+x code/utils/spotdetect.py


## ENVIRONMENT ----


# Allow piping throughout the package
`%>%` <- magrittr::`%>%`

# Check we have the right python packages and if not install
if (!reticulate::py_module_available('cv2')) {
  system('pip install opencv-python')
}

library("optparse")
source("code/utils/utils_in_out.R")
source("code/R/a_Wrangle/10X_to_Seurat_ST_utils.R")


## FIXED PARAMETERS ----


# Location of manufacturing files in GPR format
manufacturing_files_path <- "XXXX"
spaceranger_spatial_only_path <- "XXXX"


## PATHS ----


# Set default values for the two parameters
meta_file <- "XXXX"
input_object_path <- "XXXX"

# Output folders
output_path <- "XXXX"
output_folders <- create_output_subfolders_my(output_path, c("Figures_spotdetect", "Figures_spotdetect_man",
                                                             "raw_concentriq_images", "Manufacturing_out"))

# Log file to accumulate summary stats
log_file <- paste0(output_folders$Figures_spotdetect_man, "log_spotdetect.txt")


## INGEST, WRANGLE, OUTPUT ----


# Start the log file
write(paste(c("Sample_Name","Slide","Area","ConcentriqID",
              "Median_spot_center_distance_before","Median_spot_center_distance_after",
              "Spot_number_before","Spot_number_after","Difference","Pass"), collapse="\t"), log_file)

# Read sample file, filter for samples that pass QC
meta_df <- readr::read_delim(file=meta_file, delim="\t")


for(i in 1:nrow(meta_df)){
  
  ## Set parameters ----
  
  
  # Load the 10X Spatial folder and sample name - required parameters
  sample_name <- meta_df[[i, "Sample_Name"]]
  filename_in <- meta_df[[i, "Pipeline.FullPath"]]
  
  print(filename_in)
  
  # Path to a processed RDS object
  filename_rds <- paste0(input_object_path, sample_name, "_all.rds")  
  
  
  if(file.exists(filename_rds)){
    
    
    ## Run spotdetect script to measure spot placement error ---- 
    
    
    system(paste0('code/utils/spotdetect.py ', filename_in, ' -o ', output_folders$Figures_spotdetect, ' -s ', sample_name))
    
    # Calculate error metrics
    spotdetect_df <- gen_spotdetect_df_my(filename_in, 
                                          paste0(output_folders$Figures_spotdetect, "/", sample_name, "_spots_detected.csv"), 
                                          max_error = 50) 
    
    # Generate summary fig
    generate_spotdetect_summary_fig_my(spotdetect_df, sample_name, 
                                       paste0(output_folders$Figures_spotdetect, "/", sample_name,  "_spotdetect_error.png") )
    
    
    ## Rerun spaceranger spatial part to check for better placement with manufacturing files ----
    
    
    all_good <- TRUE
    
    # Check for required parameters
    if("Concentriq_Image_ID" %in% colnames(meta_df) & "Slide" %in% colnames(meta_df) & "Area" %in% colnames(meta_df)){
      concentriq_image <- meta_df[[i, "Concentriq_Image_ID"]]
      slideid <- meta_df[[i, "Slide"]]
      slidearea <- meta_df[[i, "Area"]]
    } else {
      all_good <- FALSE
      print("Missing key metadata for rerun_with_manufacturing option") 
    }
    
    if (!file.exists(paste0(manufacturing_files_path,"/", slideid, ".gpr"))) {
      all_good <- FALSE
      print("Missing manufacturing file for rerun_with_manufacturing option")
    }
    
    
    if(all_good) {
      
      # Download image from concentriq if it doesn't exist
      system(paste0(spaceranger_spatial_only_path, "/download_fig.py ", concentriq_image, " ", output_folders$raw_concentriq_image))
      
      # Rerun spaceranger spatial only with manufacturing file
      system(paste0(spaceranger_spatial_only_path, "/spaceranger count --fastqs ", spaceranger_spatial_only_path, "/fastqs/  --id ", sample_name, "  --transcriptome ", spaceranger_spatial_only_path, 
                    "/transcriptome/  --image ",output_folders$raw_concentriq_images, "/", concentriq_image, ".jpg --probe-set ", spaceranger_spatial_only_path, 
                    "/probes/probes.csv  --reorient-images true --no-probe-filter --slide ", slideid, "  --area ", slidearea, "1 --r1-length 26  --slidefile ",
                    manufacturing_files_path,"/", slideid, ".gpr --sample test"))
      
      system(paste0("cp -r ", sample_name, " ", output_folders$Manufacturing_out))
      
      new_input <- paste0(output_folders$Manufacturing_out, "/", sample_name, "/outs/")
      
      # Rerun spotdetect and compare output
      system(paste0('code/utils/spotdetect.py ', new_input, ' -o ', output_folders$Figures_spotdetect_man, ' -s ', sample_name, "_with_man"))
      
      # Calculate error metrics
      spotdetect_df_man <- gen_spotdetect_df_my(new_input,
                                                paste0(output_folders$Figures_spotdetect_man, "/", sample_name, "_with_man_spots_detected.csv"),
                                                max_error = 50)
      
      # Generate summary fig
      generate_spotdetect_summary_fig_my(spotdetect_df_man, paste0(sample_name, "_with_man"),
                                         paste0(output_folders$Figures_spotdetect_man, "/", sample_name,  "_spotdetect_error_with_man.png") )
      
      # Calculate stats
      med_before <- round(median(spotdetect_df$distance_um), 2)
      med_after <- round(median(spotdetect_df_man$distance_um), 2)
      
      # Create a log of progress
      write(paste(c(sample_name, slideid, slidearea, concentriq_image,
                    med_before, med_after,
                    nrow(spotdetect_df), nrow(spotdetect_df_man), 
                    med_before - med_after,
                    med_before > med_after), collapse="\t"), log_file, append = TRUE)
      
      if(median(spotdetect_df$distance_um) > median(spotdetect_df_man$distance_um)) {
        
      }
    }
    
  }
}