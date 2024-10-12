# Author: Andrew Fisher. Date: Aug 3, 2023
# Extract 10X spots as a HALO annotation
# Pathology XML template
# xml_template_path <- "XXXX"


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

library(Seurat)
library(xml2)

source("./code/R/Utils/utils_10X_signatures.R")
source("./code/utils/utils_signatures.R")

## PARAMETERS ----


# Sample / Cohort name
cohort_name <- "PDAC_Naive_Full_Spots"
cohort_str <- "PDAC"

# define subset of samples to run or NULL for full set
sample_list <- NULL

# Select a specific sample in case of an integrated cohort
field_split <- "user.Sample_Name"  # NULL

# Partition XML into annotation layer using the following variable
# "Classifier_PDAC_Epi" or "user.Sample_Name" or NULL
region_name_value <- NULL  # "user.Sample_Name"

# Signature Params
sig_select <- NULL # c("PDAC.U.Stroma.activated")
# sct_threshold <- 0.25 # minimum SCT value
# spot_threshold <- 5 # how many spots need to pass minimum
assay = 'SCT'
# slot = 'data'

# Spot parameters
spot_grow_factor <- 1 # multiplier to spot radius for creating training spot from Visium spot


## PATHS ----

# Input file with sample meta data
meta_file <- "XXXX"

# Input rds folder
input_path <- "XXXX"

# Output folder
output_path <- "XXXX"

# Path to a signature file
sig_path <- "/repos/P02567_TBIO-3021_10X_pilot/data/import/Signatures/signatures_220210.txt"

# Create output folders
dir.create(output_path, showWarnings = FALSE)

# For PDAC classes: c("user.Sample_Name", "Clustering_Preferred", "sig.collisson.pdac.classical.EPvalue", "Pathology.Epithelium.percent", "Classifier_PDAC_Epi")
if (is.null(sig_select)){
  meta_data_cols <- c("user.Sample_Name")
} else {
  meta_data_cols <- c("Pathology.Epithelium.percent", "user.Sample_Name", paste0("sig.",sig_select))
}


## DEFINE XML ----


# Form the XML scaffold
xml_start <- "<Annotations>\n"
xml_head_init <- "<Annotation LineColor=\"16776960\" Name=\"ANNOTATION_NAME\" Visible=\"False\">\n<Regions>"
xml_node_init <- "<Region Type=\"Ellipse\" HasEndcaps=\"0\" NegativeROA=\"0\">\n<Vertices>\n<V X=\"QQ1000QQ\" Y=\"QQ2000QQ\" />\n<V X=\"QQ3000QQ\" Y=\"QQ4000QQ\" />\n</Vertices>\n<Comments>\n
<Comment Author=\"Anna.Lyubetskaya@bms.com\" CreatedTime=\"2021-11-02T12:00:00-04:00\" ModifiedTime=\"2021-11-02T12:00:00-04:00\" Body=\"Barcode = BARCODE_ID\" />\n"
xml_tail_init <- "</Regions>\n</Annotation>\n"
xml_stop <- "</Annotations>\n"

# Add additional meta data in the comments section
if(!is.null(meta_data_cols)){
  for(m in meta_data_cols){
    xml_node_init <- paste0(xml_node_init, "<Comment Author=\"Anna.Lyubetskaya@bms.com\" CreatedTime=\"2021-11-02T12:00:00-04:00\" ModifiedTime=\"2021-11-02T12:00:00-04:00\" Body=\"", m," = ", toupper(m),"\" />\n")
  }
}
xml_node_init <- paste0(xml_node_init, "</Comments>\n</Region>")


## INGEST DATA ----

# Read sample file; select only cohort-relevant samples
meta_df <- readr::read_delim(file=meta_file, delim="\t") %>%
  dplyr::filter(grepl(cohort_str, Sample_Name) &
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
    
    # Initiate a list of spot sizes for each sample
    test_size_list <- list()
    
    # Calculate Signatures
    if(is.null(sig_select) == FALSE){
      
      # Find genes abundant in this sample
      gene_list <- data_seurat@assays$SCT@var.features 
      # Load signatures and filter them down to only well represented genes
      signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select,
                                                  sig_length_min=1, sig_length_max=1000, ratio_threshold=0)
      # Add signature scores to a seurat object
      data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.", assay=assay)
      
      # Seurat renames column names, this step finds new signature names
      sig_names_upd <- colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))])
    }
    
    ## EXTRACT SEURAT SPOT DATA ----
    
    
    # Iterate through images if necessary
    XY_list <- list()
    for(s in names(data_seurat@images)){
      # Seurat image object
      image_structure <- data_seurat@images[[s]]
      
      # Spot diameter at full resolution
      spot_radius <- image_structure@scale.factors$spot_diameter_fullres / 2
      
      # Seurat spot center coordinates in full resolution: Y = imagerow, X = imagecol
      XY_list[[s]] <- image_structure@coordinates[c("imagecol", "imagerow")] %>%
        dplyr::rename(X = imagecol, Y = imagerow) %>%
        dplyr::mutate(spot_radius = spot_radius) %>%
        tibble::rownames_to_column("Coordinate")
    }
    
    # A tibble of image coordinates and spot radii
    XY <- dplyr::bind_rows(XY_list)
    
    
    ## EXTRACT META DATA ----
    
    
    if(!is.null(meta_data_cols) || !is.null(field_split)){
      # Extract meta data request by the user
      seurat_meta_df <- data_seurat@meta.data %>%
        tibble::rownames_to_column("Coordinate") %>%
        dplyr::select(dplyr::all_of(c("Coordinate", meta_data_cols, field_split)))
      
      # Join meta data with coordinate information
      XY <- XY %>%
        dplyr::inner_join(seurat_meta_df, by="Coordinate") %>%
        tibble::column_to_rownames("Coordinate")
    }
    
    if(!is.null(region_name_value)){
      region_name_list <- unique(seurat_meta_df[[region_name_value]])
      curr_region_name <- region_name_value
    } else{
      curr_region_name <- "region_name_value"
      region_name_list <- c(sample_name)
      XY[[curr_region_name]] <- sample_name
    }
    
    # Split output by a meta data field if user-defined
    if(!is.null(field_split)){
      field_split_list <- unique(seurat_meta_df[[field_split]])
    } else{
      field_split <- "field_split"
      field_split_list <- c(sample_name)
      XY[[field_split]] <- sample_name
    }
    
    
    ## APPLY GROUPING TO METADATA ----
    # # drop out the majority epithelium spots (otherwise it biases classes)
    # XY <- XY %>%
    #   dplyr::filter(Pathology.Epithelium.percent < 50)
    # # currently set to 25% and 75% quartiles
    # for (sig.col in paste0("sig.",sig_select)){
    #   cutoffs <- quantile(XY[,sig.col], prob = c(0.25, 0.75), na.rm = TRUE)
    #   XY[[paste0(sig.col,".HIGH")]] <- XY[[sig.col]] >= cutoffs[2]
    #   XY[[paste0(sig.col,".LOW")]] <- XY[[sig.col]] <= cutoffs[1]
    # }
    
    # Update the region_name_list to split on groupings
    if(!is.null(sig_select)){
      region_name_tab <- expand.grid(region_name_list,
                                      "..",
                                      paste0("sig.",sig_select),
                                      c(".HIGH",".LOW"),
                                      KEEP.OUT.ATTRS = TRUE,
                                      stringsAsFactors = FALSE)
    } else {
      region_name_tab <- expand.grid(region_name_list,"","","",
                                      KEEP.OUT.ATTRS = TRUE,
                                      stringsAsFactors = FALSE)
    }
    
    
    
    
    ## CREATE HALO ANNOTATE FILE ---- 
    
    # Iterate over field values that define the output file
    for(f_out in field_split_list){
      
      # Iterate over desired regions in the final XML
      region_list <- c()
      for(region_row in 1:nrow(region_name_tab)){
        
        # Subset the XY tibble with spot coordinates and select meta data
        if(!is.null(sig_select)){
          XY_loc <- XY %>%
            dplyr::filter(!!rlang::sym(field_split) == f_out & 
                            !!rlang::sym(curr_region_name) == region_name_tab[region_row,1] &
                            !!rlang::sym(paste0(region_name_tab[region_row,3:4], collapse = "")) == TRUE)
        } else {
          XY_loc <- XY %>%
            dplyr::filter(!!rlang::sym(field_split) == f_out & 
                            !!rlang::sym(curr_region_name) == region_name_tab[region_row,1])
        }
        
        # Initiate all components of the section
        xml_head <- xml_head_init
        xml_node <- xml_node_init
        xml_tail <- xml_tail_init
        
        # Add image name to the XML
        xml_head <- gsub("ANNOTATION_NAME", paste0(region_name_tab[region_row,], collapse = ""), xml_head)
        
        # We are approximating spots as ellipses in HALO annotations. Ellipse shape is defined by bounding box corners.
        xml_node_list <- list()
        for(spot in rownames(XY_loc)){
          xml_node_loc <- xml_node
          
          # Spot X and Y coordinates
          X_loc <- XY_loc[[spot, "X"]]
          Y_loc <- XY_loc[[spot, "Y"]]
          spot_radius <- XY_loc[[spot, "spot_radius"]] * spot_grow_factor
          
          # Find bounding box corners and insert them into the XML
          xml_node_loc <- gsub("\\bQQ1000QQ\\b", round(X_loc - spot_radius), xml_node_loc)
          xml_node_loc <- gsub("\\bQQ2000QQ\\b", round(Y_loc - spot_radius), xml_node_loc)
          xml_node_loc <- gsub("\\bQQ3000QQ\\b", round(X_loc + spot_radius), xml_node_loc)
          xml_node_loc <- gsub("\\bQQ4000QQ\\b", round(Y_loc + spot_radius), xml_node_loc)
          xml_node_loc <- gsub("\\bBARCODE_ID\\b", spot, xml_node_loc)
          
          # Add user-requested meta data
          if(!is.null(meta_data_cols)){
            for(m in meta_data_cols){
              xml_node_loc <- gsub(paste0("\\b", toupper(m),"\\b"), XY_loc[[spot, m]], xml_node_loc)
            }
          }
          
          xml_node_list[[spot]] <- xml_node_loc
          
          test_size_list[[spot]] <- (round(X_loc + spot_radius) - round(X_loc - spot_radius)) * (round(Y_loc + spot_radius) - round(Y_loc - spot_radius))
        }
        
        region_list <- c(region_list, c(xml_head, unlist(unname(xml_node_list)), xml_tail))
      }
      
      # Combine all part of the XML
      xml_file <- paste(c(xml_start, region_list, xml_stop), collapse="\n")
      
      # Write to file
      filename <- file.path(output_path, paste0(f_out, ".annotations"))
      write(xml_file, filename)
      
      # Read and write the file again using XML2 library to adjust formatting
      xml_data <- xml2::read_xml(filename)
      xml2::write_xml(xml_data, filename)
    }
  } # conditional check for single rds file per sample
} # loop through meta_df
