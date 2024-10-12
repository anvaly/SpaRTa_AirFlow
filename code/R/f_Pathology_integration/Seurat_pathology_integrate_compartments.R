# Author: Anna Lyubetskaya. Date: 21-01-18
# Add pathology compartments to a Seurat object

# This script takes in 3 vectors of pathology annotation name vectors:
# 1. expected_classes = a list of pathology classes to read in at the XML read stage;
# ---- XML can contain a bunch of other annotations that will be ignored
# 2. reference_classes = a list of pathology classes that add up to "tissue" for % calculation;
# ---- percent coverage will be calculated for every ingested pathology annotation against the sum of these pathology classes
# 3. annotation_classes = a list of pathology classes in a specific order to decide on a single categorical label.
# ---- the order is very important, if two classes tie for % coverage, the first one on this list will be picked


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

library(Seurat)

source("code/utils/utils_ggplot.R")
source("code/utils/utils_tibble.R")

source("code/R/Utils/utils_pathology_halo.R")
source("code/R/Utils/utils_10X_image.R")


## PARAMETERS ----


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

# A list of classes that add up to all annotations we want to count towards the total surface of the spot
reference_classes <- c("BenignEpi", "Blood", "ExoEndo", "LuminalNec", "NonEpi", "Tumor")

# A list of classes to count towards the dominant class label for each spot
# The order of this vector is important: it decides how to resolve ties between dominant classes
annotation_classes <- c("IntestineAdj", "MuscleAdj", "NormalAdj", 
                        "LymphNode", "TLSAggregate", "TLSImmature", "TLSMature",
                        "Adipose", "Vessel", "Muscle", "Nerve", "NonEpi", 
                        "BenignEpi", "Blood", "ExoEndo", "LuminalNec", "Tumor")

# Perform a subset of the Seurat object to only those spots that have a pathology annotation
# Subset necessitates re-clustering
do_subset <- TRUE

# Prefix to add to column names
col_prefix <- "Pathology."

# Don't write the final RDS object
no_rds_output <- FALSE

# Set seed for clustering
set.seed <- 531

# User defined colors
cols <- NULL


## PATHS ----


# Input file with sample meta data
meta_file <- "XXXX"

# Input folder
input_path <- "XXXX"

# Output folder
output_path <- "XXXX"

# Output folder
output_init_figs <- "XXXX"

# Create output folders
dir.create(output_path, showWarnings = FALSE)
dir.create(output_init_figs, showWarnings = FALSE)


## INGEST DATA ----


# Read sample file; select only cohort-relevant samples
meta_df <- readr::read_delim(file=meta_file, delim="\t") %>%
  dplyr::select(Sample_Name, Pathology_Compartment_File) %>%
  tidyr::drop_na()


## PERFROM XML PRE-CHECK ----


# Collect all XML region names
annot_list <- list()
for(xml_file in meta_df$Pathology_Compartment_File){
  annot_list[[xml_file]] <- xml_inspect_regions_my(xml_file)
}

# List of all pathology classes and their frequencies
path_dictionary <- tibble::as_tibble(as.matrix(sort(table(c(unname(unlist(annot_list)))))), rownames="PathClass") %>%
  dplyr::rename(Frequency = V1)

filename <- paste0(output_init_figs, "pathology_classes.txt")
readr::write_delim(path_dictionary, filename, delim="\t")


## SET COLORS ----


if(is.null(cols)){
  num_classes <- length(expected_classes_renamed)
  cols <- define_cols_my(n=num_classes, col_type="jet")
  names(cols) <- sort(expected_classes_renamed)      
  
  # Add special classes: Tissue and None
  cols[["Tissue"]] <- "black"
  cols[["None"]] <- "white"
}


## WRANGLE DATA ----


filename_log <- paste0(output_init_figs, "pathology_subset_log.txt")
write(paste(c("Sample_Name", "SpotsBefore", "SpotsAfter", "SpotsUnderMedian", "SpotsUnderThirdMedian"), collapse="\t"), filename_log, append=TRUE)

for(i in 1:nrow(meta_df)){
  
  gc()
  
  # Sample name
  sample_name <- meta_df[[i, "Sample_Name"]]
  # Pathology file to integrate
  path_file <- meta_df[[i, "Pathology_Compartment_File"]]
  
  # Find the RDS object
  file_list <- dir(input_path, pattern=paste0(sample_name), full.names=TRUE, recursive=TRUE)
  
  
  if(length(file_list) == 1){
    
    # Create the output folder for all figures
    output_figs <- paste0(output_init_figs, sample_name, "/")
    dir.create(output_figs, showWarnings = FALSE)
    
    print(file_list)
    print(path_file)
    
    
    filename_check <- paste0(output_figs, "path_perc_box_", sample_name)
    
    if(!file.exists(paste0(filename_check, ".png")) || 
       (!file.exists(paste0(output_path, sample_name, "_ann_path.rds")) && no_rds_output == FALSE)){
      
      
      ## INGEST SEURAT DATA ----
      
      
      # Open a connection to the RDS object
      con <- gzfile(file_list[1])
      
      # Ingest the Seurat object
      data_seurat <- readRDS(con)
      
      # Close the connection to be able to overwrite
      close(con)
      
      
      if(!"Pathology.Group" %in% colnames(data_seurat@meta.data)){
        
        
        ## EXTRACT SEURAT SPOT DATA ----
        
        
        # Seurat image object
        image_structure <- data_seurat@images[[1]]
        
        # Seurat spot center coordinates in full resolution: Y = imagerow, X = imagecol
        X <- image_structure@coordinates[["imagecol"]]
        Y <- image_structure@coordinates[["imagerow"]]
        
        # Spot diameter at full resolution
        spot_radius <- image_structure@scale.factors$spot_diameter_fullres / 2
        
        
        ## LOAD ANNOTATIONS ----
        
        
        # Load the annotations as a list of owin-s
        # Only keep specified classes for analysis (ignoring extraneous annotation)
        polys <- halo_annot_to_poly_my(path_file, region_list=expected_classes)
        
        # Merge polys subsets in case the user defined re-naming schema requires it
        expected_classes_renamed_table <- table(expected_classes_renamed)
        expected_classes_renamed_count <- names(expected_classes_renamed_table[expected_classes_renamed_table > 1])
        if(!is.null(expected_classes_renamed_count)){
          
          for(cl in expected_classes_renamed_count){
            p <- NULL
            
            for(cl_rep in expected_classes[grep(cl, expected_classes_renamed)]){
              if(is.null(p)){
                p <- polys[[cl_rep]]
              } else{
                p <- spatstat.geom::union.owin(p, polys[[cl_rep]])
              }
              
              # Remove the poly class that has just be merged into the union
              polys <- polys[names(polys) != cl_rep]
            }
            
            # Add the new union of polys to the list under the last class name
            polys[[cl_rep]] <- p
          }
          
        }
        
        # Rename classes if user-defined
        if(!is.null(expected_classes_renamed)){
          names(polys) <- unname(name_dict[names(polys)])
        }
        
        # Get list of available classes
        pathology_classes_init <- sort(names(polys))
        print(pathology_classes_init)
        
        
        ## ADD PATHOLOGY OVERLAPS TO SEURAT OBJECT ----
        
        
        # Compute overlap of each spot with each pathology annotation region
        # This function is parallelized for speed
        spot_data_list <- list()
        for(cl in pathology_classes_init){
          print(paste0("Overlap calculation: ", cl))
          spot_data_list[[cl]] <- spot_annotation_overlap_my(X, Y, spot_radius, polys[[cl]])
        }
        
        # Add overlap annotations to the Seurat object
        for(cl in pathology_classes_init){
          col_name <- paste0(col_prefix, cl)
          data_seurat@meta.data[[col_name]] <- round(spot_data_list[[cl]])
        }
        
        
        ## CALCULATE SPOT COVERAGE FOR EACH ANNOTATION ----
        
        
        # Sum across reference classes; this number will be used for subsetting to include spots that received no pathology annotation
        total_col <- paste0(col_prefix, "Total")
        col_names <- paste0(col_prefix, reference_classes)
        data_seurat@meta.data[[total_col]] <- unname(Matrix::rowSums(data_seurat@meta.data[intersect(col_names, colnames(data_seurat@meta.data))]))
        
        # Number of spots with less than full coverage
        total_coverage <- data_seurat@meta.data[[total_col]]
        total_under1_num <- length(which(total_coverage < mean(total_coverage)))
        
        # Number of spots with less than 1/3rd coverage
        subset_threshold <- mean(total_coverage) / 3
        total_under2_num <- length(which(total_coverage < subset_threshold))
        
        # Select spots characterized by pathology using a user-defined pixel threshold
        spots_remaining <- rownames(data_seurat@meta.data)[which(data_seurat[[total_col]] > subset_threshold)]
        
        # Write a subsetting log to a file
        write(paste(c(sample_name, ncol(data_seurat), length(spots_remaining), total_under1_num, total_under2_num), collapse="\t"), 
              filename_log, append=TRUE)
        
        
        # Calculate the ratio of each class relative to the spot surface in pixels
        for(cl in pathology_classes_init){
          col_name <- paste0(col_prefix, cl, ".percent")
          data_seurat@meta.data[[col_name]] <- round(spot_data_list[[cl]] / data_seurat@meta.data[[total_col]] * 100)
          data_seurat@meta.data[which(data_seurat@meta.data[[col_name]] > 100), col_name] <- 100
          
          print(paste(c(cl, length(which(data_seurat@meta.data[[col_name]] > 0)))))
        }
        
        # Identify classes for spots about to be removed
        spots_removed <- rownames(data_seurat@meta.data)[which(total_coverage < subset_threshold)]
        spots_removed_info <- data_seurat@meta.data[spots_removed, paste0(col_prefix, pathology_classes_init, ".percent")]
        
        
        ## SUBSET SEURAT OBJECT IF NECESSARY ----
        
        
        if(do_subset == TRUE && length(spots_remaining) < ncol(data_seurat)){
          
          cat("Subsetting", length(spots_remaining), "spots of", ncol(data_seurat), "\n")
          
          # Subset Seurat object
          data_subset_seurat <- subset(data_seurat, cells = spots_remaining)
          
          # Remove any floating tissue spots using the contiguity filter
          data_subset_seurat <- tissue_contiguity_filter_my(data_subset_seurat, spot_num=2)
          
        } else{
          data_subset_seurat <- data_seurat
          spots_removed <- NULL
          
          # Reference column
          total_col <- paste0(col_prefix, "Total")
        }
        
        # Add any missing pathology columns - this would happen when a sample doesn't have one of the annotation layers
        for(c in setdiff(paste0("Pathology.", expected_classes_renamed), colnames(data_subset_seurat@meta.data))){
          data_subset_seurat@meta.data[[c]] <- 0
          data_subset_seurat@meta.data[[paste0(c, ".percent")]] <- 0
          
          data_seurat@meta.data[[c]] <- 0
          data_seurat@meta.data[[paste0(c, ".percent")]] <- 0
        }
        
        
        ## ADD PATHOLOGY DOMINANT CLASS TO SEURAT OBJECT ----
        
        
        # Seurat meta-data column names storing pathology class percentages
        # Only keep specified classes for clustering analysis (ignoring extraneous annotation)
        # Other class data will remain in the Seurat object but won't be counted towards the Pathology.Group
        path_cols <- intersect(colnames(data_subset_seurat@meta.data), paste0(col_prefix, annotation_classes, ".percent"))
        
        # Seurat meta data
        meta_seurat_df <- tibble::as_tibble(data_subset_seurat@meta.data, rownames="Coordinate")
        
        # Add a column representing "no pathology classification found for a given spot"
        col_extra <- paste0(col_prefix, "None.percent")
        meta_seurat_df[[col_extra]] <- 100 - unname(rowSums(meta_seurat_df[path_cols]))
        meta_seurat_df[which(meta_seurat_df[[col_extra]] < 0), col_extra] <- 0
        
        # Extract pathology classes and Seurat clustering data
        clust_df <- meta_seurat_df %>%
          dplyr::select(dplyr::all_of(c("Coordinate", path_cols, col_extra))) %>%
          df_wide2long_my(key="Pathology.Group", val="Percent", start_col=2) %>%
          dplyr::mutate(Pathology.Group = gsub(paste0(col_prefix, "|.percent"), "", Pathology.Group))
        
        # Find the most abundant pathology class for each spot
        clust_max_df <- clust_df %>%
          dplyr::group_by(Coordinate) %>%
          dplyr::arrange(desc(Percent), factor(Pathology.Group, levels=annotation_classes)) %>%
          dplyr::slice(1) %>%
          dplyr::ungroup() %>%
          dplyr::arrange(Pathology.Group, Percent)
        
        # Add the dominant cluster tag to the meta data tibble
        meta_seurat_df <- meta_seurat_df %>%
          dplyr::left_join(clust_max_df %>%
                             dplyr::select(Coordinate, Pathology.Group), by="Coordinate") %>%
          tibble::column_to_rownames("Coordinate") #%>%
        
        # Add the dominant cluster tag to the Seurat object
        data_subset_seurat@meta.data <- meta_seurat_df
        colnames(data_subset_seurat@meta.data) <- gsub("Pathology.Group", paste0(col_prefix, "Group"), 
                                                       colnames(data_subset_seurat@meta.data))
        
        
        ## VISUALIZE PATHOLOGY CLASSES AND SUBSET SPOTS ----
        
        
        if(length(spots_removed) > 0){
          
          # Create a spatial plot of the tissue with no overlays
          p1 <- Seurat::SpatialDimPlot(data_seurat, pt.size.factor = 0) + 
            ggplot2::theme(legend.position = "none",
                           plot.title = element_text(size=7)) +
            ggplot2::labs(x="", y="", title=paste0("Spot number = ", length(colnames(data_seurat))))
          
          # Create a spatial plot highlighting removed spots
          p2 <- spatial_dim_plot_my(data_seurat, cell_highlight=spots_removed,
                                    title=paste0("Spots removed = ", total_under2_num)) + 
            ggplot2::theme(legend.position = "none")
          
          # Plot spot coverage by pathology
          p3 <- spatial_feature_plot_my(data_seurat, feature=total_col, min.cutoff="q0", max.cutoff="q100", title="")
          
          # Write figures to file
          filename <- paste0(output_figs, "spots_removed_", sample_name)
          write_plot2file_my(patchwork::wrap_plots(list(p1, p2, p3), nrow=1, ncol=3), filename, num_row=1, num_col=3)
          
          
          # Establish a color schema
          class_list <- intersect(gsub(col_prefix, "", colnames(data_seurat@meta.data)), expected_classes_renamed)
          
          # Pathology classification image reflecting full resolution of the annotation
          # Visually check that the classification and spots make sense
          filename <- paste0(output_figs, "path_spat_overlay_", sample_name, ".png")
          plot_pathology_ingestion_check_my(polys, X, Y, spot_radius, filename, main_class="Tissue", 
                                            title=sample_name, cols=cols, class_list=class_list)
          
          # Spatial representation of each pathology layer by spot
          filename <- paste0(output_figs, "path_individ_class_", sample_name)
          p <- batch_spatial_feature_plot_my(list(sample_name = data_subset_seurat), paste0(col_prefix, pathology_classes_init), 
                                             output_file=filename, min.cutoff="q0", max.cutoff="q100")
        }
        
        
        ## PLOT PATHOLOGY VIEWS ----
        
        
        # Plot boxplots of spot pathology classes
        p <- create_box_plot_my(clust_df, x_label="Pathology.Group", y_label="Percent", 
                                fill_label="Pathology.Group", filename=filename_check,
                                labels=c("Pathology Annotation", "Percent of spot surface", sample_name), reorder_x=TRUE) 
        
        # Plot barplot of spot counts by pathology classes
        filename <- paste0(output_figs, "path_perc_hist_", sample_name)
        p <- create_bar_plot_my(clust_max_df %>% 
                                  dplyr::group_by(Pathology.Group) %>% 
                                  dplyr::summarise(Count = dplyr::n_distinct(Coordinate)), 
                                x_label="Pathology.Group", y_label="Count",
                                fill_label="Pathology.Group", filename=filename,
                                labels=c("Pathology Annotation", "Number of spots", sample_name), reorder_x=TRUE)
        
        
        ## OVERWRITE SEURAT FILE ----
        
        
        if(no_rds_output == FALSE){    
          # Write the updated Seurat object
          filename <- paste0(output_path, sample_name, "_ann_path.rds")
          saveRDS(data_subset_seurat, file=filename)
        }
        
      } else{
        print("Skipping because this samples already has the Pathology.Group column!")
      }
    } else{
      print("Skipping because this samples is already analyzed!")
    }
  } else{
    print("More than one file found!")
  }
}
