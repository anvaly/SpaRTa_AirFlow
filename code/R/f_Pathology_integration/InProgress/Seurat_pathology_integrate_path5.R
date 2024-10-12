# Author: Andrew Fisher. Date: 23-09-19
# Add pathology compartments (from PathAI) to a Seurat object

# This script takes in 3 vectors of pathology annotation name vectors:
# 1. expected_classes = a list of pathology classes to read in at the XML read stage;
# ---- path5 can contain a bunch of other annotations that will be ignored
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
cohort_name <- "IBD"

# All classes to process
expected_classes <- c('blood vessel',
                      'crypt abscess',
                      'erosion or ulceration',
                      'granulation tissue',
                      'infiltrated epithelium',
                      'inter gland lumen',
                      'lamina propria',
                      'muscularis mucosa',
                      'normal epithelium')

# Rename classes to new labels if necessary
expected_classes_renamed <- c('vessel',
                              'crypt_abscess',
                              'erosion_ulceration',
                              'granulation',
                              'infil_epi',
                              'gland_lumen',
                              'lamina_propria',
                              'muscularis_mucosa',
                              'normal_epi')

expected_cell_classes <- c("eosinophil",
                           "epithelial.non goblet cell enterocyte",
                           "goblet.cell nucleus",
                           "intraepithelial.lymphocyte",
                           "lymphocyte.non intraepithelial",
                           "neutrophil",
                           "other.cells",
                           "plasma.cell")

if(is.null(expected_classes_renamed)){
  expected_classes_renamed <- expected_classes
}
name_dict <- expected_classes_renamed
names(name_dict) <- expected_classes

# A list of classes that add up to all annotations we want to count towards the total surface of the spot
reference_classes <- c('vessel',
                       'crypt_abscess',
                       'erosion_ulceration',
                       'granulation',
                       'infil_epi',
                       'lamina_propria',
                       'muscularis_mucosa',
                       'normal_epi')

# A list of classes to count towards the dominant class label for each spot
# The order of this vector is important: it decides how to resolve ties between dominant classes
annotation_classes <- c('lamina_propria',
                        'muscularis_mucosa',
                        'granulation',
                        'normal_epi',
                        'infil_epi',
                        'erosion_ulceration',
                        'gland_lumen',
                        'crypt_abscess',
                        'vessel')

# Perform a subset of the Seurat object to only those spots that have a pathology annotation
# Subset necessitates re-clustering
do_subset <- FALSE

# Prefix to add to column names
col_prefix <- "AIM_UC."

# Don't write the final RDS object
no_rds_output <- FALSE

# Set seed for clustering
set.seed <- 428

# User defined colors
cols <- NULL


## PATHS ----


# Input file with sample meta data
meta_file <- "XXXX"
hif_file <- "XXXX"
path5_file <- "XXXX"

# Input folder
input_path <- "XXXX"

# Output folder
output_path <- "XXXX"

# Output folder
output_figs <- file.path(output_path, "figs")

# Create output folders
dir.create(output_path, showWarnings = FALSE)
dir.create(output_figs, showWarnings = FALSE)

filename_log <- file.path(output_figs, paste0("pathology_subset_log.txt"))
write(paste(c("Sample_Name", "SpotsBefore", "SpotsAfter", "SpotsUnderMedian", "SpotsUnderThirdMedian"), collapse="\t"), filename_log, append=TRUE)


## INGEST DATA ----


# Read sample file; select only cohort-relevant samples
manifest_df <- readr::read_csv(file=meta_file)
hif_df <- readxl::read_xlsx(hif_file, sheet = "FEATURE")
meta_df <- merge(manifest_df, hif_df, by = 'SAMPLE_ID')

# Read in the RData file from stash
SpaRTa_IBD_points_and_polygons <- readRDS(path5_file)

# loop through all samples
for (curr_sample in meta_df$SAMPLE_ID){
  
  # get matching seurat file
  seurat_filename <- gsub("\\.svs$", "_annotated.rds", meta_df$SVSFILENAME[meta_df$SAMPLE_ID==curr_sample])
  
  # Open a connection to the RDS object
  con <- gzfile(file.path(input_path, seurat_filename))
  
  # Ingest the Seurat object
  data_seurat <- readRDS(con)
  
  # Close the connection to be able to overwrite
  close(con)
  
  # Read in the polygons file (owin) file for the example PathAI ID
  pathai_file <- as.character(meta_df$CASE_ID[meta_df$SAMPLE_ID==curr_sample])
  
  # extract the regions
  polys_um <- SpaRTa_IBD_points_and_polygons[[pathai_file]]$polygons
  # Pipeline expects annotation dimensions in pixels, not microns
  polys <- lapply(polys_um, spatstat.geom::rescale, s=0.2529, unitname='pixel')
  
  # extract the cells
  pts_um <- SpaRTa_IBD_points_and_polygons[[pathai_file]]$points
  # Pipeline expects annotation dimensions in pixels, not microns
  pts <- spatstat.geom::rescale(pts_um, s=0.2529, unitname='pixel')
  
  ## PROCESS SPOTS
  # Seurat image object
  image_structure <- data_seurat@images[[1]]
  
  # Seurat spot center coordinates in full resolution: Y = imagerow, X = imagecol
  X <- image_structure@coordinates[["imagecol"]]
  Y <- image_structure@coordinates[["imagerow"]]
  
  # Spot diameter at full resolution
  spot_radius <- image_structure@scale.factors$spot_diameter_fullres / 2
  
  # Convert well matrix to spatstat ppp object
  spot_ppp <- spatstat.geom::as.ppp(cbind(X, Y), W = spatstat.geom::owin(range(X), range(Y)))
  
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
  write(paste(c(curr_sample, ncol(data_seurat), length(spots_remaining), total_under1_num, total_under2_num), collapse="\t"), 
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
  for(c in setdiff(paste0(col_prefix, expected_classes_renamed), colnames(data_subset_seurat@meta.data))){
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
    
    # Plot spot coverage by pathology
    p3 <- spatial_feature_plot_my(data_seurat, feature=total_col, min.cutoff="q0", max.cutoff="q100", title="")
    
    # Write figures to file
    filename <- file.path(output_figs, paste0("pathai_anno_", curr_sample))
    write_plot2file_my(patchwork::wrap_plots(list(p1, p3), nrow=1, ncol=2), filename, num_row=1, num_col=2)
    
    # Establish a color schema
    class_list <- intersect(gsub(col_prefix, "", colnames(data_seurat@meta.data)), expected_classes_renamed)
    
    # Pathology classification image reflecting full resolution of the annotation
    # Visually check that the classification and spots make sense
    filename <- file.path(output_figs, paste0("path_spat_overlay_", curr_sample, ".png"))
    plot_pathology_ingestion_check_my(polys, X, Y, spot_radius, filename, main_class="Tissue", 
                                      title=curr_sample, cols=cols, class_list=class_list)
    
    # Spatial representation of each pathology layer by spot
    filename <- file.path(output_figs, paste0("path_individ_class_", curr_sample))
    p <- batch_spatial_feature_plot_my(list(curr_sample = data_subset_seurat), paste0(col_prefix, pathology_classes_init), 
                                       output_file=filename, min.cutoff="q0", max.cutoff="q100")
  }
  
  
  ## PLOT PATHOLOGY VIEWS ----
  
  
  # Plot boxplots of spot pathology classes
  filename <- file.path(output_figs, paste0("path_class_prevalence_", curr_sample))
  p <- create_box_plot_my(clust_df, x_label="Pathology.Group", y_label="Percent", 
                          fill_label="Pathology.Group", filename=filename,
                          labels=c("Pathology Annotation", "Percent of spot surface", curr_sample), reorder_x=TRUE) 
  
  # Plot barplot of spot counts by pathology classes
  filename <- file.path(output_figs, paste0("path_perc_hist_", curr_sample))
  p <- create_bar_plot_my(clust_max_df %>% 
                            dplyr::group_by(Pathology.Group) %>% 
                            dplyr::summarise(Count = dplyr::n_distinct(Coordinate)), 
                          x_label="Pathology.Group", y_label="Count",
                          fill_label="Pathology.Group", filename=filename,
                          labels=c("Pathology Annotation", "Number of spots", curr_sample), reorder_x=TRUE)
  
  
  ## ADD CELL COUNTS
  
  # get fresh copy of metadata_df
  meta_seurat_df <- tibble::as_tibble(data_subset_seurat@meta.data, rownames = "Coordinate")
  
  # loop all cell types
  for (curr_cell_type in expected_cell_classes){
    # get cell type subset from pts ppp
    pts_split_phen <- split(pts, curr_cell_type, un = T)$`1`
    # compute spot compositions
    spot_cellcounts <- spatstat.core::crosspaircounts(spot_ppp, pts_split_phen, r = spot_radius)
    
    # add to Seurat
    meta_seurat_df[[paste0(col_prefix, "Cell.Cts.", curr_cell_type)]] <- spot_cellcounts
    
    # counts should be NA for spots not evaluated by PathAI (no predicted regions from above)
    meta_seurat_df[is.na(meta_seurat_df[[col_extra]]), paste0(col_prefix, "Cell.Cts.", curr_cell_type)] <- NA
  }
  
  # now bring data back into seurat object
  # Add the dominant cluster tag to the meta data tibble
  meta_seurat_df <- meta_seurat_df %>%
    tibble::column_to_rownames("Coordinate")
  
  # Add the dominant cluster tag to the Seurat object
  data_subset_seurat@meta.data <- meta_seurat_df
  
  
  
  ## OVERWRITE SEURAT FILE ----
  
  
  if(no_rds_output == FALSE){    
    # Write the updated Seurat object
    filename <- file.path(output_path, paste0(curr_sample, "_ann_path.rds"))
    saveRDS(data_subset_seurat, file=filename)
  }
  
  
}





