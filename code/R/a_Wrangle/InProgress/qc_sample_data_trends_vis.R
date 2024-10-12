# Author: Anna Lyubetskaya. Date: 21-02-20

# Visualize and correlate previously identified unspecific genes


## ENVIRONMENT ----


# Allow piping throughout the package
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")
source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_vis.R")


## PARAMETERS ----


# Cohort name
cohort_name <- "PDAC"


## PATHS ----


# A list of previously identified promiscuous genes
promiscuous_genes_file <- "XXXX"

# Input file with sample meta data
meta_file <- "XXXX"

# Output folder
output_path <- "XXXX"

# Create output folders
dir.create(output_path, showWarnings = FALSE)

# If probes were removed in the 10x pipeline, and you don't want to go back and re-run you can restore those probes from the raw data with this flag
return_probes <- FALSE
# File with probe offtarget information
probe_file <- "XXXX"


## INGEST DATA ----


# A list of previously identified unspecific probes
promiscuous_genes <- readr::read_delim(promiscuous_genes_file, delim="\t")

# Read sample file
cohort_meta_df <- readr::read_delim(file=meta_file, delim="\t")  %>% 
  dplyr::filter(grepl(cohort_name, Sample_Name) & Quality_Status == TRUE)


## RAW DATA CHECK ----


# Collect information into a list
# expression_list <- list()

for(i in c(1:nrow(cohort_meta_df))){
  
  
  ## INGEST 10X DATA ----
  
  
  # Input data path
  filename_in <- cohort_meta_df[[i, "FullPath"]]
  # Sample name
  sample_name <- cohort_meta_df[[i, "Sample_Name"]]
  
  # Ingest raw (unfiltered) data
  data_seurat_raw <- read_10X_spatial_folder_my(filename_in, filename="raw_feature_bc_matrix.h5")
  # Ingest filtered data
  data_seurat <- read_10X_spatial_folder_my(filename_in)
  
  # Update filtered data if probes were filtered out in the pipeline
  if(return_probes == TRUE){
    data_seurat <- restore_filtered_probes(data_seurat_raw, data_seurat, probe_file)
  }
  
  # Ingest full tissue positions list (because Seurat filters it)
  tissue_position_list <- readr::read_delim(paste0(filename_in, "/spatial/tissue_positions_list.csv"), 
                                            delim=",", col_names = c("tissue", "id", "row", "col", "imagerow", "imagecol"))
  
  # Put a full image into the raw Seurat data object instead of the filtered one
  data_seurat_raw@images$slice1@coordinates <- tissue_position_list %>%
    tibble::column_to_rownames("tissue") %>%
    dplyr::rename(tissue = id)
  
  # Fix the issue with an extra coordinates in the image slot
  cells <- intersect(rownames(data_seurat_raw@images$slice1@coordinates), colnames(data_seurat_raw))
  data_seurat_raw@images$slice1@coordinates <- data_seurat_raw@images$slice1@coordinates[cells,]
  
  
  ## EXCLUDE TISSUE ADJACENT SPOTS ----
  
  
  # 10X filtered spot list
  spot_under_tissue <- colnames(data_seurat)
  
  # Extract all Seurat coordinates
  xy <- tissue_position_list[, c("row", "col")]
  barcode_list <- tissue_position_list$tissue
  xy <- mapply(xy, FUN=as.numeric)
  
  # Adjacency distance
  distance <- 2
  
  # Find the matrix with the indices of points belonging to the set of the k nearest neighbours of each other
  spot_nb <- spdep::dnearneigh(xy, 0, distance, row.names=barcode_list)
  names(spot_nb) <- barcode_list
  
  # Find all spots within or adjacent to tissue
  spot_under_adjacent_tissue <- unique(unlist(unname(sapply(spot_nb[spot_under_tissue], function(x) barcode_list[x]))))
  # Check that all spots under tissue got included in the neighborhood analysis
  setdiff(spot_under_tissue, spot_under_adjacent_tissue)
  
  # Spots outside of tissue
  spot_outside_tissue <- setdiff(colnames(data_seurat_raw), spot_under_adjacent_tissue)
  # Spots adjacent to tissue
  spot_adjacent_tissue <- setdiff(spot_under_adjacent_tissue, spot_under_tissue)
  
  
  ## ANALYZE COORDINATE PROPERTIES ----
  
  
  # Expression long tibble
  data_meta_df <- tibble::as_tibble(data_seurat_raw@meta.data, rownames = "Coordinate")
  
  # UMI sum per spot, including filtered spots
  coordinate_df <- data_meta_df %>%
    dplyr::mutate(UnderTissue = Coordinate,
                  nCount_Spatial_log10 = log10(nCount_Spatial)) %>%
    dplyr::select(-orig.ident)
  
  # Assign spots into 3 buckets - under / adjacent / outside of tissue
  coordinate_df[which(coordinate_df$Coordinate %in% spot_under_tissue), "UnderTissue"] <- "UnderTissue"
  coordinate_df[which(coordinate_df$Coordinate %in% spot_adjacent_tissue), "UnderTissue"] <- "AdjacentTissue"
  coordinate_df[which(coordinate_df$Coordinate %in% spot_outside_tissue), "UnderTissue"] <- "OutsideTissue"
  
  # Factorize the tissue localization label
  coordinate_df["UnderTissue"] <- factor(coordinate_df[["UnderTissue"]], levels=c("UnderTissue", "AdjacentTissue", "OutsideTissue"))
  
  
  ## EXTRACT RELEVANT EXPRESSION DATA ----
  
  
  # Extract expression matrix, transform it into a long tibble, and calculate z-score for each gene across all barcodes
  expression_df <- seurat_expression_to_long_tibble_my(data_seurat_raw, assay="Spatial", slot="data") %>%
    dplyr::filter(Symbol %in% promiscuous_genes$Symbol) %>%
    dplyr::mutate(Spatial_data_log2 = round(log2(Spatial_data + 1), 2),
                  Sample = sample_name) %>%
    dplyr::inner_join(coordinate_df, by="Coordinate")

  # expression_list[[sample_name]] <- expression_df

  
  ## PLOT SPATIAL UMI DISTRIBUTIONS ----
  
  
  filename <- paste0(output_path, "box_promiscuous_", sample_name)
  create_box_plot_my(expression_df, x_label="Symbol", y_label="Spatial_data_log2", fill_label="Symbol",
                     facet_var=c("UnderTissue", "fixed"), filename=filename)
  
  filename <- paste0(output_path, "spatial_promiscuous_", sample_name)
  batch_spatial_feature_plot_my(list(object = data_seurat_raw), promiscuous_genes$Symbol, output_file=filename, 
                                min.cutoff="q1", max.cutoff="q99")


  ## CORRELATE GENES ----
  
  
  
}
