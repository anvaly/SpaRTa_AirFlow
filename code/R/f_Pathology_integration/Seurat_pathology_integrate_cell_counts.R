# Author: Anna Lyubetskaya. Date: 21-03-07
# Add pathology cell counts to a Seurat object


## SETUP ENVIRONMENT ----


# Allow piping throughout the package
`%>%` <- magrittr::`%>%`

library(Seurat)

source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_pathology_halo.R")
source("code/R/Utils/utils_10X_vis.R")


## PARAMETERS ----


# Sample / Cohort name
cohort_name <- "PDAC"

# Prefix to add to column names
col_prefix <- "CellCounts"

# Don't write the final RDS object
no_rds_output <- FALSE


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
  dplyr::select(Sample_Name, Pathology_CellCount_File) %>%
  tidyr::drop_na()


## WRANGLE DATA ----


for(i in 1:nrow(meta_df)){
  
  print(i)
  
  # Sample name
  sample_name <- meta_df[[i, "Sample_Name"]]
  # Pathology file to integrate
  path_file <- meta_df[[i, "Pathology_CellCount_File"]]

  # Find the RDS object
  file_list <- dir(input_path, pattern=paste0(sample_name), full.names=TRUE, recursive=TRUE)
  
  filename_check <- paste0(output_path, sample_name, "_annotated_pathology.rds")
  
  if(length(file_list) == 1 & !file.exists(filename_check)){
    
    output_figs <- paste0(output_init_figs, sample_name, "/")
    dir.create(output_figs, showWarnings = FALSE)
    
    print(file_list[1])
    print(path_file)
    
    
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
    
    # Spot diameter at full resolution
    spot_radius <- image_structure@scale.factors$spot_diameter_fullres / 2

    # Convert well matrix to spatstat ppp object
    spot_ppp <- spatstat.geom::as.ppp(cbind(X, Y), W = spatstat.geom::owin(range(X), range(Y)))
    
    
    ## ADD PATHOLOGY CELL COUNTS TO SEURAT OBJECT ----
    
    
    # Transform HALO nuclei segmentation output file to a spatstat ppp object
    nuclei_ppp <- halo_to_ppp_my(path_file, in_Anno=NULL, in_mpp=NULL, coreg_prefix="", in_anno_str="Tissue", 
                                 drop_marks=NULL, add_marks=NULL)
    
    # Filter the spatstat ppp object to remove background objects
    nuclei_ppp <- spatstat.geom::subset.ppp(nuclei_ppp, Classifier.Label=="Nuclei")

    # Compute nuclei within radius of each well
    spot_cellcounts <- spatstat.core::crosspaircounts(spot_ppp, nuclei_ppp, r = spot_radius)
    
    # Calculate the ratio of each class relative to the total of all non-reference classes
    data_seurat@meta.data[[col_prefix]] <- spot_cellcounts

    
    ## WRANGLE DATA ----
    
    
    # Wide tibble of signature scores + nCount_Spatial and seurat_clusters
    meta_loc_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")
    
    
    ## VISUALIZE CELL COUNTS: HISTOGRAM ----
    
    
    # Cell count stats to plot
    cell_count_min <- min(meta_loc_df[[col_prefix]])
    cell_count_mean <- round(mean(meta_loc_df[[col_prefix]]))
    cell_count_sd <- round(sd(meta_loc_df[[col_prefix]]))
    cell_count_max <- max(meta_loc_df[[col_prefix]])
    
    # A list of cell count stats
    stat_names <- c("NAME\tMIN\tMEAN-SD\tMEAN\tMEAN+SD\tMAX")
    stat_list <- c(cell_count_min, cell_count_mean-cell_count_sd, cell_count_mean, 
                   cell_count_mean+cell_count_sd, cell_count_max)

    # Write cell count stats to file
    filename <- paste0(output_figs, sample_name, "_cell_count_stat.txt")
    write(paste(stat_names, collapse="\t"), filename)
    write(paste(c(sample_name, stat_list), collapse="\t"), filename, append=TRUE)
    
    # Plot a cell count histogram
    filename <- paste0(output_figs, sample_name, "_cell_count_hist")
    p <- create_hist_plot_my(meta_loc_df, x_label=col_prefix, fill_label="orig.ident", 
                             intercept=stat_list, binwidth=1, filename=filename,
                             labels = c("Number of cells per spot", "Number of spots", sample_name))
    
    
    ## VISUALIZE DATA: SPATIAL PLOTS ----
    
    
    # Create a spatial plot of the tissue with no overlays
    p1 <- Seurat::SpatialDimPlot(data_seurat, pt.size.factor = 0) + 
      theme(legend.position = "none")
    
    # Cell density spatial plot
    p2 <- spatial_feature_plot_my(data_seurat, feature=col_prefix, min.cutoff="q1", max.cutoff="q99", 
                                  title=paste0(sample_name, "\nCell counts"))

    # UMI count spatial plot
    p3 <- spatial_feature_plot_my(data_seurat, feature="nCount_Spatial", min.cutoff="q1", max.cutoff="q99", 
                                  title=paste0(sample_name, "\nUMIs per spot"))

    # Feature count spatial plot
    p4 <- spatial_feature_plot_my(data_seurat, feature="nFeature_Spatial", min.cutoff="q1", max.cutoff="q99", 
                                  title=paste0(sample_name, "\nFeatures per spot"))
    
    # Write figure to file
    filename <- paste0(output_figs, sample_name, "_cell_count_tissue")
    write_plot2file_my(patchwork::wrap_plots(list(p1, p2, p3, p4), nrow=1), filename, num_col=4, num_row=1)

    
    ## VISUALIZE DATA: SCATTER ----
    
    
    # List of Seurat parameters of interest
    feature_list <- c("nCount_Spatial", "nFeature_Spatial")

    # Compare the number of cells with a few other parameters
    plot_list <- list()
    for(f in feature_list){
      plot_list[[f]] <- create_scatter_plot_my(meta_loc_df, x_label=col_prefix, y_label=f, 
                                               fill_label="orig.ident", filename=NULL, stroke=0.1,
                                               labels=c("Number of cells per spot", f, sample_name), do_fit=TRUE)
    }
    
    filename <- paste0(output_figs, sample_name, "_cell_count_scatter")
    write_plot2file_my(patchwork::wrap_plots(plot_list, nrow=1), filename, num_col=length(feature_list), num_row=1)

    
    ## OVERWRITE SEURAT FILE ----
    
    
    if(no_rds_output == FALSE){
      # Write the updated Seurat object
      filename <- paste0(output_path, sample_name, "_annotated_pathology.rds")
      saveRDS(data_seurat, file=filename)
    }
    
    gc()
    
  } else{
    warning("More than one file found!")
  }
}
