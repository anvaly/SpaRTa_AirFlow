# Authors: Anna Lyubetskaya, Andrew Kavran, Hannah Pliner, Yulong Bai. Date: 20-11-27

# Important:
# Keep Sample_Name values unique at least within your major cohorts - these are used extensively throughout
# This script names the output objects based on the meta data file Sample_Name column.
# Only use underscores or fullstops in the Sample_Name nomenclature! (E.g., no dashes) Sample_Name must start with letters!


# Create a normalized, filtered Seurat object from Spatial Transcriptomics 10X folder data

# For a list of samples from a meta data file, 
# 1. load raw and processed SpaceRanger data (upper case all gene symbols)
# 2. plot UMI/feature distributions outside of the tissue to demonstrate diffusion
# 3. remove unspecific probe genes identified by data trend QCs
# 4. evaluate the precision of spot centers
# 5. apply SpotClean
# 6. apply two filters using ST images
# --- color value filter - remove spots located in too-light-to-be-tissue areas (input Color_trim_sd parameter);
# --- contiguity filter - remove floater spots not connected to the main tissue slice;
# 7. augment meta data
# --- data from the input meta data file, with prefix "user." in the misc slot
# --- data from the 10X metrics file, with prefix "10X." in the misc slot
# --- spot diameter in original resolution from 10X image JSON, in data_seurat@images$slice1@scale.factors
# 8. calculate ribo / mito content if applicable
# 9. calculate additional QC metrics
# 10. SCT transform using vst.flavor = "v2"
# 11. cluster and run DEA for a selected resolution
# 12. save to RDS as a normal and diet Seurat objects

# Input meta data file contain the following columns:
# - Required:
# --- Sample_Name and Pipeline.FullPath: Sample_Name is the RDS object output name, FullPath is the location of the 10X spatial report
# - Optional:
# --- Pipeline.Quality_Status - boolean if the data passes first pass QC and should be analyzed by the script
# --- Pipeline.pt.size.factor - size of the spot for Seurat spatial visualizations
# --- Pipeline.Split - to indicate whether the dataset should be split into multiple samples according to the number of tissue slices identified
# --- Pipeline.Color_trim_sd - number of SD to cleanup spots by color
# --- Pipeline.Unspecific_Genes_File - a list of genes to remove from the object due to the concern over lack of specificity of the probe
# --- Pipeline.Feature_threshold - a feature threashold to remove bad spots
# --- Pipeline.Contiguity_spot_num - a number of spots if they create a floater island away from the main tissue


## ENVIRONMENT ----

print('--- Running 10X_to_Seurat.R')

# Allow piping throughout the package
`%>%` <- magrittr::`%>%`

# Update and install packages
if(packageVersion("Seurat") < 4.1){
  install.packages("Seurat", repos=c('https://pm.rdcloud.bms.com/prod-cran/latest','https://pm.rdcloud.bms.com/prod-cran/2022-11-02', 'https://pm.rdcloud.bms.com/bms-cg-biogit-bran/latest'))
}
print('--- Installed Seurat')

if(packageVersion("glmGamPoi") < 1.6){
  BiocManager::install("glmGamPoi")
}

if(!"SpotClean" %in% installed.packages()){
  devtools::install_github("zijianni/SpotClean", build_manual = TRUE, build_vignettes = TRUE)
}

# Check we have the right python packages and if not install
if (!reticulate::py_module_available('cv2')) {
  system('pip install opencv-python')
}

print('--- Installed opencv-python')

# library(Seurat)
# library("optparse")

source("code/utils/utils_in_out.R")

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_image.R")
source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_dea.R")

source("code/R/a_Wrangle/10X_to_Seurat_ST_utils.R")
source("code/R/b_Cluster_Biomarkers/Seurat_cluster_utils.R") # check

print('--- Sourced all R scripts')

## FIXED PARAMETERS ----


# A list of genes to spot check SpotClean decontamination
# When set to a number, top n most changed genes before and after spotclean will be used. 
# Example: spotclean_gene_list = 5
spotclean_gene_list <- 5

# Set defaults for optional parameters
split_area <- FALSE
color_trim_sd <- 1
feature_threshold <- 100
contiguity_spot_num <- 5
unspecific_genes_file <- NULL

# Setup parameters
params <- cluster_params_my()
params[["test_use"]] <- "wilcox"  # wilcox; "LR", "negbinom", "poisson", or "MAST"
params[["assay"]] <- "SCT"  # Spatial or SCT


## PATHS ----


# Set default values for the two parameters
default_meta_file <- "XXXX"
default_output_path <- "XXXX"

# Parse the arguments passed in from the cli (sets to defaults if none provided)
option_list <- list(
  optparse::make_option(c("-m", "--meta_file"), action="store", 
                        default=default_meta_file, help="Parse the path to the meta file/factor sheet."),
  optparse::make_option(c("-o", "--output_path"), action="store", 
                        default=default_output_path, help="Parse the path to the output location."))

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list)) 

# Meta data file listing and annotating all samples for processing
if (opt$meta_file == 'default'){
  meta_file <- default_meta_file
} else {
  meta_file <- opt$meta_file    
}


# Check if provided meta file exists
if (!file.exists(meta_file)){
  stop(paste0("ERROR: The following file does not exists: ", meta_file))
}

# Output location for Seurat objects
output_path <- opt$output_path

if (!endsWith(output_path, '/')){
  output_path <- paste0(output_path, '/')
}

# if (dir.exists(output_path)){
#   unlink(output_path, recursive=TRUE)
# }

# Create output subfolders
dir.create(output_path, recursive=TRUE)
output_folders <- create_output_subfolders_my(output_path, c("2a_Figures_spatial_coverage",
                                                             "2b_Figures_filter_spots_removed", 
                                                             "2c_Figures_check_gene",
                                                             "2d_Figures_spotclean", "2d_Matrices_spotclean_topgene",
                                                             "2e_Figures_spotdetect", 
                                                             "2f_Clustering", 
                                                             "2g_Seurat_object", "2g_Seurat_diet"))

# Write meta data file to specified output location 
file.copy(meta_file, output_path)


## INGEST, WRANGLE, OUTPUT ----


# Read sample file, filter for samples that pass QC
meta_df <- readr::read_delim(file=meta_file, delim="\t")

# Check that sample names are unique
if(length(unique(meta_df$Sample_Name)) != length(meta_df$Sample_Name)){
  stop("ERROR: Sample_Name column doesn't contain unique values. Pipeline terminated")
}

# Check that sample names are unique
if(stringr::str_detect(paste(meta_df$Sample_Name, collapse=""), "[^a-zA-Z0-9_\\.]+")){
  stop("ERROR: Sample_Name column should only contain letters, nubers, underscores, or periods. Pipeline terminated")
}

# Check if the user defined Pipeline.Quality_Status and if yes, remove samples that don't pass
# Without this parameter, all samples in the meta file will be processed
if("Pipeline.Quality_Status" %in% colnames(meta_df)){
  meta_df <- meta_df %>% 
    dplyr::filter(Pipeline.Quality_Status != FALSE)
}


for(i in 1:nrow(meta_df)){
  
  print(i)
  
  
  ## Set parameters ----
  
  
  # Load the 10X Spatial folder and sample name - required parameters
  sample_name <- meta_df[[i, "Sample_Name"]]
  filename_in <- meta_df[[i, "Pipeline.FullPath"]]
  
  # Load optional parameters
  if("Pipeline.Split" %in% colnames(meta_df)){
    split_area <- meta_df[[i, "Pipeline.Split"]]
  }
  
  if("Pipeline.Color_trim_sd" %in% colnames(meta_df)){
    color_trim_sd <- meta_df[[i, "Pipeline.Color_trim_sd"]]
  }
  
  if("Pipeline.Unspecific_Genes_File" %in% colnames(meta_df)){
    unspecific_genes_file <- meta_df[[i, "Pipeline.Unspecific_Genes_File"]]
  }
  
  if("Pipeline.Feature_threshold" %in% colnames(meta_df)){
    feature_threshold <- meta_df[[i, "Pipeline.Feature_threshold"]]
  }
  
  if("Pipeline.Contiguity_spot_num" %in% colnames(meta_df)){
    contiguity_spot_num <- meta_df[[i, "Pipeline.Contiguity_spot_num"]]
  }
  
  
  filename_out3 <- paste0(output_folders[["2b_Figures_filter_spots_removed"]], sample_name, "_spots_removed")
  filename_out4 <- paste0(output_folders[["2d_Figures_spotclean"]], sample_name, "_spotclean_qc")
  
  if(!file.exists(paste0(filename_out3, ".png"))){
    
    print(filename_in)
    
    
    ## Load data ----
    
    
    # Ingest raw (unfiltered) data
    data_seurat_raw <- read_10X_spatial_folder_my(filename_in, filename="raw_feature_bc_matrix.h5")
    
    # Ingest filtered data
    data_seurat_filt <- read_10X_spatial_folder_my(filename_in)
    
    ## For data generated by CytAssist, coordinates features in image slot will be converted to character by Read10X_Image
    ## Line below makes sure coordinates are integer.
    data_seurat_filt@images$slice1@coordinates[colnames(data_seurat_filt@images$slice1@coordinates)] <- lapply(data_seurat_filt@images$slice1@coordinates[colnames(data_seurat_filt@images$slice1@coordinates)], as.integer)
    data_seurat_raw@images$slice1@coordinates[colnames(data_seurat_raw@images$slice1@coordinates)] <- lapply(data_seurat_raw@images$slice1@coordinates[colnames(data_seurat_raw@images$slice1@coordinates)], as.integer)
    
    
    ## Plot spatial data distributions ----
    
    
    # UMI count spatial plot - all spots
    p1 <- spatial_feature_plot_my(data_seurat_raw, feature="nCount_Spatial",
                                  title=paste0(sample_name, "\nUMIs per spot"), crop=FALSE, color = Turbo())
    
    # Feature count spatial plot - all spots
    p2 <- spatial_feature_plot_my(data_seurat_raw, feature="nFeature_Spatial",
                                  title=paste0(sample_name, "\nFeatures per spot"), crop=FALSE, color = Turbo())
    
    # UMI count spatial plot - tissue spots
    p3 <- spatial_feature_plot_my(data_seurat_filt, feature="nCount_Spatial",
                                  title=paste0(sample_name, "\nUMIs per spot"), crop=FALSE, color = Turbo())
    
    # Feature count spatial plot - tissue spots
    p4 <- spatial_feature_plot_my(data_seurat_filt, feature="nFeature_Spatial",
                                  title=paste0(sample_name, "\nFeatures per spot"), crop=FALSE, color = Turbo())
    
    # Create a spatial plot of the tissue with no overlays
    p5 <- Seurat::SpatialDimPlot(data_seurat_raw, pt.size.factor = 0, crop=FALSE) + 
      theme(legend.position = "none", plot.title = element_text(size=8)) +
      ggplot2::labs(x="", y="", title=paste0("Spot number = ", length(colnames(data_seurat_raw))))
    
    # Create a spatial plot of the tissue with no overlays
    p6 <- Seurat::SpatialDimPlot(data_seurat_filt, pt.size.factor = 0, crop=FALSE) + 
      theme(legend.position = "none", plot.title = element_text(size=8)) +
      ggplot2::labs(x="", y="", title=paste0("Spot number = ", length(colnames(data_seurat_filt))))
    
    # Write figure to file
    filename <- paste0(output_folders[["2a_Figures_spatial_coverage"]], "spatial_coverage_", sample_name)
    write_plot2file_my(patchwork::wrap_plots(list(p5, p1, p2, p6, p3, p4), nrow=2, ncol=3), filename, num_row=1, num_col=2)
    
    
    ## Run spotdetect script to measure spot placement error ---- 
    
    
    system(paste0('code/utils/spotdetect.py ', filename_in, ' -o ', output_folders[["2e_Figures_spotdetect"]], ' -s ', sample_name))
    
    # Calculate error metrics
    spotdetect_df <- gen_spotdetect_df_my(filename_in, 
                                          paste0(output_folders[["2e_Figures_spotdetect"]], "/", sample_name, "_spots_detected.csv"), 
                                          max_error = 50) 
    
    # Generate summary fig
    generate_spotdetect_summary_fig_my(spotdetect_df, sample_name, 
                                       paste0(output_folders[["2e_Figures_spotdetect"]], "/", sample_name,  "_spotdetect_error.png") )
    
    
    ## SpotClean for diffusion ----
    
    
    # Attempt to deploy SpotClean, it errors out if the data is too bad of a quality so catch it
    data_seurat_init <- tryCatch({
      
      if("SpotClean" %in% rownames(installed.packages())){
        data_seurat_init <- apply_spotclean_my(filename_in, filename_out4)
      } else{
        data_seurat_init <- data_seurat_filt
      }
      
    },
    error = function(err) {
      print(paste("MY ERROR: ", err))
    })
    
    if(is.null(dim(data_seurat_init))){
      next
    }
    
    
    ## Read image information ----
    
    
    # Read a JSON with 10X defined image properties
    data_seurat_init <- update_seurat_image_properties_from_json_my(data_seurat_init, filename_in)
    # Read a 10X defined metrics fine
    qc_df <- read_10X_qc_metrics_my(filename_in)
    
    # Temp Fix: recast @images$slice1@coordinates as integers
    these_cols <- colnames(data_seurat_init@images$slice1@coordinates)
    data_seurat_init@images$slice1@coordinates[these_cols] <- lapply(data_seurat_init@images$slice1@coordinates[these_cols], as.integer)
    
    
    ## Remove promiscuous genes ----
    
    
    # Remove / correct for genes that were identified by data trends QC as having expression outside of the tissue
    if(!is.null(unspecific_genes_file)){
      if(!is.na(unspecific_genes_file)){
        unspecific_genes <- readr::read_delim(unspecific_genes_file, delim="\t")
        data_seurat_init <- subset(data_seurat_init, features=setdiff(rownames(data_seurat_init), unspecific_genes$Symbol))
      }
    }
    
    
    ## Tissue color filter ----
    
    
    # Calculate the mean color value of each spot in Seurat image matrix
    spot_color_df <- image_mean_color_value_my(data_seurat_init, num_sd=color_trim_sd)
    
    # Find barcodes corresponding to spots with unusually light color
    barcode_remove <- spot_color_df %>%
      dplyr::filter(Exclude == TRUE) %>%
      dplyr::pull(Barcode)
    
    # Find barcodes corresponding to spots under tissue
    barcode_keep <- spot_color_df %>%
      dplyr::filter(Exclude == FALSE) %>%
      dplyr::pull(Barcode)
    
    # Create a histogram of mean color values for each spot
    p_density <- create_hist_plot_my(spot_color_df, x_label="Value_mean", fill_label="Exclude", 
                                     binwidth=max(spot_color_df$Value_mean)/100, filename=NULL)
    
    # Subset the Seurat object
    if(length(barcode_remove) > 0){
      
      # Error the pipeline out if no good spots are found
      if(length(barcode_keep) < 20){
        cat("ERROR: Less than 20 spots found that pass the color intensity condition of", color_trim_sd, "standard deviations\n")
        break()
      }
      
      data_seurat_filter1 <- subset(data_seurat_init, cells=barcode_keep)
    } else{
      data_seurat_filter1 <- data_seurat_init
    }
    
    
    ## Apply UMI threshold ----
    
    
    # Select spots that have at least 300 genes
    barcodes_select <- rownames(data_seurat_filter1@meta.data)[which(data_seurat_filter1@meta.data$nFeature_Spatial >= feature_threshold)]
    
    # Error the pipeline out if no good spots are found
    if(length(barcodes_select) < 20){
      cat("ERROR: Less than 20 spots found that satisfy the UMI threshold (nFeature_Spatial) of", feature_threshold, "\n")
      break()
    }
    
    data_seurat_filter2 <- subset(data_seurat_filter1, cells = barcodes_select)
    
    
    ## Contiguity filter ----
    
    
    # Remove any floating tissue spots (if environment contains a working spdep library)
    if("spdep" %in% rownames(installed.packages()) && contiguity_spot_num > 0){      
      data_seurat_filter3 <- tissue_contiguity_filter_my(data_seurat_filter2, spot_num=contiguity_spot_num)
    } else{
      data_seurat_filter3 <- data_seurat_filter2
    }
    
    
    ## Update meta data ----
    
    
    # Add user-defined meta-data to a Seurat object misc slot tagged with "user." prefix
    for(v in colnames(meta_df)){
      data_seurat_filter3 <- add_misc_to_seurat_object_my(data_seurat_filter3, field=paste0("user.", v), dict_value=meta_df[[i, v]])
    }
    # Add the sample name to the meta.data slot for downstream compatibility
    data_seurat_filter3@meta.data$user.Sample_Name <- sample_name
    
    # Add 10X default QC metrics to a Seurat object tagged with "10X." prefix
    if(!is.null(qc_df)){
      for(col in 1:ncol(qc_df)){
        data_seurat_filter3 <- add_misc_to_seurat_object_my(data_seurat_filter3, field=paste0("10X.", colnames(qc_df)[col]), dict_value=qc_df[[1, col]])
      }
    }
    
    # Select a well expressed gene
    gene_list <- seurat_select_abundant_genes_my(data_seurat_filter3, sct_threshold=3, spot_threshold=50, assay="Spatial")
    if(length(gene_list) == 0){
      gene_list <- seurat_select_abundant_genes_my(data_seurat_filter3, sct_threshold=1, spot_threshold=1, assay="Spatial")
    }
    
    
    ## Visualize image filtering steps ----
    
    
    # Create a spatial plot of the tissue with no overlays
    p1 <- Seurat::SpatialDimPlot(data_seurat_init, pt.size.factor = 0) + 
      theme(legend.position = "none") +
      ggplot2::labs(x="", y="", title=paste0("Spot number = ", length(colnames(data_seurat_init))))
    
    # Create a spatial plot highlighting unusually light spots
    p2 <- spatial_dim_plot_my(data_seurat_init, setdiff(colnames(data_seurat_init), colnames(data_seurat_filter1))) + 
      theme(legend.position = "none") +
      ggplot2::labs(x="", y="", title=paste0("Spots removed by color = ", length(setdiff(colnames(data_seurat_init), colnames(data_seurat_filter1)))))
    
    # Visualize spots removed by the UMI threshold
    p3 <- spatial_dim_plot_my(data_seurat_filter1, setdiff(colnames(data_seurat_filter1), colnames(data_seurat_filter2))) + 
      theme(legend.position = "none") +
      ggplot2::labs(x="", y="", title=paste0("Spots removed by feature number = ", length(setdiff(colnames(data_seurat_filter1), colnames(data_seurat_filter2)))))
    
    # Visualize spots removed by the contiguity filter
    p4 <- spatial_dim_plot_my(data_seurat_filter2, setdiff(colnames(data_seurat_filter2), colnames(data_seurat_filter3))) + 
      theme(legend.position = "none") +
      ggplot2::labs(x="", y="", title=paste0("Spots removed by contiguity = ", length(setdiff(colnames(data_seurat_filter2), colnames(data_seurat_filter3)))))
    
    # Visualize the data before and after the contiguity filter using a well expressed gene
    p5 <- spatial_feature_plot_my(data_seurat_filter3, gene_list[1], 
                                  title=paste0(sample_name, "\n", gene_list[1], "\nSpot number = ", length(colnames(data_seurat_filter3))))
    
    # Write figures to file
    write_plot2file_my(patchwork::wrap_plots(list(p1, p_density, p2, p3, p4, p5), nrow=1), filename_out3, num_row=1, num_col=6)
    
    
    ## Update barcode and image names ----
    
    
    # Change the image name from "slice1" to the sample name - important for integrative scripts
    names(data_seurat_filter3@images) <- sample_name
    
    # Add the sample name to barcode names - important for integrative scripts
    data_seurat_filter3 <- SeuratObject::RenameCells(data_seurat_filter3,
                                                     new.names = paste0(sample_name, ":", colnames(data_seurat_filter3)))
    
    
    ## Split a Visium area into tissue slices ----
    
    
    if(split_area == "TRUE"){
      # Identify tissue slices
      spot_clusters <- seurat_image_to_clusters_my(data_seurat_filter3, spot_distance=2)
      
      seurat_list <- list()
      for(cluster in unique(spot_clusters$membership)){
        # Select barcodes corresponding to the cluster
        barcode_list <- names(spot_clusters$membership)[which(spot_clusters$membership == cluster)]
        
        # Subset the Seurat object
        data_seurat <- subset(data_seurat_filter3, cells=barcode_list)
        
        # Add the partial Seurat object to a list
        seurat_list[[as.character(cluster)]] <- data_seurat
      }
    } else{
      seurat_list <- list(NONE = data_seurat_filter3)
    }
    
    
    for(name in names(seurat_list)){
      # Select Seurat dataset
      data_seurat <- seurat_list[[name]]
      
      
      ## Calculate ribosomal and mitochondrial content ----
      
      
      # Calculate mitochondrial and ribosomal content
      data_seurat <- calculate_mt_ribo_my(data_seurat)
      
      
      ## Cycle break if data gone ----
      
      
      if(min(dim(data_seurat)) < 100){
        next
      }
      
      
      ## SCTransform ----
      
      
      # Variables to regress out
      vars_to_regress <- c("mito_percent", "ribo_percent")    
      
      # Perform SCTransform normalization
      data_seurat <- Seurat::SCTransform(data_seurat, assay="Spatial", vars.to.regress=vars_to_regress,
                                         residual.features=NULL, variable.features.n=5000, 
                                         return.only.var.genes=FALSE, verbose=FALSE, vst.flavor="v2")
      
      
      ## Perform clustering and DEA analysis ----
      
      
      # Perform PCA, UMAP, FindNeighbors, FindClusters, and FindMarkers
      data_seurat <- cluster_analysis_my(data_seurat, params, sample_name, output_folders[["2f_Clustering"]])
      
      # Identify the best clustering resolution as the smallest resolution corresponding to the mean number of clusters
      if(!"user.Pipeline.Clustering" %in% names(data_seurat@misc) || 
         is.na(nchar(data_seurat@misc$user.Pipeline.Clustering)) ||
         nchar(data_seurat@misc$user.Pipeline.Clustering) < 13){
        
        clust_num <- sapply(params$resolution_list, function(x) length(unique(data_seurat@meta.data[[paste0("SCT_snn_res.", x)]])))
        data_seurat@misc$user.Pipeline.Clustering <- paste0("SCT_snn_res.", params$resolution_list[which.min(abs(clust_num - floor(mean(clust_num))))])
        
      }
      
      # Find markers for the selected resolution
      data_seurat <- cluster_analysis_markers_my(data_seurat, params, sample_name, output_folders[["2f_Clustering"]])
      
      
      ## Calculate mean QC metrics on final objects ----
      
      
      # List of meta data fields to parse
      field_list <- c("nCount_Spatial", "nFeature_Spatial", "nCount_SCT", "nFeature_SCT", vars_to_regress)
      
      # Add mean QC metrics to the misc slot
      for(f in field_list){
        val <- mean(data_seurat@meta.data[[f]])
        data_seurat <- add_misc_to_seurat_object_my(data_seurat, field=paste0("qc.mean.", f), dict_value=val)
      }
      
      
      ## Write data ----
      
      
      if(name == "NONE"){
        filename_out1 <- paste0(output_folders[["2g_Seurat_object"]], sample_name, "_all.rds")
        filename_out2 <- paste0(output_folders[["2g_Seurat_diet"]], sample_name, "_diet.rds")
        
      } else{
        # Update Sample_Name field within the Seurat object if slice suffix was added
        data_seurat@misc$user.Sample_Name <- paste0(sample_name, "_slice", name)
        
        filename_out1 <- paste0(output_folders[["2g_Seurat_object"]], sample_name, "_slice", name, "_all.rds")
        filename_out2 <- paste0(output_folders[["2g_Seurat_diet"]], sample_name, "_slice", name, "_diet.rds")
      }
      
      # Save the full dataset to RDS
      saveRDS(data_seurat, file = filename_out1, compress = FALSE)
      
      # Strip away all unnecessary parts
      data_seurat_diet <- Seurat::DietSeurat(data_seurat, assays = "SCT")
      
      # Dummy up pathology values for the BMS Portal App
      data_seurat_diet@meta.data[["Pathology.Group"]] <- "Tissue"
      
      # Save the diet dataset to RDS
      saveRDS(data_seurat_diet, file = gsub("NONE", "", gsub(".rds", paste0(name, ".rds"), filename_out2)), compress = FALSE)
    }
    
    
    ## Visualize a random gene ----
    
    
    # Visualize a random gene using the last version of the data
    p_list <- batch_spatial_feature_plot_my(seurat_list, gene_list[1])
    
    # Write figures to file
    filename_out4 <- paste0(output_folders[["2c_Figures_check_gene"]], sample_name, "_gene_check")
    write_plot2file_my(patchwork::wrap_plots(p_list, nrow=1), filename_out4, num_row=1, num_col=length(names(seurat_list)))

  }
  
  gc()
}
