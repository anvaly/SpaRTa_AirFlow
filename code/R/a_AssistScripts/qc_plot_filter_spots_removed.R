# Author: Andy Kavran. Date 21-11-13

# Re-plots the image filtering figures after manually updating the pt.size.factor metadata fields



# Allow piping throughout the package
`%>%` <- magrittr::`%>%`

library(Seurat)

source("code/utils/utils_in_out.R")

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_image.R")
source("code/R/Utils/utils_10X_vis.R")


## PARAMETERS ----


group_name <- c("PDAC_759110QB_ROI1_s2","PDAC_759110QB_ROI3_s1","PDAC_759110QB_ROI3_s2","PDAC_759110QB_ROI4_s1","PDAC_759110QB_ROI4_s2","PDAC_1255880B_ROI1_s1","PDAC_1255880B_ROI2_s1","PDAC_1255880B_ROI2_s2","PDAC_1275301B_s3","PDAC_1255908B_s4","PDAC_1274423B_ROI1_v2_r1","PDAC_1274423B_ROI2_v2_r1","PDAC_1274423B_ROI2_v2_r2","PDAC_1275233B","PDAC_1275301B_s4","PDAC_1275301B_v1_r1","PDAC_1275301B_v1_r2","PDAC_1275301B_v2_r1","PDAC_1275301B_v2_r2","PDAC_1275309B_s1","PDAC_1275309B_s2","PDAC_E1265_ROI1","PDAC_E1265_ROI2","PDAC_E1265_ROI4","PDAC_E2547_ROI1","PDAC_E2547_ROI2","PDAC_E2547_ROI3","PDAC_E2547_ROI4","PDAC_E5058_ROI1","PDAC_E5058_ROI2","PDAC_E5058_ROI3","PDAC_E5058_ROI4","PDAC_E27038_ROI4","PDAC_ILS52188PT1","PDAC_Pt2_ROI2","PDAC_Pt2_ROI3","PDAC_Pt3_ROI1_s2","PDAC_Pt3_ROI3_s1","PDAC_Pt3_ROI3_s2","PDAC_Pt3_ROI4_s2","PDAC_Pt4_ROI4","PDAC_TXG_v1_r1","PDAC_TXG_v1_r2","PDAC_TXG_v2_r1","PDAC_TXG_v2_r2")
# group_name <- "PDAC_1275301B_v1_r2"
remove_promiscuous_genes <- FALSE
apply_feature_threshold <- FALSE

# Input location containing 10X folders
input_path <- "XXXX"

# Output location for Seurat objects
output_path <- "XXXX"
output_folders <- create_output_subfolders_my(output_path, c("Seurat_object", "Seurat_diet",
                                                             "Figures_filter_spots_removed", "Figures_check_gene"))
# Input file with sample meta data
meta_file <- paste0(input_path, "full_pdac_meta_data.txt")

## INGEST, WRANGLE, OUTPUT ----
# Read sample file, filter for samples that pass QC
meta_df <- readr::read_delim(file=meta_file, delim="\t") %>% dplyr::filter(Quality_Status == TRUE)

if(!is.null(group_name)){
  meta_df <- meta_df %>%
    dplyr::filter(Sample_Name %in% group_name)
}


for(i in c(1:nrow(meta_df))){
#for(i in 1){
  # Load the 10X Spatial folder
  filename_in <- meta_df[[i, "FullPath"]]
  sample_name <- meta_df[[i, "Sample_Name"]]
  split_area <- meta_df[[i, "Split"]]
  color_trim_sd <- meta_df[[i, "Color_trim_sd"]]
  
  filename_out3 <- paste0(output_folders$Figures_filter_spots_removed, sample_name, "_spots_removed")
  
  print(filename_in)
  
  
  ## Load data ----
  
  # Read a 10X Spatial folder
  data_seurat_init <- read_10X_spatial_folder_my(filename_in)
  # Read a JSON with 10X defined image properties
  data_seurat_init <- update_seurat_image_properties_from_json_my(data_seurat_init, filename_in)
  # Read a 10X defined metrics fine
  qc_df <- read_10X_qc_metrics_my(filename_in)
  
  #Temp Fix recast  data_seurat_raw@images$slice1@coordinates as integers
  these_cols <-colnames(data_seurat_init@images$slice1@coordinates)
  data_seurat_init@images$slice1@coordinates[these_cols] <- lapply(data_seurat_init@images$slice1@coordinates[these_cols], as.integer)
  
  ## Remove promiscuous genes ----
  
  # Remove genes that were identified by data trends QC as having expression outside of the tissue
  if(remove_promiscuous_genes == TRUE){
    promiscuous_genes <- readr::read_delim(paste0(filename_in, "promiscuious_genes.txt"), delim="\t", col_names="Symbol")
    data_seurat_init <- subset(data_seurat_init, features=setdiff(rownames(data_seurat_init), promiscuous_genes$Symbol))
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
    data_seurat_filter1 <- subset(data_seurat_init, cells=barcode_keep)
  } else{
    data_seurat_filter1 <- data_seurat_init
  }
  
  
  ## Apply UMI threshold ----
  
  if(apply_feature_threshold == TRUE){
    # Find UMI mean and stdev
    val_mean <- mean(data_seurat_filter1@meta.data$nFeature_Spatial)
    val_std <- sd(data_seurat_filter1@meta.data$nFeature_Spatial)
    val_threshold <- val_mean - val_std * 3
    
    # Select spots that have at sufficient amount of UMIs
    barcodes_select <- rownames(data_seurat_filter1@meta.data)[which(data_seurat_filter1@meta.data$nFeature_Spatial >= val_threshold)]
    data_seurat_filter2 <- subset(data_seurat_filter1, cells = barcodes_select)
  } else{
    # Select spots that have at least 100 Features
    barcodes_select <- rownames(data_seurat_filter1@meta.data)[which(data_seurat_filter1@meta.data$nFeature_Spatial >= 100)]
    data_seurat_filter2 <- subset(data_seurat_filter1, cells = barcodes_select)
  }
  
  
  ## Contiguity filter ----
  
  # Remove any floating tissue spots
  data_seurat_filter3 <- tissue_contiguity_filter_my(data_seurat_filter2, spot_num=200)
  
  
  ## Update meta data ----
  
  # Add user-defined meta-data to a Seurat object misc slot tagged with "user." prefix
  for(v in colnames(meta_df)){
    data_seurat_filter3 <- add_misc_to_seurat_object_my(data_seurat_filter3, field=paste0("user.", v), dict_value=meta_df[[i, v]])
  }
  
  # Add 10X default QC metrics to a Seurat object tagged with "10X." prefix
  if(!is.null(qc_df)){
    for(col in 1:ncol(qc_df)){
      data_seurat_filter3 <- add_misc_to_seurat_object_my(data_seurat_filter3, field=paste0("10X.", colnames(qc_df)[col]), dict_value=qc_df[[1, col]])
    }
  }
  
  # Select a well expressed gene
  #gene_list <- seurat_select_abundant_genes_my(data_seurat_filter3, sct_threshold=3, spot_threshold=50)
  gene_list <- seurat_select_abundant_genes_my(data_seurat_filter3, sct_threshold=3, spot_threshold=50, assay = "Spatial")
  
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
                                name=paste0(sample_name, "\n", gene_list[1], "\nSpot number = ", length(colnames(data_seurat_filter3))))
  
  # Write figures to file
  write_plot2file_my(patchwork::wrap_plots(list(p1, p_density, p2, p3, p4, p5), nrow=1), filename_out3, num_row=1, num_col=6)
}