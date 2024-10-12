# Author: Anna Lyubetskaya. Date: 21-04-22
# Select specific spots based on meta data and create H&E "swatches"

# Important concern: The swatches extracted from different images are of a different size - matches the original microscope seetting of the image


## ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_image.R")
source("code/R/Utils/utils_10X_signatures.R")
source("code/R/Utils/utils_pathology_spots.R")

source("code/R/b_Cluster_Biomarkers/Seurat_cluster_utils.R")


## PARAMETERS ----


# Path to processed Seurat data
sample_name <- "PDAC108_path14_5K_harmony"
sample_exclude <- NULL

# Name for the signature group being plotted
sig_select <- "PDAC.P19.Bcell"  # e.g., NULL or "PDAC.collisson.classical", "PDAC.moffitt.activatedstroma"
sig_empirical <- FALSE  # If want to swap scores for bootstrapped p-values
sig_pval_threshold1 <- 0  # score, e.g. 0.5
sig_pval_threshold2 <- 0  # score, e.g. 0

# Name of the pathology field to filter by groups
pathology_select <- "Pathology.TLSMature.percent"  # e.g., "Tumor"
pathology_threshold1 <- 50  # upper, 75
pathology_threshold2 <- 0  # lower, 25

# SCT and spot number thresholds for gene filtering
sct_threshold <- 0.5
spot_threshold <- 5

# Assay and slot to use to calculate signature score
assay <- "SCT"
slot <- "data"

# User clustering resolution
resolution_user <- "integrated_snn_res.0.1"

# Group swatches by
swatch_group <- "user.Sample_Name"  # e.g., user.Sample_Name

# Number of spots to return
topn_param <- 40

# Don't write the final RDS object
no_rds_output <- TRUE

# Name of the classifier to use in the Seurat object meta data
classifier_name <- "Swatches_TLS"

# Subset to specific clusters for modeling if needed
cluster_select <- NULL  # c(1, 3, 2, 9)

# Remove spots with a lot of white space or blur
remove_bad_spots <- FALSE


## PATHS ----


# Path to a signature file
sig_path <- "data/import/Signatures/signatures_pdac_20230720.txt"

# Location of pre-processed data
input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, classifier_name, "_", sample_name, "_", sig_select, "_", pathology_select, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Seurat data
data_seurat <- readRDS(input_path)

# Remove samples if necessary
if(!is.null(sample_exclude)){
  barcode_list <- rownames(data_seurat@meta.data[which(!data_seurat@meta.data[["user.Sample_Name"]] %in% sample_exclude),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# Remove clusters if necessary
if(!is.null(cluster_select)){
  barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data[[resolution_user]] %in% cluster_select),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# Remove bad spots by pathology if necessary
if(remove_bad_spots == TRUE){
  barcode_list <- rownames(data_seurat@meta.data)[which(data_seurat@meta.data$Pathology.White_Space.percent <= 10 & data_seurat@meta.data$Pathology.Blur_Tissue.percent <= 10)]
  if(length(barcode_list) > 0){
    data_seurat <- subset(data_seurat, cells=barcode_list)
  }
}


## CALCULATE TARGET SIGNATURE SCORES ----


## Calculate signature scores and add them to Seurat meta data

# Find genes abundant in this sample
gene_list <- seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold, assay=assay, slot=slot)

# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=c(sig_select), sig_length_min=1, sig_length_max=1000, ratio_threshold=0)
cat("Target signature length =", length(signature_list[[1]]))

# Add signature scores to a seurat object
data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.", assay=assay)

# Find column names in the Seurat object - Seurat can rename original signatures
sig_names <- colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))])

# Column name for the signature empirical p-value if available
if(sig_empirical == FALSE){
  sig_name_pval <- sig_names
}else {
  sig_name_pval <- paste0(sig_names, ".EPvalue")  
}


## CREATE A RANDOM SIGNATURE SCORE BACKGROUND ----


if(sig_empirical == TRUE){
  # Write/read random gene expression scores to file
  filename <- paste0(output_path, "/sig_random_", sample_name, "_", sig_names, "_", pathology_select, ".txt")
  
  # Create a random signature background distribution and score the target signature against it
  data_seurat <- signature_empirical_pvalue_my(data_seurat, signature_list[[1]], output_path, sample_name, sig_names[1], sig_random_filename=filename,
                                               n_simulations=1000, sct_threshold=1, spot_threshold=10, assay="SCT", slot="data", output_col_name=sig_name_pval)
}


## WRANGLE META DATA ----


# Extract coordinate meta data
meta_init_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

# Add user.Sample_Name to meta data if absent: happens in stand-alone samples
if(!"user.Sample_Name" %in% colnames(meta_init_df)){
  meta_init_df["user.Sample_Name"] = data_seurat@misc$user.Sample_Name
}

# Select relevant meta data
meta_df <- meta_init_df  %>%
  dplyr::select(dplyr::all_of(unique(c("Coordinate", "user.Sample_Name", sig_name_pval, pathology_select, resolution_user, swatch_group))))

# Anonymize two main parameters by which spots are selected
colnames(meta_df) <- gsub(sig_name_pval, "Signature", colnames(meta_df))
colnames(meta_df) <- gsub(pathology_select, "Pathology", colnames(meta_df))


## VISUALIZE DISTRIBUTIONS ----


# Find mean signature and pathology score by cluster
clust_df <- meta_df %>%
  dplyr::group_by(user.Sample_Name, !!rlang::sym(swatch_group)) %>%
  dplyr::summarise(SignatureMean = mean(Signature),
                   PathologyMean = mean(Pathology))

# Define custom colors
cols <- define_cols_my(n=length(unique(clust_df[[swatch_group]])))
names(cols) <- sort(unique(clust_df[[swatch_group]]))

# Plot cluster mean signature score v cluster mean pathology score  
filename <- paste0(output_path, "scatter_clust_", sample_name, "_", sig_select, "_", pathology_select)
p <- create_scatter_plot_my(clust_df, x_label="PathologyMean", y_label="SignatureMean", 
                            fill_label=resolution_user, facet_var=c("user.Sample_Name", "free_y"), 
                            filename=filename, size=5, labels=c(pathology_select, sig_select, sample_name), cols=cols)

# Plot spot signature score v spot pathology score
filename <- paste0(output_path, "scatter_spot_", sample_name, "_", sig_select, "_", pathology_select)
p <- create_scatter_plot_my(meta_df, x_label="Pathology", y_label="Signature", 
                            fill_label=resolution_user, facet_var=c("user.Sample_Name", "fixed"), 
                            filename=filename, size=0.5, labels=c(pathology_select, sig_select, sample_name), stroke=0, cols=cols)


## SELECT SPOTS ----


extract_topn_my <- function(meta_df, topn_param=5, category=""){
  # Extract top N hits in a tibble based on two parameters
  
  meta_df <- meta_df %>%
    dplyr::slice_head(n = topn_param) %>%
    dplyr::mutate(Category = category)
  
  cat(category, nrow(spot_lists[[category]]), "\n")
  
  return(meta_df)
  
}


# Swatch group names
sample_list <- sort(unique(meta_df[[swatch_group]]))


spot_lists <- list()

# Threshold each sample individually
for(s in sample_list){
  
  meta_loc_df <- meta_df %>%
    dplyr::filter(!!rlang::sym(swatch_group) == s)
  
  # Pathology high, signature high
  category <- "Hi-Hi"

  spot_lists[[paste(category, s)]] <- meta_loc_df %>%
    dplyr::filter(Pathology >= pathology_threshold1 & Signature >= sig_pval_threshold1) %>%
    dplyr::arrange(desc(Pathology), desc(round(Signature, 1)), user.Sample_Name) %>%
    extract_topn_my(topn_param=topn_param, category=category)
  
  
  # Pathology low, signature high
  category <- "Lo-Hi"
  spot_lists[[paste(category, s)]] <- meta_loc_df %>%
    dplyr::filter(Pathology < pathology_threshold2 & Signature >= sig_pval_threshold1) %>%
    dplyr::arrange(Pathology, desc(round(Signature, 1)), user.Sample_Name) %>%
    extract_topn_my(topn_param=topn_param, category=category)
  
  # Pathology high, signature low
  category <- "Hi-Lo"
  spot_lists[[paste(category, s)]] <- meta_loc_df %>%
    dplyr::filter(Pathology >= pathology_threshold1 & Signature < sig_pval_threshold2) %>%
    dplyr::arrange(desc(Pathology), round(Signature, 1), user.Sample_Name) %>%
    extract_topn_my(topn_param=topn_param, category=category)
  
  # Pathology low, signature low
  category <- "Lo-Lo"
  spot_lists[[paste(category, s)]] <- meta_loc_df %>%
    dplyr::filter(Pathology < pathology_threshold2 & Signature < sig_pval_threshold2) %>%
    dplyr::arrange(Pathology, round(Signature, 1), user.Sample_Name) %>%
    extract_topn_my(topn_param=topn_param, category=category)
}


# Tibble of selected spots
spot_df <- dplyr::bind_rows(spot_lists)

# Write to file
filename <- paste0(output_path, paste(c("table", sig_select, pathology_select), collapse="_"), ".txt")
readr::write_delim(spot_df, filename, delim="\t")


## VISUALIZE SPOTS ----


# Add categories to the Seurat object
data_seurat@meta.data[classifier_name] <- "-"
for(spot_name in names(spot_lists)){
  data_seurat@meta.data[spot_lists[[spot_name]]$Coordinate, classifier_name] <- gsub(" .+", "", spot_name)
}
table(data_seurat@meta.data[classifier_name])

# Assign these colors to the categories
cols <- list("Hi-Hi" = "red", "Hi-Lo" = "orange", "Lo-Hi" = "darkblue", "Lo-Lo" = "grey", "-" = "white")

# Plot selected spots
# For a given variable, plot PCA, UMAP, and spatial distributions
filename <- paste0(output_path, paste(c("spatial", sig_select, pathology_select), collapse="_"))
Seurat_pca_umap_spatial_my(data_seurat, classifier_name, filename, cols=cols)


## CREATE SWATCHES ----


# Concentriq image information
image_df <- meta_init_df %>% 
  dplyr::select(user.Sample_Name, user.Concentriq_Image_ID, user.Concentriq_Repo_ID) %>% 
  unique()


# Keep log of swatches
log_file <- paste0(output_path, "log.txt")
write("Cluster\tBarcode\tRow (Y)\tCol (X)", log_file, append=FALSE)

# Go through samples and extract swatches
for(group in sample_list){
  
  for(category in unique(spot_df$Category)){
    
    # Grab swatches in the group
    spot_loc_df <- spot_df %>%
      dplyr::filter(Category == category & !!rlang::sym(swatch_group) == group)
    
    
    swatch_list <- list()
    for(sample in unique(spot_loc_df$user.Sample_Name)){
      
      # Extract spot information for a specific group and specific sample
      coordinate_list <- spot_loc_df %>%
        dplyr::filter(user.Sample_Name == sample) %>%
        dplyr::pull(Coordinate)
      
      # Concentriq image ID
      concetriq_id <- image_df %>%
        dplyr::filter(user.Sample_Name == sample) %>%
        dplyr::pull(user.Concentriq_Image_ID)
      
      if(length(coordinate_list) >= 1){
        
        # Collect a patch for every coordinate
        for(coord in coordinate_list){
          xy <- data_seurat@images[[sample]]@coordinates[coord, c("imagerow", "imagecol")]
          rad <- data_seurat@images[[sample]]@scale.factors$spot_diameter_fullres / 2
          
          # Create patches
          swatch_list[[coord]] <- fetch_image_patch(in_concentriq_id=concetriq_id, 
                                                    in_coord_px=c(xy[[2]], xy[[1]]), 
                                                    in_radius_px=rad)
          
          write(paste(c(group, coord, xy), collapse="\t"), log_file, append=TRUE)
        }
      }
    }
    
    # Visualize patches
    if(length(swatch_list) > 0){
      title <- paste0(group, ".  ", paste(pathology_select, sig_select, collapse="; "), ".  ", category)
      filename <- paste0(output_path, paste(c("swatch", group, pathology_select, sig_select, category), collapse="_"), ".png")
      
      patch_image_list <- combine_and_vis_patches_my(swatch_list, title=title, filename=filename)
      
    }
    
  }
}


## UPDATE THE SEURAT OBJECT ----


if(no_rds_output == FALSE){
  # Write the updated Seurat object
  saveRDS(data_seurat, file=input_path)
}
