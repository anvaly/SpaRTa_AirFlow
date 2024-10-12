# Author: Anna Lyubetskaya. Date: 22-09-09
# For a select set of spots find neighborhoods of different sizes and summarize their properties


## ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_image.R")
source("code/R/Utils/utils_10X_signatures.R")


## PARAMETERS ----


# Name of the analysis to use in folder/file names
run_name <- "pub"

# Path to processed Seurat data
# PDAC108_path14_5K_harmony, PDAC108_path14_harmonyepi_rpca_sct
sample_name <- "PDAC108_path14_5K_harmony"
sample_exclude <- NULL

# Name of the pathology field to filter by groups
pathology_criteria <- c("Pathology.BenignEpi.percent", 
                        "Pathology.LuminalNec.percent",
                        "Pathology.Tumor.percent")
pathology_name <- "Region.Epi"
pathology_threshold <- 0

# Name for the signatures to be plotted
sig_select <- c("PDAC.collisson.classical", "PDAC.collisson.exocrine","PDAC.collisson.quasimesenchymal",
                "PDAC.moffitt.activatedstroma","PDAC.moffitt.basal", "PDAC.moffitt.classical","PDAC.moffitt.normalstroma",
                "PDAC.P19.Acinar","PDAC.P19.Bcell","PDAC.P19.Ductal_1","PDAC.P19.Ductal_2",
                "PDAC.P19.Endocrine","PDAC.P19.Endothelial","PDAC.P19.Fibroblast",
                "PDAC.P19.Macrophage","PDAC.P19.Stellate","PDAC.P19.Tcell",
                "PDAC.CosMx.Mast","PDAC.CosMx.Plasma","PDAC.U.Nervous",
                "BMS.Pathway.IFNa","BMS.Pathway.IFNg",
                "BMS.Pathway.TGFB","BMS.Pathway.TNFa",
                "BMS.CL.Hypoxia", "Syng.U.State.Hypoxia_Hallmark",
                "PDAC.Elyada19.panCAF","PDAC.Elyada19.iCAF","PDAC.Elyada19.myCAF"
)

# Names of pathology fields to be plotted
pathology_select <- paste0("Pathology.", c("IntestineAdj", "MuscleAdj", "NormalAdj",  "LymphNode", "TLSAggregate", "TLSImmature", "TLSMature", 
                                           "Adipose", "Vessel", "Muscle", "Nerve", "NonEpi", "BenignEpi", "Blood", "ExoEndo", "LuminalNec", "Tumor"), ".percent")

# SCT and spot number thresholds for gene filtering
sct_threshold <- 0.5
spot_threshold <- 100

# Assay and slot to use to calculate signature score
assay <- "SCT"
slot <- "data"

# Feature # in a spot threshold
feature_threshold <- NULL


## PATHS ----


# Path to a signature file
sig_path <- "data/import/Signatures/signatures_pdac_20230720.txt"

# Location of pre-processed data
input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, run_name, "_", sample_name)

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

# Remove spots if necessary
if(!is.null(feature_threshold)){
  barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data$nFeature_SCT >= feature_threshold),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}


## CALCULATE SIGNATURE SCORES ----


## Calculate signature scores and add them to Seurat meta data

# Find genes abundant in this sample
gene_list <- unique(seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold, 
                                                    assay=assay, slot=slot))

# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select, 
                                            sig_length_min=1, sig_length_max=1000, ratio_threshold=0)

# Find all genes in signatures and their signature assignations
sig_invert_df <- invert_list_my(signature_list)

# Add signature scores to a Seurat object
data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.", assay=assay)

# Find column names in the Seurat object - Seurat can rename original signatures
sig_names <- colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))])


## WRANGLE META DATA ----


# Extract coordinate meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

# Create a joint pathology label
meta_df[[pathology_name]] <- matrixStats::rowMaxs(as.matrix(meta_df[pathology_criteria]))

# Add user.Sample_Name to meta data if absent: happens in stand-alone samples
if(!"user.Sample_Name" %in% colnames(meta_df)){
  meta_df["user.Sample_Name"] = data_seurat@misc$user.Sample_Name
}

# Select the compartment of interest
meta_df <- meta_df %>%
  dplyr::mutate(PathologyCompartment = !!rlang::sym(pathology_name) >= pathology_threshold)

# Threshold and plot gene expression signatures
for(s in sig_names){
  # Global signature stats
  sig_mean_global <- mean(meta_df[[s]])
  sig_median_global <- median(meta_df[[s]]) 
  sig_sd_global <- sd(meta_df[[s]])
  
  sig_threshold1 <- sig_median_global + sig_sd_global
  sig_threshold2 <- sig_median_global + sig_sd_global * 2
  
  cat(s, sig_mean_global, sig_median_global, sig_sd_global, sig_threshold1, sig_threshold2, "\n")
  
  # ID tumor categories
  meta_df <- meta_df %>%
    dplyr::mutate(SigTemp = !!rlang::sym(s) >= sig_threshold1)
  
  colnames(meta_df) <- gsub("^SigTemp$", paste0(s, "_LabelTrue"), colnames(meta_df))
  
  table(meta_df[[paste0(s, "_LabelTrue")]])
  spot_num <- table(meta_df[[paste0(s, "_LabelTrue")]])[["TRUE"]]
  spot_perc <- round(table(meta_df[[paste0(s, "_LabelTrue")]])[["TRUE"]] / nrow(meta_df) * 100)
  
  # Plot signature score distribution with global thresholds
  plot_title <- paste0(s, "\nSignature score histogram\nSpots TRUE = ", spot_num, " (", spot_perc, "%)")
  filename <- paste0(output_path, "/hist_sig_scores_", s, "_", sample_name)
  create_hist_plot_my(meta_df, x_label=s, fill_label=paste0(s, "_LabelTrue"),
                      intercept=c(sig_median_global, sig_mean_global, sig_threshold1, sig_threshold2), binwidth=0.01, add_density=FALSE,
                      filename=filename, labels=c(s, "Spot Number", plot_title))
}


# Threshold and plot pathology features
for(s in pathology_select){
  # Global pathology stats
  sig_mean_global <- mean(meta_df[[s]])
  sig_median_global <- median(meta_df[[s]]) 
  sig_sd_global <- sd(meta_df[[s]])
  sig_threshold1 <- 50
  
  cat(s, sig_mean_global, sig_median_global, sig_sd_global, sig_threshold1, "\n")
  
  # ID tumor categories
  meta_df <- meta_df %>%
    dplyr::mutate(SigTemp = !!rlang::sym(s) >= sig_threshold1)
  
  colnames(meta_df) <- gsub("^SigTemp$", paste0(s, "_LabelTrue"), colnames(meta_df))
  
  table(meta_df[[paste0(s, "_LabelTrue")]])
  
  if(TRUE %in% names(table(meta_df[[paste0(s, "_LabelTrue")]]))){
    spot_num <- table(meta_df[[paste0(s, "_LabelTrue")]])[["TRUE"]]
    spot_perc <- round(table(meta_df[[paste0(s, "_LabelTrue")]])[["TRUE"]] / nrow(meta_df) * 100)
    
    # Plot signature score distribution with global thresholds
    plot_title <- paste0(s, "\nPathology percent histogram\nSpots TRUE = ", spot_num, " (", spot_perc, "%)")
    filename <- paste0(output_path, "/hist_path_percent_", s, "_", sample_name)
    create_hist_plot_my(meta_df, x_label=s, fill_label=paste0(s, "_LabelTrue"),
                        intercept=c(sig_median_global, sig_mean_global, sig_threshold1), binwidth=1, add_density=FALSE,
                        filename=filename, labels=c(s, "Spot Number, log2", plot_title), log_scale=TRUE)
  }
}


## FIND CONCENTRIC NEIGHBORHOODS AND CALCULATE THEIR PROPERTIES ----


# Move coordinate to row names
meta_df <- meta_df %>%
  tibble::column_to_rownames("Coordinate")

# Find target barcode list to make calculations below faster
barcode_target_list <- meta_df %>%
  dplyr::filter(PathologyCompartment == TRUE) %>%
  rownames()


# Neighborhood radius
nb_size <- 1.5

# Tibble list of neighborhood properties
nb_global_list <- list()

cl <- parallel::makeCluster(parallel::detectCores() / 2)
doParallel::registerDoParallel(cl)

# Loop through distances and calculate mean signature scores and total pathology composition by feature by neighborhood
for(distance in seq(nb_size, nb_size*5, nb_size)){
  print(distance)
  
  barcode_list <- list()
  
  # Loop through images
  barcode_list <- foreach::foreach(im_loc = data_seurat@images, .combine='c') %dopar% {
    # for(im_loc in data_seurat@images){
    barcode_names <- rownames(im_loc@coordinates)
    
    spot_clusters <- seurat_spot_neigbors_my(image=im_loc, spot_start=distance-nb_size, 
                                             spot_distance=distance)
    
    # barcode_list <- c(barcode_list, sapply(spot_clusters, function(x) barcode_names[x]))
    barcode_loc_list <- sapply(spot_clusters, function(x) barcode_names[x])
    names(barcode_loc_list) <- barcode_names
    
    barcode_loc_list
  }
  
  # Subset to only target barcodes - saves time on calculation
  barcode_list <- barcode_list[barcode_target_list]
  
  # Check the average size of the neighborhoods
  barcode_num <- sapply(barcode_list, function(x) length(x))
  mean(barcode_num)
  
  # Calculate neighborhood pathology / signature scores
  nb_list <- list()
  nb_list <- foreach::foreach(feature = c(sig_names, pathology_select), .combine='c') %dopar% {
    # for(feature in c(sig_names, pathology_select)){
    nb_loc_list <- list(a = sapply(barcode_list, function(x) mean(meta_df[x, feature])))
    names(nb_loc_list) <- feature
    
    nb_loc_list
  }
  names(nb_list) <- c(sig_names, pathology_select)
  
  
  # Signature information of this nb to tibble
  nb_global_list[[distance]] <- tibble::as_tibble(nb_list)
  # Add coordinate names
  nb_global_list[[distance]][["Coordinate"]] <- barcode_target_list
  # Add neighborhood distance
  nb_global_list[[distance]][["NbDistance"]] <- distance
  
  ## Print output on every loop
  
  # Combine all neighborhood information into tibble
  nb_global_df <- dplyr::bind_rows(nb_global_list) %>%
    dplyr::relocate(NbDistance, Coordinate) %>%
    df_wide2long_my(key="Feature", val="ScoreMean", start_col=3)
  
  # Write tibbles to file
  filename <- paste0(output_path, "/table_features_", run_name, "_", distance, ".txt")
  readr::write_delim(nb_global_df, filename, delim="\t")
  
  gc()
}

parallel::stopCluster(cl)
gc()


## VISUALIZE NEIGHBORHOOD COMPOSITIONS ----


# Box plots of neighborhood properties - pathology
filename <- paste0(output_path, "/box_sig_", run_name)
create_box_plot_my(nb_global_df %>%
                     dplyr::filter(grepl("sig.", Feature)), 
                   x_label="Feature", y_label="ScoreMean", fill_label="NbDistance", 
                   facet_var=c("NbDistance", "fixed"), filename=filename, labels=NULL)

# Box plots of neighborhood properties - signatures
filename <- paste0(output_path, "/box_path_", run_name)
create_box_plot_my(nb_global_df %>%
                     dplyr::filter(grepl("Pathology.", Feature)), 
                   x_label="Feature", y_label="ScoreMean", fill_label="NbDistance", 
                   facet_var=c("NbDistance", "fixed"), filename=filename, labels=NULL)
