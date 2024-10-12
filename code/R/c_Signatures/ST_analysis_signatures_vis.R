# Author: Anna Lyubetskaya. Date: 22-09-02
# Create various visualizations of a set of signatures as well as specific genes within them


## ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_heatmap.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_signatures.R")


## PARAMETERS ----


# Name of the processed Seurat object
cohort_name <- "PDAC108_path14_harmonyepi_rpca_sct"
remove_string <- ""

# Name for the signature group being plotted
sig_group_name <- "pub_tumor"

# SCT and spot number thresholds for gene filtering
sct_threshold <- 0.5
spot_threshold <- 5

# Assay and slot to use to calculate signature score
assay <- "SCT"
slot <- "data"

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
                "PDAC.Elyada19.panCAF","PDAC.Elyada19.iCAF","PDAC.Elyada19.myCAF")

# Meta parameters select
# Names of pathology fields to be plotted
meta_select <-  paste0("Pathology.", c("IntestineAdj", "MuscleAdj", "NormalAdj",  "LymphNode", "TLSAggregate", "TLSImmature", "TLSMature", 
                                       "Adipose", "Vessel", "Muscle", "Nerve", "NonEpi", "BenignEpi", "Blood", "ExoEndo", "LuminalNec", "Tumor"), ".percent")

do_heatmap <- FALSE


## Plotting cluster order controls


# Set defaults
resolution <- "integrated_snn_res.0.1"
pt.size <- 0.02
cluster_order <- NULL
cols <- NULL


## Harmony full integration, res 0.1


if(cohort_name == "PDAC108_path14_5K_harmony"){
  
  # Pt size for PCA / UMAP
  pt.size <- 0.02
  
  # User defined cluster order
  cluster_order <- c(3, 6, 4, 1, 8, 0, 5, 2, 7, 10, 11, 9)
  
  # User defined color schema - when full control of vis is desired
  cols <- c("#056DB5", "#153C65", "#4CB0E1", "#9A2626", "#96257D", "#E5E5AA",
            "#A0D5B5", "#FF9F2C", "#DF7126", "#AEE0EA", "#D4EEF5", "#939393")
  
  resolution <- "integrated_snn_res.0.1"
}


## Harmony-defined epi niche integrated by RPCA, res 0.4


if(cohort_name == "PDAC108_path14_harmonyepi_rpca_sct"){
  
  # Pt size for PCA / UMAP
  pt.size <- 0.1
  
  # User defined cluster order
  # cluster_order <- c(12,11,3,10,0,4,6,14,9,1,8,5,2,7,13)
  cluster_order <- c(3, 10, 0, 4, 6, 9)
  
  # User defined color schema - when full control of vis is desired
  # cols <- c("#056DB5","#153C65",
  #           "#4CB0E1","#D08EB3","#F16666","#9A2626","#96257D","#360F2E","#681D59",
  #           "#E5E5AA","#5E5E39","#358E5B","#FF9F2C","#D66100","#939393")
  cols <- c("#4CB0E1","#D08EB3","#F16666","#9A2626","#96257D","#681D59")
  
  resolution <- "integrated_snn_res.0.4"
}


## Harmony-defined stroma niche integrated by RPCA, res 0.2


if(cohort_name == "PDAC108_path14_harmonystr_rpca_sct"){
  
  # Pt size for PCA / UMAP
  pt.size <- 0.1
  
  # User defined cluster order
  cluster_order <- c(11,7,0,3,1,4,8,5,2,6,10,9)
  
  # User defined color schema - when full control of vis is desired
  cols <- c("#056DB5","#9A2626","#E5E5AA","#B5B572","#5E5E39","#184928",
            "#A0D5B5","#358E5B","#FF9F2C","#D66100","#49392E","#939393")
  
  resolution <- "integrated_snn_res.0.2"
}


## Harmony-defined immune niche integrated by RPCA, res 0.2


if(cohort_name == "PDAC108_path14_harmonyimm_rpca_sct"){
  
  # Pt size for PCA / UMAP
  pt.size <- 0.1
  
  # User defined cluster order
  cluster_order <- c(10,16,9,0,1,4,3,7,14,2,11,12,13,8,5,6,15)
  
  # User defined color schema - when full control of vis is desired
  cols <- c("#056DB5","#153C65","#9A2626","#E5E5AA","#B5B572",
            "#5E5E39","#878545","#358E5B","#184928","#FBE99E","#FEDC4B",
            "#FAE573","#FFCE07","#FF9F2C","#D66100","#A45024","#939393")
  
  resolution <- "integrated_snn_res.0.4"
}


names(cols) <- cluster_order


## PATHS ----


# Path to a signature file
sig_path <- "data/import/Signatures/signatures_pdac_20230720.txt"

# Location of pre-processed data
input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_name, "_", sig_group_name, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Seurat data
data_seurat <- readRDS(input_path)

# Remove clusters if necessary
if(!is.null(cluster_order)){
  barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data[[resolution]] %in% cluster_order),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}


## CALCULATE SIGNATURE SCORES ----


## Calculate signature scores and add them to Seurat meta data

# Find genes abundant in this sample
gene_list <- seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold,
                                             assay=assay, slot=slot)

# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select,
                                            sig_length_min=1, sig_length_max=1000, ratio_threshold=0)

# Add signature scores to a seurat object
data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.", assay=assay)

# Find column names in the Seurat object - Seurat can rename original signatures
sig_names <- colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))])

# Find all genes in signatures and their signature assignations
sig_invert_df <- invert_list_my(signature_list)


## WRANGLE DATA ----


# User defined resolution from the object
if(is.null(resolution)){
  resolution <- data_seurat@misc$user.Clustering
}

# If user didn't define the cluster order, just use the default cluster order
if(is.null(cluster_order)){
  cluster_order <- unique(sort(data_seurat@meta.data[[resolution]]))
}

# Factorize the meta data to plot in a specific order
meta_order <- match(cluster_order, unique(data_seurat@meta.data[[resolution]]))
data_seurat@meta.data[[resolution]] <- factor(data_seurat@meta.data[[resolution]], levels=cluster_order)


# Extract coordinate meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

# Add user.Sample_Name to meta data if absent: happens in stand-alone samples
if(!"user.Sample_Name" %in% colnames(meta_df)){
  meta_df["user.Sample_Name"] = data_seurat@misc$user.Sample_Name
}

# Crete a long tibble for meta data
meta_long_df <- meta_df %>%
  dplyr::select(dplyr::all_of(c("Coordinate", "user.Sample_Name", resolution, sig_names, meta_select))) %>%
  df_wide2long_my(key="Sig_Name", val="Score", start_col=4) %>%
  dplyr::mutate(NbDistance = 0)


# Write tibbles to file
filename <- paste0(output_path, "/table_features_", sig_group_name, ".txt")
readr::write_delim(meta_long_df %>%
                     dplyr::rename(Feature = Sig_Name, ScoreMean = Score) %>%
                     dplyr::select(-user.Sample_Name, -!!rlang::sym(resolution)), 
                   filename, delim="\t")


## CLUSTERING PLOT ----


# Define color schema
if(is.null(cols)){
  cols <- define_cols_for_var_my(data_seurat, resolution, col_names=NULL)
}

# Plot a set of variables spatially
filename <- paste0(output_path, "/cluster_", cohort_name, "_", resolution)
Seurat_pca_umap_spatial_my(data_seurat, resolution, filename, col_names=NULL, cols=cols, 
                           title=cohort_name, pt.size=pt.size)


## DIM AND SPATIAL PLOT ----


# Cycle through signatures
for(sig in c(sig_names, meta_select)){
  
  data_seurat@meta.data[which(is.na(data_seurat@meta.data[sig])), sig] <- 0
  
  # Create a UMAP
  p <- feature_plot_my(data_seurat, sig, pt.size=pt.size, reduction=NULL, min.cutoff="q5", max.cutoff="q95")
  filename <- paste0(output_path, "/umap_", sig, "_", cohort_name)
  write_plot2file_my(p, filename, num_row=1, num_col=1)
  
  # Create spatial plots
  filename <- paste0(output_path, "/spatial_", sig, "_", cohort_name)
  batch_spatial_feature_plot_my(list(Cohort = data_seurat), sig, output_file=filename)
  
  
  # Select data for a particular feature
  meta_long_select_df <- meta_long_df %>%
    dplyr::filter(Sig_Name == sig)
  
  # Signature histogram
  filename <- paste0(output_path, "/hist_", sig, "_", cohort_name)
  create_hist_plot_my(meta_long_select_df, x_label="Score", fill_label="NbDistance", 
                      filename=filename, labels=c(sig, "Score", sig), 
                      log_scale=FALSE, add_density=TRUE,
                      binwidth=0.01, intercept=c(0,1))
  
  # Average signature score vs cluster
  filename <- paste0(output_path, "/box_", sig, "_", cohort_name)
  create_box_plot_my(meta_long_select_df, x_label=resolution, y_label="Score", fill_label=resolution,
                     facet_var=c("Sig_Name", "free_x"), filename=filename, labels=c(resolution, "Score", sig),
                     outlier_shape=NA, cols=cols)
  
  gc()
}


## SIGNATURE BY CLUSTER PLOT ----


# Breakdown boxplots by samples if less than 5 samples in the cohort
if(length(unique(meta_long_df$user.Sample_Name)) <= 5){
  # Average signature score vs cluster
  filename <- paste0(output_path, "/box_", cohort_name, "_", sig_group_name, "_sample")
  create_box_plot_my(meta_long_df, x_label=resolution, y_label="Score", fill_label="user.Sample_Name", 
                     facet_var=c("Sig_Name", "free_x"), filename=filename, labels=NULL, outlier_shape=NA)
}


## DOT PLOT ----


for(sig in names(signature_list)){
  if(!is.null(cluster_order)){
    cluster.idents <- FALSE
  } else{
    cluster.idents <- TRUE
  }
  
  p <- Seurat::DotPlot(data_seurat, features = signature_list[[sig]], group.by = resolution,
                       assay = assay, col.min = 0, col.max = 1, cluster.idents = cluster.idents) +
    Seurat::RotatedAxis()
  
  num_col <- round(length(signature_list[[sig]]) / 10)
  if(num_col < 1){
    num_col <- 1
  }
  
  filename <- paste0(output_path, "/dot_", sig, "_", cohort_name)
  write_plot2file_my(p, filename, num_row=1, num_col=num_col)
}


## HEATMAP PLOT ----


if(do_heatmap == TRUE){
  # Switch slots for the heatmap to use scaled data
  slot <- "scale.data"
  
  # Values to plot
  cell_value <- paste0(assay, "_", slot)
  
  # Extract expression matrix, transform it into a long tibble, and calculate z-score for each gene across all barcodes
  data_df <- seurat_expression_to_long_tibble_my(data_seurat, assay=assay, slot=slot)
  
  # Setup heatmap parameters
  if(slot == "scale.data"){
    range <- c(0, 1)
    colors <- c("#55185D", "#FFD524")
  } else{
    mean_val <- mean(data_df[[cell_value]])
    sd_val <- sd(data_df[[cell_value]])
    min_val <- min(data_df[[cell_value]])
    max_val <- max(data_df[[cell_value]])
    
    range <- c(ceiling(min_val), round(mean_val + sd_val, 1))
    colors <- c("#55185D", "#FFD524")
  }
  
  # Setup heatmap parameters
  params <- list(
    cell_value = cell_value,
    row_label = "Coordinate", 
    col_label = "Symbol", 
    distance = "pearson",
    row_annotation = c("user.Sample_Name", resolution, "Pathology.Group"),
    col_annotation = c("Sig_Name"),
    range = range,
    colors = colors
    #column_km = length(unique(meta_df[[resolution]])),
    #row_km = length(signature_list)
  )
  
  # Create a heatmap of DE biomarkers
  filename <- paste0(output_path, "/hm_", cohort_name, "_", sig_group_name, ".png")
  
  # Plot the heatmap
  hm <- create_heatmap_my(data_df, params, row_list=NULL, col_list=sig_invert_df$Symbol, 
                          row_meta_df=meta_df, col_meta_df=sig_invert_df, filename=filename)
}
