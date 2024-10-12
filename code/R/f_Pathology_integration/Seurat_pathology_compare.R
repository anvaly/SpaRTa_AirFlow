# Author: Anna Lyubetskaya. Date: 21-02-22

# Compare pathology and clustering annotations
# "single-identity spot" is a spot that falls into a single pathology compartment


## SETUP ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_signatures.R")
source("code/utils/utils_heatmap.R")

source("code/R/Utils/utils_10X_vis.R")

source("code/R/b_Cluster_Biomarkers/Seurat_cluster_utils.R")


## PARAMETERS ----


# Sample / Cohort name
sample_name <- "PDAC108_path14_5K_harmony"

# Which field to use for clustering
clust_res_col <- "integrated_snn_res.0.1"  # integrated_snn_res.0.3

# Field to use in plotting
field_list <- c("orig.ident")  # "user.Sample_Name"

# A list of classes to count towards the dominant class label for each spot
# The order of this vector is important: it decides how to resolve ties between dominant classes
annotation_classes <- paste0("Pathology.", c("IntestineAdj", "MuscleAdj", "NormalAdj",  "LymphNode", "TLSAggregate", "TLSImmature", "TLSMature", 
                                             "Adipose", "Vessel", "Muscle", "Nerve", "NonEpi", "BenignEpi", "Blood", "ExoEndo", "LuminalNec", "Tumor"), ".percent")

# Plot additional figures that only include "single identity spots"
extra_plots <- FALSE
# Threshold to identify a single-identity spot
pathology_percent <- 90
# Min number of spots in a group
spot_num <- 10

# User defined class order for plotting
# NULL or c(...)
class_order <- c("Tumor", "BenignEpi", "NormalAdj", "IntestineAdj", "ExoEndo",
                 "LuminalNec", "NonEpi", "Muscle", "MuscleAdj", "Nerve",
                 "Vessel", "Blood", "LymphNode",
                 "TLSAggregate", "TLSImmature","TLSMature", "Adipose")

# User defined colors when it has to be fully prescriptive
cols <- c("#9A2626", "#4CB0E1", "#AEE0EB", "#AEE0EA", "#056DB5", "#939393",
          "#E5E5AA", "#A0D5B5", "#A0D5B6", "#FAE573", "#9B5E2C",
          "#9B5E2D", "#DF7126", "#FF9F2C", "#FF9F2D", "#FF9F2E", "#CCCCCB")
names(cols) <- class_order


# Cluster res order
# Harmony
cluster_order <- as.character(c(3, 6, 4, 1, 8, 0, 5, 2, 7, 10, 11, 9))


## PATHS ----


# Input folder
input_init_path <- "XXXX"
rds_file <- paste0(input_init_path, sample_name, ".rds")

# Output folder
output_init_figs <- "XXXX"
output_figs <- paste0(output_init_figs, sample_name, "/")

# Create output folders
dir.create(output_init_figs, showWarnings = FALSE)
dir.create(output_figs, showWarnings = FALSE)


## INGEST DATA ----


# Ingest Seurat data
data_seurat <- readRDS(rds_file)

# Add absent parameters in case non-integrated cohort
for(f in field_list){
  if(!f %in% colnames(data_seurat@meta.data)){
    data_seurat@meta.data[f] <- data_seurat@misc[[f]]
  }
}

# Define pathology class order if not defined by user
if(is.null(class_order)){
  class_order <- gsub("Pathology.|.percent", "", annotation_classes)
}


## WRANGLE DATA ----


# Seurat meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

# Preferred clustering resolution
if(is.null(clust_res_col)){
  clust_res <- gsub(".0$", "", data_seurat@misc[["user.Clustering"]])
} else{
  clust_res <- clust_res_col
}

# # Factorize the meta data to plot in a specific order if defined by user
# if(!is.null(class_order)){
#   meta_order <- match(class_order, unique(meta_df[[clust_res]]))
#   meta_df[[clust_res]] <- factor(meta_df[[clust_res]], 
#                                  levels=unique(meta_df[[clust_res]])[meta_order])
# }

# Pathology classes names
path_cols <- intersect(colnames(meta_df)[grep("^Pathology.*percent$", colnames(meta_df))], annotation_classes)

# Extract pathology classes and Seurat clustering data
clust_wide_df <- meta_df %>%
  dplyr::select(dplyr::all_of(c("Coordinate", field_list, clust_res, path_cols)))

# Wide to long cluster information
clust_df <- clust_wide_df %>%
  df_wide2long_my(key="Pathology.Group", val="Percent", start_col=3+length(field_list)) %>%
  dplyr::mutate(Pathology.Group = gsub("Pathology.|.percent", "", Pathology.Group))

# Sort coordinates by the majority pathology group and corresponding percentage
# Find the most abundant pathology class for each spot
clust_max_df <- clust_df %>%
  dplyr::group_by(Coordinate) %>%
  dplyr::arrange(desc(Percent), factor(Pathology.Group, levels=class_order)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Pathology.Group, Percent)


# Basic checks
setdiff(unique(clust_df$Pathology.Group), class_order)
setdiff(class_order, unique(clust_df$Pathology.Group))
unique(clust_max_df$Pathology.Group)


# Write dominant spot identity to file
filename <- paste0(output_figs, "table_path_id_", sample_name, ".txt")
readr::write_delim(clust_max_df, filename, delim="\t")


# Sort coordinates by the majority pathology group and corresponding percentage
clust_df$Coordinate <- factor(clust_df$Coordinate, levels=clust_max_df$Coordinate)

# Select spots dominated by a single pathology class (single-identity spots)
spots_clean_df <- clust_max_df %>% 
  dplyr::filter(Percent >= pathology_percent)

# Summarize groups of single-identity spots
spots_clean_sum_df <- spots_clean_df %>%
  dplyr::group_by(!!rlang::sym(clust_res), Pathology.Group) %>%
  dplyr::summarise(Count = dplyr::n_distinct(Coordinate)) %>%
  dplyr::filter(Count >= spot_num)

# Barcodes for single-identity spots
barcode_list <- spots_clean_df %>%
  dplyr::inner_join(spots_clean_sum_df, by=c(clust_res, "Pathology.Group")) %>%
  dplyr::pull(Coordinate)


## PLOT CLUSTERS VS PATHOLOGY ----


# Variable of interest
var <- "Pathology.Group"

# For a given variable, plot PCA, UMAP, and spatial distributions
output_name <- paste0(output_figs, "/cluster_", sample_name, "_", var)
cols <- Seurat_pca_umap_spatial_my(data_seurat, var, output_name, col_names=intersect(class_order, unique(data_seurat@meta.data[[var]])), 
                                   cols=cols, title=sample_name)

# Make sure that the barplot category order matches the spatial plot category order
clust_df["Pathology.Group"] <- factor(clust_df[["Pathology.Group"]], levels=names(cols))


if(nrow(clust_df) <= 10000){
  
  # Plot spots by pathology classes and Seurat clusters
  p <- create_bar_plot_my(clust_df, x_label="Coordinate", y_label="Percent", fill_label="Pathology.Group", 
                          facet_var=c(clust_res, "free"), position="stack", filename=NULL, 
                          labels=c("Spot", "Pathology Class, %", ""), reorder_x=TRUE) +
    scale_fill_manual(values=cols)
  
  filename <- paste0(output_figs, "pathology_", sample_name, "_bar_spot_clust_by_path")
  write_plot2file_my(p, filename)
  
  
  if(extra_plots == TRUE){
    # Plot spots by pathology classes and Seurat clusters - only include single-identity spots
    p <- create_bar_plot_my(clust_df %>%
                              dplyr::filter(Coordinate %in% barcode_list), 
                            x_label="Coordinate", y_label="Percent", fill_label="Pathology.Group", 
                            facet_var=c(clust_res, "free"), position="stack", filename=NULL, 
                            labels=c("Spot", "Pathology Class, %", ""), reorder_x=TRUE) +
      scale_fill_manual(values=cols)
    
    filename <- paste0(output_figs, "pathology_select_", sample_name, "_bar_spot_clust_by_path")
    write_plot2file_my(p, filename)
  }
}


# Plot additional views broken by up by user-defined fields
for(field in c(field_list)){
  
  # Aggregate pathology classes by cluster
  clust_sum_df <- clust_df %>%
    tidyr::drop_na() %>%
    dplyr::group_by(Pathology.Group, !!rlang::sym(clust_res), !!rlang::sym(field)) %>%
    dplyr::summarise(Percent_Sum = sum(as.numeric(Percent))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!rlang::sym(clust_res), !!rlang::sym(field)) %>%
    dplyr::mutate(Cluster_Sum = sum(Percent_Sum)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Percent = round(Percent_Sum / Cluster_Sum * 100))
  
  filename <- paste0(output_figs, "pathology_", sample_name, "_", field, "_clust_by_path.txt")
  readr::write_delim(clust_sum_df, filename, delim="\t")
  
  # Plot pathology classes by Seurat clusters
  p <- create_bar_plot_my(clust_sum_df, x_label=clust_res, y_label="Percent", fill_label="Pathology.Group", 
                          facet_var=c(field, "free"), position="stack", filename=NULL, 
                          labels=c("Cluster", "Pathology Class, %", field)) +
    scale_fill_manual(values=cols)
  
  filename <- paste0(output_figs, "pathology_", sample_name, "_", field, "_bar_clust_by_path")
  write_plot2file_my(p, filename, 
                     num_row=length(unique(clust_sum_df[[field]])), 
                     num_col=round(length(unique(clust_sum_df[[clust_res]]))/10))
  
  
  if(extra_plots == TRUE){
    # Aggregate pathology classes by cluster - only include single-identity spots
    clust_sum_unique_df <- clust_df %>%
      dplyr::filter(Coordinate %in% barcode_list) %>%
      tidyr::drop_na() %>%
      dplyr::group_by(Pathology.Group, !!rlang::sym(clust_res), !!rlang::sym(field)) %>%
      dplyr::summarise(Percent_Sum = sum(as.numeric(Percent))) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(!!rlang::sym(clust_res), !!rlang::sym(field)) %>%
      dplyr::mutate(Cluster_Sum = sum(Percent_Sum)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Percent = round(Percent_Sum / Cluster_Sum * 100))
    
    filename <- paste0(output_figs, "pathology_select_", sample_name, "_", field, "_clust_by_path.txt")
    readr::write_delim(clust_sum_unique_df, filename, delim="\t")
    
    # Plot pathology classes by Seurat clusters
    p <- create_bar_plot_my(clust_sum_unique_df, x_label=clust_res, y_label="Percent", fill_label="Pathology.Group", 
                            facet_var=c(field, "free"), position="stack", filename=NULL, 
                            labels=c("Cluster", "Pathology Class, %", field)) +
      scale_fill_manual(values=cols)
    
    filename <- paste0(output_figs, "pathology_select_", sample_name, "_", field, "_bar_clust_by_path")
    write_plot2file_my(p, filename, 
                       num_row=length(unique(clust_sum_unique_df[[field]])), 
                       num_col=round(length(unique(clust_sum_df[[clust_res]]))/10))
  }
  
}


## Create a heatmap of groups by % pathology niche
params <- list(cell_value = "Percent",
               row_label = var, 
               col_label = clust_res_col, 
               distance = "pearson",
               row_annotation = NULL,
               col_annotation = NULL,
               range = seq(0,50,10),
               colors = RColorBrewer::brewer.pal(n=6, name="YlGnBu"),
               row_order = class_order,
               column_order = cluster_order)

clust_sum_df[[clust_res_col]] <- as.character(clust_sum_df[[clust_res_col]])
clust_sum_df[[var]] <- as.character(clust_sum_df[[var]])

filename <- paste0(output_figs, "hm_pathology_", sample_name, ".png")
create_heatmap_my(clust_sum_df, params, row_list=NULL, col_list=NULL, 
                  col_meta_df=NULL, row_meta_df=NULL, filename=filename,
                  width=6, height=6)
