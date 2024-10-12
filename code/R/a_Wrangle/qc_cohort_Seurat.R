# Author: Anna Lyubetskaya. Date: 20-12-28

# Create post-filtering QC visualizations using the final Seurat object

# Plot mean post-filter QCs and corresponding distributions:
# "nCount_Spatial", "nFeature_Spatial", "nCount_SCT", "nFeature_SCT", "mito_percent", "ribo_percent"


## SETUP ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_tibble.R")
source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_10X_vis.R")


## PARAMETERS ----


# Select this cohort
cohort_name <- "PDAC108_path14_samples"
cohort_regex <- "PDAC_"

# Plotting parameters
# user.Block_ID
meta_col_x <- "user.Block_ID"
meta_col_f <- "user.Block_ID"

# List of meta data fields to parse
field_list <- c("nCount_Spatial", "nFeature_Spatial", "nCount_SCT", "nFeature_SCT",
                paste0("Pathology.", c("IntestineAdj", "MuscleAdj", "NormalAdj",  "LymphNode", "TLSAggregate", "TLSImmature", "TLSMature", 
                                       "Adipose", "Vessel", "Muscle", "Nerve", "NonEpi", "BenignEpi", "Blood", "ExoEndo", "LuminalNec", "Tumor"), ".percent"))
col_names <- paste0("qc.mean.", field_list)

# Boolean that controls whether to do variable genes analysis or not (this analysis takes a while)
do_var_genes <- FALSE

# Parameter to use as clusters
resolution <- "Pathology.Group"


## PATHS ----


# Input folder
input_path <- "XXXX"

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_name, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Find all files following the regex
file_list <- dir(input_path, pattern=paste0(cohort_regex, ".*.rds"), full.names=TRUE)

# Ingest a set of RDS Seurat objects
seurat_list <- read_rds_list_simple_my(file_list)
names(seurat_list) <- gsub("_annotated", "", names(seurat_list))


## VISUALIZE UMI DISTRIBUTIONS ----


if(do_var_genes == TRUE){
  for(name in names(seurat_list)){
    
    data_seurat <- seurat_list[[name]]
    
    # Identify variable genes
    data_seurat <- Seurat::FindVariableFeatures(data_seurat, selection.method = "vst", nfeatures = 2000, assay="Spatial")
    
    # Identify the 10 most highly variable genes
    top10 <- head(Seurat::VariableFeatures(data_seurat), 20)
    
    # Plot variable features with and without labels
    plot1 <- Seurat::VariableFeaturePlot(data_seurat, assay="Spatial")
    plot2 <- Seurat::LabelPoints(plot = plot1, points = top10, repel = TRUE)
    
    # Save to file
    filename <- paste0(output_path, "scatter_variable_genes_", name)
    write_plot2file_my(plot2, filename, num_row=2, num_col=2)
    
    
    # Extract UMI expression data as a long tibble
    data_df <- seurat_expression_to_long_tibble_my(data_seurat, assay="Spatial", slot="data") %>%
      dplyr::filter(Spatial_data > 0) %>%
      dplyr::mutate(IsVariable = Symbol %in% Seurat::VariableFeatures(data_seurat))
    
    # Histogram of counts
    filename <- paste0(output_path, "hist_umi_", name)
    create_hist_plot_my(data_df, x_label="Spatial_data", fill_label="IsVariable", intercept=c(0, 1), 
                        binwidth=1, filename=filename, labels=NULL, add_density=FALSE, log_scale=TRUE)
  }
}


## VISUALIZE UMI DISTRIBUTIONS BY CLUSTER ----


for(name in names(seurat_list)){
  
  data_seurat <- seurat_list[[name]]
  
  if(resolution %in% colnames(data_seurat@meta.data)){
    # UMI count spatial plot
    p1 <- spatial_feature_plot_my(data_seurat, feature="nCount_Spatial", min.cutoff="q1", max.cutoff="q99", 
                                  title=paste0(name, "\nUMIs per spot"))
    
    # Feature count spatial plot
    p2 <- spatial_feature_plot_my(data_seurat, feature="nFeature_Spatial", min.cutoff="q1", max.cutoff="q99", 
                                  title=paste0(name, "\nFeatures per spot"))
    
    # Cluster spatial plot
    p3 <- spatial_dim_plot_my(data_seurat, group.by=resolution, title=paste0(name, "\nClustering"))
    
    # Identify colors from the previous plot to match in the subsequent plots
    num_colors <- find_col_num_my(data_seurat, resolution)
    cols <- define_cols_my(n=num_colors)
    
    # Boxplot of UMIs by cluster
    p4 <- create_box_plot_my(data_seurat@meta.data, x_label=resolution, y_label="nCount_Spatial", 
                             fill_label=resolution, filename=NULL) +
      scale_fill_manual(values=cols)
    
    # Boxplot of UMIs by cluster
    p5 <- create_box_plot_my(data_seurat@meta.data, x_label=resolution, y_label="nFeature_Spatial", 
                             fill_label=resolution, filename=NULL) +
      scale_fill_manual(values=cols)
    
    # Write figure to file
    filename <- paste0(output_path, "cluster_", name)
    write_plot2file_my(patchwork::wrap_plots(list(p1, p2, p3, p4, p5), nrow=1), filename, num_col=5, num_row=1)
  }
}


## WRANGLE QC DATA ----


qc_mean_list <- list()
qc_dist_list <- list()
for(name in names(seurat_list)){
  # Select Seurat dataset
  data_seurat <- seurat_list[[name]]
  
  # Count number of spots in each object
  data_seurat@misc[["SpotCountFinal"]] <- ncol(data_seurat)
  
  
  # Extract data from the misc slot of a Seurat object
  qc_mean_list[[name]] <- tibble::tibble(Sample = name,
                                         Parameter = names(data_seurat@misc),
                                         Value = unlist(unname(data_seurat@misc)))
  
  # Extract select spot meta data fields
  qc_dist_list[[name]] <- tibble::as_tibble(data_seurat@meta.data[,field_list]) %>%
    dplyr::mutate(Sample = name) %>%
    dplyr::select(dplyr::all_of(c("Sample", field_list)))
  
}

# Tibble of all values extracted from the misc slot of Seurat objects
qc_mean_wide_df <- dplyr::bind_rows(qc_mean_list) %>%
  df_long2wide_my(rows="Sample", cols="Parameter", value="Value") #%>%
#dplyr::mutate_at(vars(dplyr::all_of(col_names)), as.numeric)

# Save stats to a file
filename <- paste0(output_path, "stats_misc.txt")
readr::write_delim(qc_mean_wide_df, filename, delim="\t")

# Tibble of select spot meta data fields
qc_dist_long_df <- dplyr::bind_rows(qc_dist_list) %>%
  df_wide2long_my(key="Parameter", val="Value") %>%
  dplyr::inner_join(qc_mean_wide_df, by="Sample") %>%
  dplyr::arrange(Parameter)

# Calculate additional stats from the meta.data
qc_dist_stat_df <- qc_dist_long_df %>%
  dplyr::group_by(Sample, Parameter) %>%
  dplyr::summarise(Mean = mean(Value)) %>%
  df_long2wide_my(rows="Sample", cols="Parameter", value="Mean") %>%
  dplyr::inner_join(qc_mean_wide_df %>%
                      dplyr::select(Sample, SpotCountFinal), by="Sample")

# Save stats to a file
filename <- paste0(output_path, "stats_meta_mean.txt")
readr::write_delim(qc_dist_stat_df, filename, delim="\t")


## VISUALIZE QC MEANs ----


# Go through QC field groups and bar plot them
for(metric in paste0("qc.mean.", field_list)){
  
  labels <- c(meta_col_x, metric, metric)
  filename <- paste0(output_path, "box_", metric)
  
  if(metric %in% colnames(qc_mean_wide_df)){
    qc_mean_wide_df[[metric]] <- as.numeric(qc_mean_wide_df[[metric]])
    
    p <- create_box_plot_my(qc_mean_wide_df, x_label=meta_col_x, y_label=metric, fill_label=meta_col_f, 
                            labels=labels, with_dots=TRUE, filename=NULL, reorder_x=FALSE)
    
    write_plot2file_my(p, filename, num_row=1, num_col=2)
  }
}


## VISUALIZE QC DISTRIBUTIONS ----


# Go through metrics and box plot them
for(metric in unique(qc_dist_long_df$Parameter)){
  
  filename <- paste0(output_path, "cohort_qc_box_", metric)
  p <- create_box_plot_my(qc_dist_long_df %>%
                            dplyr::filter(Parameter == metric) %>%
                            dplyr::mutate(Value = as.numeric(Value)), 
                          x_label="Sample", y_label="Value", fill_label=meta_col_f, filename=NULL)
  
  write_plot2file_my(p, filename, num_row=1, num_col=2)
}
