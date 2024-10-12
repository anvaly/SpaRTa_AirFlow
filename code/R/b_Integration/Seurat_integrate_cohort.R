# Author: Anna Lyubetskaya. Date: 21-01-18
# Script for cohort integration:
# Gather a list of Seurat RDS objects, merge, recalculate SCT, output to RDS


## SETUP ENVIRONMENT ----


# Install MAST
if(!"MAST" %in% rownames(installed.packages()) || packageVersion("MAST") < 1.2){
  BiocManager::install("MAST")
}

# Update Seurat and glmGamPoi
if(packageVersion("Seurat") < 4.1){
  update.packages("Seurat")
}

if(packageVersion("glmGamPoi") < 1.6){
  install.packages("glmGamPoi")
}

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_signatures.R")

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_10X_cohort.R")
source("code/R/Utils/utils_10X_dea.R")
source("code/R/Utils/utils_10X_vis.R")

source("code/R/b_Cluster_Biomarkers/Seurat_cluster_utils.R")

set.seed(777)


## PARAMETERS INPUT ----


# Cohort name
cohort_name <- "PDAC108_path14_harmonyimm"
# Cohort regex to find inputs
cohort_regex <- ""

# Exclude specific samples from the analysis
samples_exclude <- c("PDAC_1275301B_s2", "PDAC_1275301B_s4", "PDAC_1275301B_s8", "PDAC_1275301B_s6", "PDAC_1275321B_s1", 
                     "PDAC_1275301B_s5", "PDAC_1275233B", "PDAC_1275301B_s7",
                     "PDAC_1275321B_s2", "PDAC_S158925_ROI2_s1", "PDAC_759104PB", "PDAC_1275309B_s2",
                     "PDAC_Pt1_ROI2", "PDAC_S158926_ROI3_s1", "PDAC_1275309B_s1")

# List of samples to use for integration reference
integration_ref_list <- c("PDAC_1255880B_ROI1_s1",
                          "PDAC_759110QB_ROI3_s2",
                          "PDAC_E2547_ROI2",
                          "PDAC_Pt1_ROI3",
                          "PDAC_Pt11_ROI3",
                          "PDAC_Pt2_ROI3",
                          "PDAC_Pt3_ROI2_s1",
                          "PDAC_Pt3_ROI3_s1",
                          "PDAC_Pt5_ROI2",
                          "PDAC_Pt6_ROI3",
                          "PDAC_Pt6_ROI4",
                          "PDAC_Pt8_ROI2",
                          "PDAC_S158915_ROI2_s1",
                          "PDAC_S158924_ROI2_s1",
                          "PDAC_S158926_ROI1_s1",
                          "PDAC_S158928_ROI1_s1")

# User clustering resolution
resolution_user <- "integrated_snn_res.0.4"

# List of misc fields to keep in meta data
misc_field_list <- NULL

# Color UMAPs by these categorical and continuous variables
color_field_list <- c("user.Sample_Name", "user.Region_ID", "user.Block_ID", 
                      "user.Area", "user.Slide", "user.Tissue_Source", 
                      "Pathology.Group")


## PARAMETERS TECHNICAL ----


# Filter spots that have too few features after SCT normalization on ingestion
feature_filter <- 1000

# Select integration method: simple merge without any change to SCT slot or sct integration
# Tried setting up SCVI - can't make it work
integration_method <- "rpca_sct"  # merge, cca_sct, rpca_sct, precast, harmony (scvi attempted but not achieved)

# Perform or skip FindMarkers step
do_findmarkers <- TRUE
# Perform or skip Find Conservative Markers step
do_findmarkers_conservative <- FALSE

# Save resulting RDS
do_save_rds <- TRUE

# Set parameters for the clustering analysis
params <- cluster_params_my()

# Add latent variables to DEA
params[["latent_vars"]] <- c("user.Sample_Name")

# Variables to regress with SCTransform
vars_to_regress <- c("mito_percent", "ribo_percent", params[["latent_vars"]])  # "mito_percent", "ribo_percent"

# Number of anchor features in the SCT CCA/RPCA integration method
feature_num <- 5000

# PRECAST requires an explicit setting of the number of clusters expected
precast_K <- 30

# Harmony parameters
# https://github.com/immunogenomics/harmony/issues/24
theta <- rep(4, length(params[["latent_vars"]]))  # cluster diversity enforcement. larger values encourage more diverse clusters. default theta = 2
lambda <- rep(1, length(params[["latent_vars"]]))  # ridge regression penalty. Smaller values result in more aggressive correction. default lambda = 1
nclust_harmony <- 2000  # number of clusters harmony uses as anchors. default nclust = 100


## FINISH SETUP ENVIRONMENT ----


# Install PRECAST - need at least version 1.5
if(!"harmony" %in% rownames(installed.packages()) && integration_method == "harmony"){
  BiocManager::install("harmony")
}

# Install PRECAST - need at least version 1.5
if(!"PRECAST" %in% rownames(installed.packages()) && integration_method == "precast"){
  remotes::install_github("feiyoung/PRECAST")  # upgrade="always"
  packageVersion("PRECAST")
}


## PATHS ----


# File of barcodes to use for integration
# E.g.: NULL or "XXXX"
barcode_selection_file <- "XXXX"

# Extended cohort name
cohort_name_full <- paste0(cohort_name, "_", integration_method)

# Location of annotated Seurat objects
input_paths <- c("XXXX")

# Output information
output_path <- "XXXX"
output_init_figs <- "XXXX"
output_figs <- paste0(output_init_figs, cohort_name_full)

# The integrated RDS output location
cohort_filename <- paste0(output_path, cohort_name_full, ".rds")

# Create output folder
dir.create(output_path, showWarnings = FALSE)
dir.create(output_init_figs, showWarnings = FALSE)
dir.create(output_figs, showWarnings = FALSE)


# Define unspecific probes
unspecific_probes_file <- NULL


## INTEGRATE DATA ----


if(!file.exists(cohort_filename)){
  
  ## INGEST DATA ----
  
  
  # Read barcode file if defined
  barcode_list <- NULL
  if(!is.null(barcode_selection_file)){
    barcode_list <- readr::read_delim(barcode_selection_file, delim="\t") %>%
      dplyr::pull(Coordinate)
  }
  
  # Read unspecific probes if defined
  if(!is.null(unspecific_probes_file)){
    unspecific_probes_list <- readr::read_delim(unspecific_probes_file, delim="\t")
  }
  
  # Find all files following the regex
  file_list <- unname(unlist(sapply(input_paths, function(p) dir(p, pattern=paste0(cohort_regex, ".*.rds"), full.names=TRUE))))
  
  # Remove user defined samples
  if(!is.null(samples_exclude)){
    file_list <- file_list[which(!grepl(paste0(samples_exclude, collapse="|"), file_list))]
  }
  
  cat("Seurat RDS objects found =", length(file_list), "\n")
  print(file_list)
  
  # Ingest a set of RDS Seurat objects
  seurat_list <- read_rds_list_simple_my(file_list, feature_filter=feature_filter, barcode_list=barcode_list)
  print(sum(sapply(seurat_list, function(x) ncol(x))))
  
  
  ## WRANGLE DATA ----
  
  
  # Remove unspecific probes - to be moved earlier in the logic  
  if(!is.null(unspecific_probes_file)){
    for(s in names(seurat_list)){
      seurat_list[[s]] <- subset(seurat_list[[s]], features=setdiff(rownames(seurat_list[[s]]@assays$Spatial@counts), unspecific_probes_list))
    }
  }
  
  # Update image names and remove previous clustering
  for(name in names(seurat_list)){
    # names(seurat_list[[name]]@images) <- gsub("-", ".", name)
    seurat_list[[name]]@meta.data <- seurat_list[[name]]@meta.data[which(!grepl("SCT_snn_res.", names(seurat_list[[name]]@meta.data)))]
  }
  
  # Move misc parameters to meta data slot
  seurat_list <- seurat_misc_to_meta(seurat_list, misc_field_list=misc_field_list)
  
  # Integrate all data in the cohort
  if(integration_method == "merge"){
    # Simply merge data into a single object
    data_combo_seurat <- seurat_merge_my(seurat_list)
  } else if(integration_method == "precast"){
    
    # Run PRECAST
    seuInt <- precast_integrate_my(seurat_list, feature_num=feature_num, K=precast_K)
    
    # Simply merge data into a single object
    data_combo_seurat <- seurat_merge_my(seurat_list)
    
    # Add PRECAST data to the combo object
    data_combo_seurat@assays[["integrated"]] <- seuInt@assays[["PRE_CAST"]]
    
    # Set var.features because PRECAST doesn't
    data_combo_seurat@assays[["integrated"]]@var.features <- rownames(seuInt@assays[["PRE_CAST"]])
    
    # Add the final clusters
    data_combo_seurat@meta.data[[paste0("precast_snn_res.", precast_K)]] <- seuInt@meta.data$cluster
    
    # Add PRECAST reductions
    data_combo_seurat@reductions <- seuInt@reductions
    
  } else if(integration_method == "scvi"){
    
    # seurat_scvi_integration <- function()
    
  } else if(integration_method == "harmony"){
    
    # Simply merge data into a single object
    data_combo_seurat <- seurat_merge_my(seurat_list)
    rm(seurat_list)
    gc()
    
    # Perform SCTransform normalization
    data_combo_seurat <- Seurat::SCTransform(data_combo_seurat, assay="Spatial", variable.features.n=feature_num,
                                             return.only.var.genes = TRUE, verbose = FALSE, vst.flavor = "v2")
    
    # merged_seurat_data <- Seurat::FindVariableFeatures(data_combo_seurat, nfeatures = feature_num)
    
    # PCA
    data_combo_seurat <- Seurat::RunPCA(data_combo_seurat, npcs=50, verbose = FALSE)
    
    # Run Harmony
    data_combo_seurat <- harmony::RunHarmony(data_combo_seurat, group.by.vars = params[["latent_vars"]], reduction = "pca",
                                             assay.use = "SCT", theta = theta, nclust = nclust_harmony, lambda = lambda,
                                             reduction.save = "harmony")
    
    gc()
    
  } else{
    data_combo_seurat <- seurat_sct_integrate_my(seurat_list, feature_num=feature_num, 
                                                 integration_method=integration_method, reference=integration_ref_list)
  }
  
  # Cleanup
  rm(seurat_list)
  
  # Update point size for visualizations
  data_combo_seurat@misc[["user.pt.size.factor"]] <- round(mean(unique(data_combo_seurat@meta.data$user.pt.size.factor)), 1)
  
} else{
  data_combo_seurat <- readRDS(cohort_filename)
}

gc()


## MANUALLY SET RESOLUTION ----


# Manually add cluster resolution
if(!is.null(resolution_user)){
  data_combo_seurat@misc[["user.Clustering"]] <- resolution_user
}


## NORMALIZE AND CLUSTER DATA ----


if(!file.exists(cohort_filename)){
  
  
  ## CLUSTER ANALYSIS ----
  
  
  # Perform SCTransform normalization
  if(integration_method != "harmony"){
    data_combo_seurat <- Seurat::SCTransform(data_combo_seurat, assay="Spatial", 
                                             vars.to.regress = vars_to_regress, variable.features.n=feature_num,
                                             return.only.var.genes = FALSE, verbose = FALSE, vst.flavor = "v2")
  }
  
  if(integration_method != "merge" && integration_method != "harmony"){
    # Switch default assay to integrated
    Seurat::DefaultAssay(data_combo_seurat) <- "integrated"
  }
  
  # Perform PCA, UMAP, FindNeighbors, FindClusters, and FindMarkers
  if(integration_method != "precast"){
    data_combo_seurat <- cluster_analysis_my(data_combo_seurat, params, cohort_name_full, output_figs)
    
    if(integration_method == "harmony"){
      colnames(data_combo_seurat@meta.data) <- gsub("SCT_snn_res.", "integrated_snn_res.", colnames(data_combo_seurat@meta.data))
    }
  } else{
    # Visualize PCA, UMAP, and tissue cluster distribution
    var <- paste0("precast_snn_res.", precast_K)
    filename_prefix <- paste0(cohort_name_full, "_", var)
    
    Seurat_pca_umap_spatial_my(data_combo_seurat, var, paste0(output_figs, "/", filename_prefix), 
                               col_names=NULL, cols=NULL, title=filename_prefix)
  }
}


## OUTPUT DATA ----


if(do_save_rds == TRUE){
  
  Seurat::DefaultAssay(data_combo_seurat) <- "SCT"
  
  # Save the full dataset to RDS
  saveRDS(data_combo_seurat, file = cohort_filename)
}


## SAMPLE MIXING FOR EACH RESOLUTION ----


# Column names
clust_col_names <- colnames(data_combo_seurat@meta.data[grep("_snn_res.", colnames(data_combo_seurat@meta.data))])


if(length(clust_col_names) > 0){
  
  # Extract clustering data and sample names
  meta_df <- tibble::as_tibble(data_combo_seurat@meta.data, rownames="Coordinate") %>%
    dplyr::select(dplyr::all_of(c("Coordinate", "user.Sample_Name", clust_col_names)))
  
  # Count number of spots by sample
  meta_sample_df <- meta_df %>%
    dplyr::group_by(user.Sample_Name) %>%
    dplyr::summarise(CountSample = dplyr::n_distinct(Coordinate))
  
  for(cl in clust_col_names){
    
    # Count number of spots by cluster
    meta_cluster_df <- meta_df %>%
      dplyr::group_by(!!rlang::sym(cl)) %>%
      dplyr::summarise(CountCluster = dplyr::n_distinct(Coordinate))
    
    # Count number of spots by cluster and sample
    # Number of spots by cluster and sample / Number of spots by cluster
    # Number of spots by cluster and sample / Number of spots by sample
    meta_loc_df <- meta_df %>%
      dplyr::group_by(user.Sample_Name, !!rlang::sym(cl)) %>%
      dplyr::summarise(Count = dplyr::n_distinct(Coordinate)) %>%
      dplyr::ungroup() %>%
      dplyr::inner_join(meta_cluster_df, by=cl) %>%
      dplyr::mutate(PercentCluster = round(Count / CountCluster * 100, 1)) %>%
      dplyr::inner_join(meta_sample_df, by="user.Sample_Name") %>%
      dplyr::mutate(PercentSample = round(Count / CountSample * 100, 1),
                    ClusterName = paste0(!!rlang::sym(cl), " (", CountCluster, ")"),
                    SampleName = paste0(user.Sample_Name, " (", CountSample, ")")) %>%
      dplyr::arrange(-CountCluster, -PercentSample, -PercentCluster)
    
    # Define custom colors
    cols <- define_cols_my(n=length(unique(meta_loc_df$ClusterName)))
    names(cols) <- sort(unique(meta_loc_df$ClusterName))
    
    filename <- paste0(output_figs, "/breakdown_table_", gsub("integrated_snn_res.", "", cl), "_", cohort_name_full, ".txt")
    readr::write_delim(meta_loc_df, filename, delim="\t")
    
    if(length(cols) <= 20){
      # Number of spots by cluster and sample / Number of spots by cluster
      filename <- paste0(output_figs, "/breakdown_clust_", gsub("integrated_snn_res.", "", cl), "_", cohort_name_full)
      p <- create_bar_plot_my(meta_loc_df, x_label="SampleName", y_label="PercentCluster", fill_label="ClusterName", 
                              facet_var=c("ClusterName", "free_y"), filename=filename, reorder_x=FALSE, cols=cols)
      
      # Number of spots by cluster and sample / Number of spots by sample
      p <- create_bar_plot_my(meta_loc_df, x_label="SampleName", y_label="PercentSample", fill_label="ClusterName", 
                              filename=NULL, reorder_x=FALSE, cols=cols)
      
      filename <- paste0(output_figs, "/breakdown_sample_", gsub("integrated_snn_res.", "", cl), "_", cohort_name_full)  
      write_plot2file_my(p, filename, num_row=1, num_col=8)
    }
    
  }
  
}


## VISUALIZE DATA ----


# Go through meta data variables of interest
for(var in color_field_list){
  
  # Define custom colors
  cols <- define_cols_for_var_my(data_combo_seurat, var)
  
  if(length(cols) > 1){
    # Turn pathology % into intervals with 10 increment
    if(grepl("^Pathology.*percent$|^CellCounts$", var)){
      # Cap values at a 100
      data_combo_seurat@meta.data[which(data_combo_seurat@meta.data[[var]] > 100), var] <- 100
      
      # Transform variable from continuous into categorical intervals
      interval_factors <- as.character(sort(as.numeric(unique(round(data_combo_seurat@meta.data[[var]]/10)*10))))
      data_combo_seurat@meta.data[[paste0(var, "Intervals")]] <- factor(round(data_combo_seurat@meta.data[[var]]/10)*10, 
                                                                        levels=interval_factors)
      
      var <- paste0(var, "Intervals")
      cols <- define_cols_for_var_my(data_combo_seurat, var, col_names=interval_factors, col_type="viridis")
    }
    # Turn signature score into intervals with 0.25 increment
    else if(grepl("^sig\\.", var)){
      interval_factors <- as.character(sort(as.numeric(unique(round(data_combo_seurat@meta.data[[var]]/2.5, 1)*2.5))))
      data_combo_seurat@meta.data[[paste0(var, "Intervals")]] <- factor(round(data_combo_seurat@meta.data[[var]]/2.5, 1)*2.5, 
                                                                        levels=interval_factors)
      
      var <- paste0(var, "Intervals")
      cols <- define_cols_for_var_my(data_combo_seurat, var, col_names=interval_factors, col_type="viridis")
    }
    
    if(length(cols) <= 100 && length(cols) >= 1){
      # For a given variable, plot PCA, UMAP, and spatial distributions
      output_name <- paste0(output_figs, "/cluster_", cohort_name_full, "_", var)
      Seurat_pca_umap_spatial_my(data_combo_seurat, var, output_name, cols=cols)
    }
  }
  
}


## PERFORM FIND MARKERS ----


if(do_findmarkers == TRUE){
  
  # Manually set the assay
  if(!is.null(resolution_user)){
    # Seurat::DefaultAssay(data_combo_seurat) <- "Spatial"
    data_combo_seurat <- Seurat::SetIdent(data_combo_seurat, value=resolution_user)
  }
  
  # Test every approach that allows latent variables
  for(test in c("wilcox")){  # c("wilcox", "LR", "negbinom", "poisson", "MAST)
    
    params[["test_use"]] <- test
    params[["assay"]] <- "SCT"
    params[["latent_vars"]] <- c("user.Sample_Name")
    
    filename_prefix <- paste0(cohort_name_full, "_", gsub("integrated_snn_res.", "", resolution_user), "_", test)
    
    
    ## FIND MARKERS ----
    
    
    var <- resolution_user
    
    # Find markers of each cluster against the rest
    markers_df <- seurat_find_markers_my(data_combo_seurat, assay=params[["assay"]], 
                                         find_all=TRUE, group.by=var, 
                                         min.pct=params[["pct_min"]], 
                                         logfc.threshold=params[["logfc_threshold"]],
                                         test.use=params[["test_use"]], 
                                         latent.vars=params[["latent_vars"]])
    
    markers_filt_df <- marker_analysis_my(markers_df, params, var, cohort_name_full, output_figs, filename_prefix)
    
    
    ## FIND CONSERVATIVE MARKERS ----
    
    
    if(do_findmarkers_conservative == TRUE){
      
      # Find the number of spots in each group where group is defined as a unique combination of Sample_Name and cluster ID
      spot_by_group_num <- table(data_combo_seurat@meta.data[c("user.Sample_Name", resolution_user)])
      
      # Find conserved markers for each cluster at user defined resolution
      markers_list <- list()
      for(res in unique(data_combo_seurat@meta.data[[resolution_user]])){
        
        # Find samples that have at least 3 spots for this cluster
        sample_names <- rownames(spot_by_group_num)[which(spot_by_group_num[,res] >= 3)]
        barcode_names <- rownames(data_combo_seurat@meta.data)[which(data_combo_seurat@meta.data$user.Sample_Name %in% sample_names)]
        
        # Perform DEA for each sample relative to each other sample
        dea_list <- Seurat::FindConservedMarkers(subset(data_combo_seurat, cells=barcode_names), 
                                                 ident.1 = res, grouping.var = "user.Sample_Name", 
                                                 assay = params$assay, test.use = params$test_use,
                                                 min.pct = params$pct_min, logfc.threshold = params$logfc_threshold)
        
        # Add cluster name to the result of DEA
        markers_list[[res]] <- tibble::as_tibble(dea_list, rownames="Symbol") %>%
          dplyr::mutate(cluster = res)
        
      }
      
      # Create a single tibble with all clusters
      # *_p_val *_p_val_adj	*_avg_log2FC	*_pct.1	*_pct.2
      markers_con_df <- dplyr::bind_rows(markers_list) %>%
        dplyr::relocate(cluster)
      
      # Find specific column groups
      cols_pval <- grep("p_val_adj", colnames(markers_con_df))
      cols_avfc <- grep("avg_log2FC", colnames(markers_con_df))
      cols_pct1 <- grep("pct.1", colnames(markers_con_df))
      cols_pct2 <- grep("pct.2", colnames(markers_con_df))
      
      # Modify pct columns to be percentages
      markers_con_df[, c(cols_pct1, cols_pct2)] <- round(markers_con_df[, c(cols_pct1, cols_pct2)] * 100, 1)
      # Modify adjusted p-value columns to be -log10 adjusted p-values
      markers_con_df[, cols_pval] <- round(-log10(markers_con_df[, cols_pval]), 1)
      # Modify log2 FC columns to be rounded
      markers_con_df[, cols_avfc] <- round(markers_con_df[, cols_avfc], 3)
      
      # Select only relevant columns
      markers_con_df <- markers_con_df[, c(c(1, 2), cols_avfc, cols_pval, cols_pct1, cols_pct2)]
      
      # Write conserved markers to file
      filename <- paste0(output_figs, "/", filename_prefix, "_markers_conserved.txt")
      readr::write_delim(markers_con_df, filename, delim="\t")
      
      
      # Find specific column groups
      cols_pval <- grep("p_val_adj", colnames(markers_con_df))
      cols_avfc <- grep("avg_log2FC", colnames(markers_con_df))
      
      # Find genes with same directionality
      dir_pass_cols <- which(rowSums(markers_con_df[, cols_avfc] > 0) == length(cols_avfc) |
                               rowSums(markers_con_df[, cols_avfc] < 0) == length(cols_avfc))
      # Find genes that have FC higher than the threshold in all samples
      fc_pass_cols <- which(rowSums(abs(markers_con_df[, cols_avfc]) >= params$sign_fc_threshold) == length(cols_avfc))
      # Find genes that are significant in all samples
      pv_pass_cols <- which(rowSums(markers_con_df[, cols_pval] >= params$sign_pval_adj_neglog10_threshold) == length(cols_pval))
      
      # Assign direct to genes
      markers_con_df[which(rowSums(markers_con_df[, cols_avfc] > 0) == length(cols_avfc)), "direction"] <- "UP"
      markers_con_df[which(rowSums(markers_con_df[, cols_avfc] < 0) == length(cols_avfc)), "direction"] <- "DN"
      
      # Select rows with significantly upregulated genes
      # Alternative implementation: Reduce(intersect, list(a,b,c))
      markers_con_filt_df <- markers_con_df[intersect(fc_pass_cols, intersect(pv_pass_cols, dir_pass_cols)), ] %>%
        dplyr::relocate(direction, cluster, Symbol)
      
      # Write filtered conserved markers to file
      filename <- paste0(output_figs, "/", filename_prefix, "_markers_conserved_filt.txt")
      readr::write_delim(markers_con_filt_df, filename, delim="\t")
      
      
      # Gather stats for DE markers for each clustering
      stat_df <- markers_con_filt_df %>%
        dplyr::group_by(cluster) %>% 
        dplyr::summarise(UP = sum(direction == "UP"),
                         DN = sum(direction == "DN"),
                         COUNT = dplyr::n_distinct(Symbol)) %>% 
        dplyr::mutate(Resolution = var,
                      Sample_Name = cohort_name) %>%
        dplyr::select(Sample_Name, Resolution, cluster, COUNT, UP, DN)
      
      # Write marker stats to a file
      filename <- paste0(output_figs, "/", filename_prefix, "_markers_conserved_stat.txt")
      readr::write_delim(stat_df, delim="\t", filename)
      
    }    
  }
}
