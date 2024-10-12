# Author: Anna Lyubetskaya. Date: 20-06-10


source("code/utils/utils_ggplot.R")


seurat_sct_my <- function(data_seurat, vars_to_regress=c("mito_percent", "ribo_percent")){
  ## Perform SCTransform
  
  # Perform SCTransform normalization
  data_seurat <- Seurat::SCTransform(data_seurat, assay="Spatial", vars.to.regress=vars_to_regress,
                                     return.only.var.genes = FALSE, verbose = FALSE)
  
  return(data_seurat)
}


seurat_pca_my <- function(data_seurat, npcs=50, var_threshold=0.1, output_file=NULL, features = NULL){
  ## Perform PCA using Seurat and understand variance distribution
  ## Provide the number of components that explain more than X data variance
  ## See: https://www.biostars.org/p/423306/
  
  # Perform PCA
  data_seurat <- Seurat::RunPCA(data_seurat, npcs=npcs, verbose = FALSE, features = features)
  
  # Extract PCA information
  pca_data <- data_seurat@reductions$pca
  # Calculate variance from standard deviations
  pca_variance <- (pca_data@stdev)^2
  # Extract total variance of the data
  total_variance <- pca_data@misc$total.variance
  
  # % variance explained by each component
  variance_explained <- round((pca_variance / total_variance)*100, 1)
  
  # Number of PCs explaining above-threshold variance
  num_pcs <- max(10, length(which(variance_explained >= var_threshold)))
  # Round the number of PCs for convenience
  # num_pcs_rounded <- round(num_pcs + 4.1, -1)
  
  if(!is.null(output_file)){
    df <- tibble::tibble(PC = 1:length(variance_explained),
                         Percent_Variance = variance_explained,
                         UMAP_Include = variance_explained >= var_threshold)
    
    create_bar_plot_my(df, x_label="PC", y_label="Percent_Variance", fill_label="UMAP_Include", filename=output_file,
                       labels=c("PC", "Percent Variance Explained", paste("PCs selected =", num_pcs)))
    
    # Jackstraw plot for PCs
    # Seurat::JackStrawPlot(object, dims = 1:5, reduction = "pca", xmax = 0.1, ymax = 0.3)
    
  }
  
  return(list(data=data_seurat,
              num_pcs=num_pcs))
}


seurat_umap_nb_my <- function(data_seurat, num_dimensions=50){
  ## Perform UMAP and find neighbors on previously PCA Seurat data
  
  ## RunUMAP warning
  ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
  ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' (requires package "reticulate") and metric to 'correlation'
  
  ## BuildClusterTree warning
  ## The function seems to fail, specifically with reorder = TRUE, reorder.numeric = T
  ## https://github.com/satijalab/seurat/issues/1252
  
  if("harmony" %in% names(data_seurat@reductions)){
    reduction <- "harmony"
  } else{
    reduction <- "pca"
  }
  
  # Run UMAP
  data_seurat <- Seurat::RunUMAP(data_seurat, dims = 1:num_dimensions, verbose = FALSE, reduction=reduction)  # umap.method = "umap-learn", metric = "correlation"
  
  # Find neighbors
  data_seurat <- Seurat::FindNeighbors(data_seurat, dims = 1:num_dimensions, verbose = FALSE, force.recalc = TRUE, reduction=reduction)
  
  return(data_seurat)
}


seurat_cluster_my <- function(data_seurat, resolution=0.8){
  ## Find clusters on previously PCA-UMAP-FindNeighbors Seurat data
  
  # Define clusters; clusters are in the column named as follows: col_name <- paste0("SCT_snn_res.", resolution)
  data_seurat <- Seurat::FindClusters(data_seurat, resolution=resolution, verbose = FALSE)
  
  cat("Clusters identified:\n", levels(data_seurat@meta.data$seurat_clusters), "\n")
  
  return(data_seurat)
}


seurat_find_markers_my <- function(data_seurat, assay="Spatial", find_all=TRUE, group.by="seurat_clusters", cluster_pairs=NULL,
                                   test.use="wilcox", min.pct=0.1, logfc.threshold=0.25, latent.vars=NULL){
  
  ## Find genes that vary the most between every pair of clusters
  ## Execute the test as either (1) each cluster v other (find_all=TRUE) or (2) for each pair of clusters
  ## FindMarkers discussion: https://www.biostars.org/p/409790/
  
  ## As a result, Seurat outputs: p-value, avg_log2FC, pct.1, pct.2, p_val_adj
  ## pct.1 and pct.2 = The percentage of cells where the gene is detected in the first/second group
  ## Some tests available: "wilcoxon", "MAST"
  ## min.pct and logfc.threshold filter input features and speed up performance
  
  
  # Put the clustering of interest into the Seurat identity
  data_seurat <- Seurat::SetIdent(data_seurat, value=group.by)
  
  # Execute the test as either each cluster v other or for each pair of clusters
  if(find_all == TRUE){
    markers_matrix <- Seurat::FindAllMarkers(object = data_seurat, assay=assay,
                                             min.pct = min.pct, logfc.threshold = logfc.threshold,
                                             test.use = test.use, latent.vars = latent.vars)
    
    if(nrow(markers_matrix) > 0){
      markers_df <- markers_matrix  %>%
        format_markers_output_my() %>%
        dplyr::arrange(cluster, desc(p_val_adj_neg_log10), desc(avg_logFC))
    } else{
      markers_df <- NULL
    }
  } else{
    if(is.null(cluster_pairs)){
      # Generate all combinations of clusters
      cluster_pairs <- combn(levels(data_seurat@active.ident), 2)
      print(cluster_pairs)
    }
    
    # Find markers for each pair of clusters
    markers_list <- list()
    for(i in 1:ncol(cluster_pairs)){
      cluster1 <- cluster_pairs[1, i]
      cluster2 <- cluster_pairs[2, i]
      pair_name <- paste0(cluster1, "-", cluster2)
      
      markers_list[[pair_name]] <- Seurat::FindMarkers(data_seurat, assay=assay,
                                                       ident.1 = cluster1, ident.2 = cluster2, group.by=group.by,
                                                       min.pct = min.pct, logfc.threshold = logfc.threshold,
                                                       test.use = test.use, latent.vars = latent.vars)
      
      markers_list[[pair_name]][["cluster"]] <- pair_name
      
      markers_list[[pair_name]] <- markers_list[[pair_name]] %>%
        format_markers_output_my()
    }
    
    # Create a volcano plot of all pair comparisons
    markers_df <- do.call(rbind, markers_list)
  }
  
  
  return(markers_df)
}


format_markers_output_my <- function(markers_matrix){
  ## Standardize Seurat FindMarkers output
  ## Seurat outputs: p-val, avg_log2FC, pct.1, pct.2, p_val_adj, and cluster in case of FindAllMarkers
  
  # Store formatted variability data as a tibble
  markers_df <- markers_matrix %>%
    tibble::as_tibble(rownames = "Symbol") %>%
    dplyr::mutate(Symbol = gsub("\\.\\d+$", "", Symbol)) %>%
    dplyr::mutate(pct_1 = round(pct.1 * 100),
                  pct_2 = round(pct.2 * 100),
                  avg_logFC = round(avg_log2FC, 3),
                  p_val_adj_neg_log10 = round(-log10(p_val_adj), 1),
                  direction = ifelse(avg_logFC >= 0, "UP", "DN")) %>%
    dplyr::select(c(Symbol, cluster, direction, avg_logFC, p_val_adj_neg_log10, pct_1, pct_2))
  
  return(markers_df)
}


filter_and_vis_markers_my <- function(markers_df, fc_threshold=1, pval_adj_neglog10_threshold=2, pct_min=0.2,
                                      output_file=NULL, make_plot=TRUE){
  ## Select most variable genes and create a volcano plot
  
  # Mark genes as signficant according to thresholds
  markers_df <- markers_df %>%
    dplyr::mutate(IsSignificant = abs(avg_logFC) >= fc_threshold & 
                    p_val_adj_neg_log10 >= pval_adj_neglog10_threshold &
                    pct_1 >= pct_min*100)
  
  # Filter markers down to signficant only
  markers_filt_df <- markers_df %>%
    dplyr::filter(IsSignificant == TRUE) %>%
    dplyr::select(-IsSignificant)
  
  if(!is.null(output_file) && nrow(markers_filt_df) > 0){
    # Write the data to a text file
    readr::write_delim(markers_filt_df, delim="\t", paste0(output_file, ".txt"))
    
    if(make_plot == TRUE){
      # Create a scatter plot illustrating significant gene distribution
      create_scatter_plot_my(markers_df, x_label="avg_logFC", y_label="p_val_adj_neg_log10", 
                             fill_label="IsSignificant", facet_var=c("cluster", "free_y"),
                             shape=19, size=0.75, filename=output_file, labels=NULL, do_fit=NULL)
    }
  }
  
  return(markers_filt_df)
}


find_spatial_markers_my <- function(){
  ## Find spatially variable feature
  
  # Find spatially variable features
  data_seurat <- Seurat::FindSpatiallyVariableFeatures(data_seurat, assay = "SCT", 
                                                       features = Seurat::VariableFeatures(data_seurat)[1:num_var_genes], 
                                                       selection.method = "markvariogram")
  
  top_features <- head(Seurat::SpatiallyVariableFeatures(data_seurat, selection.method = "markvariogram"), 6)
  
  p1 <- Seurat::SpatialFeaturePlot(data_seurat, features = top_features, ncol = 3, alpha = c(0.1, 1))
  write_plot2file_my(p1, paste0(output_folders$DEA, "variable_features_top"))
}
