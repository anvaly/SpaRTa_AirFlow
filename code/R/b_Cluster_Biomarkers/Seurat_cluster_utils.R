# Author: Anna Lyubetskaya. Date: 21-02-11
# Perform PCA, UMAP, FindNeighbors, FindClusters, and FindMarkers


source("code/R/Utils/utils_10X_dea.R")
source("code/R/Utils/utils_10X_vis.R")


if(!"clustree" %in% rownames(installed.packages())){
  install.packages("clustree", repos=c('https://pm.rdcloud.bms.com/prod-cran/2024-01-18','https://pm.rdcloud.bms.com/prod-cran/2022-11-02', 'https://pm.rdcloud.bms.com/bms-cg-biogit-bran/latest'))
}

# The library doesn't work without loading because it has dependency on ggraph
library(clustree)


cluster_params_my <- function(){
  ## Set parameters for the clustering analysis
  
  params <- list()
  
  # PCA, minimum variance explained
  params[["pca_var_threshold"]] <- 0.1
  # List of resolutions to test in clustering
  params[["resolution_list"]] <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5)
  
  # Latent variables to regress out
  params[["latent_vars"]] <- NULL  # Only works with test_use "LR", "negbinom", "poisson", or "MAST"
  
  # Biomarker method and data source
  params[["test_use"]] <- "MAST"  # wilcox; "LR", "negbinom", "poisson", or "MAST"
  params[["assay"]] <- "Spatial"  # Spatial or SCT
  
  # Minimum % in spots
  params[["pct_min"]] <- 0.1
  # Minimum FC difference to test
  params[["logfc_threshold"]] <- 0.1
  
  # FC significance threshold
  params[["sign_fc_threshold"]] <- 0.25
  # Significance p-value threshold
  params[["sign_pval_adj_neglog10_threshold"]] <- 2
  
  # Color scheme
  params[["col_type"]] <- "jet"
  
  
  return(params)
}


cluster_analysis_my <- function(data_seurat, params, name, output_path){
  ## Perform PCA, UMAP, FindNeighbors, FindClusters
  
  gc()
  
  
  ## PCA, UMAP, FindNeighbors, FindClusters ----
  
  # Perform PCA, visualize variance explained by each PC, and report the number of useful PCs
  filename <- paste0(output_path, "/pca_barplot_", name, "")
  seurat_pca <- seurat_pca_my(data_seurat, var_threshold=params[["pca_var_threshold"]], output_file=filename)
  
  # Perform UMAP and find neighbors on previously PCA Seurat data
  data_seurat <- seurat_umap_nb_my(seurat_pca$data, num_dimensions=seurat_pca$num_pcs)
  
  # Cycle over resolution values
  for(resolution in params[["resolution_list"]]){
    # Cluster Seurat data
    data_seurat <- seurat_cluster_my(data_seurat, resolution=resolution)
  }
  
  # All cluster resolutions tested
  params[["cluster_names"]] <- names(data_seurat@meta.data)[grep("snn_res", names(data_seurat@meta.data))]
  prefix <- unique(gsub("\\..+", ".", params[["cluster_names"]]))
  
  
  ## Analyze clusters using clustree ----
  
  # Extract meta data
  meta_data <- data_seurat@meta.data %>%
    dplyr::select(-seurat_clusters)
  
  # Update Seurat meta data
  data_seurat@meta.data <- meta_data
  
  # Visualize clustering behavior
  if(length(params[["resolution_list"]]) >= 2){
    p <- clustree::clustree(meta_data, prefix=prefix)
    
    filename <- paste0(output_path, "/clustree_", name, "")
    write_plot2file_my(p, filename, num_row=3.5, num_col=2)
  }
  
  
  ## Visualization for every resolution ----
  
  
  # Reference gene list for gene over-representation analysis
  filename <- paste0(output_path, "/", name, "_ref_gene_list.txt")
  write("Symbol", file=filename)
  write(paste(sort(rownames(data_seurat@assays$SCT)), collapse="\n"), file=filename, append=TRUE)
  
  for(var in params[["cluster_names"]]){
    
    # Short resolution for file names
    var_name <- gsub("SCT_snn_res.|integrated_snn_res.", "", var)
    # Correct resolution "1" to "1.0" for file order
    if(var_name == "1"){
      var_name <- "1.0"
    }
    
    # File names start with sample name, than resolution, than specifix suffix
    filename_prefix <- paste0(name, "_", var_name)
    
    # Table of cluster names and number of spots in each
    cluster_table <- table(droplevels.data.frame(data_seurat[[var]]))
    print(cluster_table)
    
    
    ## Visualize unique cluster sets ----
    
    if(length(names(cluster_table) <= 50)){
      # Visualize PCA, UMAP, and tissue cluster distribution
      Seurat_pca_umap_spatial_my(data_seurat, var, paste0(output_path, "/", filename_prefix), col_names=NULL, cols=NULL, 
                                 title=filename_prefix)
    }
    
  }
  
  
  return(data_seurat)
}


cluster_analysis_markers_my <- function(data_seurat, params, name, output_path){
  ## After a clustering resolution is settled upon, use this to find markers of each cluster
  ## Written by Andy Kavran 12/2/2021
  
  gc()
  
  # Get preferred clustering resolution from meta_data
  cluster_res <- data_seurat@misc$user.Pipeline.Clustering
  
  # Short resolution for file names
  var_name <- gsub("SCT_snn_res.|integrated_snn_res.", "", cluster_res)
  
  # Correct resolution "1" to "1.0" for file order
  if(var_name == "1"){
    var_name <- "1.0"
  }
  
  # File names start with sample name, than resolution, than specifix suffix
  filename_prefix <- paste0(name, "_", var_name)
  
  # Table of cluster names and number of spots in each
  cluster_table <- table(droplevels.data.frame(data_seurat[[cluster_res]]))
  print(cluster_table)
  
  # Put the clustering of interest into seurat clusters
  data_seurat <- Seurat::SetIdent(data_seurat, value=cluster_res)
  
  # Find markers of each cluster against the rest
  markers_df <- seurat_find_markers_my(data_seurat, assay=params[["assay"]], 
                                       find_all=TRUE, group.by=cluster_res, 
                                       min.pct=params[["pct_min"]], 
                                       logfc.threshold=params[["logfc_threshold"]],
                                       test.use=params[["test_use"]], 
                                       latent.vars=params[["latent_vars"]])
  
  ## Find signficant markers ----
  
  markers_filt_df <- marker_analysis_my(markers_df, params, cluster_res, name, output_path, filename_prefix)
  
  
  return(data_seurat)
}


marker_analysis_my <- function(markers_df, params, var, name, output_path, filename_prefix){
  ## Find markers, filter them, and summarize the result
  
  # If at least one marker detected
  if(!is.null(markers_df)){
    
    # Write all marker data to a file
    filename <- paste0(output_path, "/markers_", filename_prefix)
    readr::write_delim(markers_df, paste0(filename, ".txt"), delim="\t")
    
    # Filter marker tibble to significant markers only
    filename <- paste0(output_path, "/markers_significant_", filename_prefix)
    markers_filt_df <- filter_and_vis_markers_my(markers_df, fc_threshold=params[["sign_fc_threshold"]], 
                                                 pval_adj_neglog10_threshold=params[["sign_pval_adj_neglog10_threshold"]], 
                                                 pct_min=params[["pct_min"]], output_file=filename)
    
    # Gather stats for DE markers for each clustering
    stat_df <- markers_filt_df %>%
      dplyr::group_by(cluster) %>% 
      dplyr::summarise(UP = sum(avg_logFC > 0),
                       DN = sum(avg_logFC < 0),
                       COUNT = dplyr::n_distinct(Symbol)) %>% 
      dplyr::mutate(Resolution = var,
                    Sample_Name = name) %>%
      dplyr::select(Sample_Name, Resolution, cluster, COUNT, UP, DN)
    
    # Write marker stats to a file
    filename <- paste0(output_path, "/stat_", filename_prefix, ".txt")
    readr::write_delim(stat_df, delim="\t", filename, append=TRUE, col_names=TRUE)
  }
  
  return(markers_filt_df)
}
