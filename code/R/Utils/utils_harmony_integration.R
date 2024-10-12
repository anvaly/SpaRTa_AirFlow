#Harmony Util Functions

harmony_cluster_analysis_my <- function(data_seurat, params, output_path, plot_together = "resolution"){
  
  resolution <- params[["resolution_list"]]
  seed <- params[["seed"]]
  sample_names <- data_seurat@meta.data$user.Sample_Name %>% unique()
  block_names <- data_seurat@meta.data$user.Block_ID
  ## Find Clusters
  for(res in resolution){
    data_seurat <- Seurat::FindClusters(data_seurat, resolution=res, verbose = FALSE, random.seed = seed)
    cat(paste0("Clustering finished for resolution ", res, "\n"))
  }
  
  # Name of all cluster resolutions tested
  params[["cluster_names"]] <- names(data_seurat@meta.data)[grep("snn_res", names(data_seurat@meta.data))]
  prefix <- unique(gsub("\\..+", ".", params[["cluster_names"]]))
  
  # Extract meta data
  meta_data <- data_seurat@meta.data %>%
    dplyr::select(-seurat_clusters)
  
  # Update Seurat meta data
  data_seurat@meta.data <- meta_data
  
 if(plot_together == "resolution") {
   num_samples <- length(unique(data_seurat@meta.data$user.Sample_Name))
   num_columns <- ceiling(sqrt(num_samples))
     for(var in params[["cluster_names"]]){
       fig_outputs <- paste0(output_path, "res_", var)
       harmony_cohort_plot(data_seurat, var=var, split = "user.Sample_Name", output_name=fig_outputs, num_columns = num_columns, image = NULL)
     }
 }else{
   for(section in sample_names){
     name <- section
     sample_output_dir <- paste0(output_path, "/", name)
     dir.create(sample_output_dir, showWarnings = FALSE)
     sample_subset <- subset(data_seurat, user.Sample_Name == section)
     sample_subset
     for(var in params[["cluster_names"]]){
       filename_prefix <- paste0(name, "_", var)
       harmony_plot_spatial_clusters_my(sample_subset, var, paste0(sample_output_dir, "/", filename_prefix), col_names=NULL, cols=NULL,
                                        title=filename_prefix, image = name)
     }
   }
 }
 
  
  return(data_seurat)
}
