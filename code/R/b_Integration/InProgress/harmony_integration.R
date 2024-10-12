# Andy Kavran 25 July 2022
# Integrating Visium datasets using Harmony
#
# Read Merged Seurat Object
# Run PCA
# Run Harmony
#
# QCs: UMAP with different features highlighted

## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

library(Seurat)
# install.packages(harmony)
library(harmony)
library(clustree)

source("code/R/Utils/utils_10X_signatures.R")
source("code/utils/utils_signatures.R")
source("code/R/Utils/utils_10X_in_out.R")
source("code/R/b_Cluster_Biomarkers/Seurat_cluster_utils.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_10X_cohort.R")
source("code/R/Utils/utils_harmony_integration.R")
source("code/R/Utils/utils_10X_vis.R")
source("code/utils/utils_ggplot.R")

## Define paths ----
seurat_path <- c("XXXX")

output_path <- "XXXX"
dir.create(output_path, showWarnings = FALSE)

output_seurat <- "XXXX"
dir.create(output_seurat, showWarnings = FALSE)

save_integrated_rds = TRUE # save seurat_object after running harmony
save_integrated_name <- paste0(output_seurat, "pdac_cohort_spotclean_harmony_integrated.rds")

## PARAMETERS ----
cohort_regex <- "PDAC"
n_var_features <- 7000
nclust_harmony <- 100 ## number of clusters harmony uses as anchors. default nclust = 100
theta <- 2 ## cluster diversity enforcement. larger values encourage more diverse clusters. default theta = 2
lambda <- 1 ## ridge regresion penalty. Smaller values result in more aggressive correction. default lambda = 1
batch_vars = "user.Sample_Name"

paramter_string <- paste0(n_var_features, "_features_", nclust_harmony, "_clust_", theta, "_theta_", lambda, "_lambda/")
output_path <- paste0(output_path,paramter_string)
dir.create(output_path, showWarnings = FALSE)
dir.create(paste0(output_path, "/spatial/"), showWarnings = FALSE)

seed <- 531

# gene signature parameters
sct_threshold <- .4
spot_threshold <- 5
sig_select <- c("PDAC.collisson.classical","PDAC.moffitt.basal", "PDAC.P19.Fibroblast", "PDAC.U.Exocrine.Acinar", "PDAC.U.Endocrine.Beta", "PDAC.moffitt.normalstroma", "PDAC.U.Immune.Tcell.CD8", "PDAC.U.Immune.Bcell", "PDAC.U.Immune.Macrophage")

# Path to a signature file
sig_path <- "/repos/P02567_TBIO-3021_10X_pilot/data/import/Signatures/signatures_220210.txt"

## Import list of seurat files and create merged object ----

# List of misc fields to keep in meta data
misc_field_list <- c("user.Sample_Name", "user.Group_Name", "user.Tissue",
                     "user.Date", "user.Area", "user.Project_ID", "user.Primary",
                     "user.Protocol", "user.Organism", "user.pt.size.factor",
                     "user.Clustering", "user.Slide", "Pathology.Tissue", "user.Image_Filename",
                     "user.Concentriq_Image_ID", "user.Concentriq_Repo_ID", "user.KRAS", "user.Block_ID")

file_list <-  file_list <- dir(seurat_path, pattern = cohort_regex, full.names = TRUE)

# Ingest a set of RDS Seurat objects
start_time <- proc.time()
seurat_list <- read_rds_list_simple_my(file_list) #full cohort = 12 minutes
cat("Finish importing inidivdual objects\n")
proc.time() - start_time

for(name in names(seurat_list)){
  sample_name <- seurat_list[[name]]@misc$user.Sample_Name #get sample name
  seurat_list[[name]]@active.assay <- "Spatial"
  seurat_list[[name]][["SCT"]] <- NULL
  names(seurat_list[[name]]@images) <- sample_name # rename image slice with sample name
  seurat_list[[name]]@meta.data <- seurat_list[[name]]@meta.data[which(!grepl("SCT_snn_res.", names(seurat_list[[name]]@meta.data)))] # remove old clustering
}

# get_sample_name <- function(seurat_object){
#   sample_name <- seurat_object@misc$user.Sample_Name
#   return(sample_name)
# }
#
# sample_names <- unname(unlist(lapply(seurat_list, get_sample_name)))
# names(seurat_list) <-sample_names

# Move misc parameters to meta data slot
seurat_list <- seurat_misc_to_meta(seurat_list, misc_field_list=misc_field_list)

# Merge individual objects together
start_time <- proc.time()
merged_seurat_data <- seurat_merge_my(seurat_list)
cat("Finished merging objects.\n")
proc.time() - start_time # whole cohort 35 minutes
rm(seurat_list)

# Temp: Correct for NAs in these path annotations
# TODO: Fix this elsewhere
tls_var <- c("Pathology.TLS-Aggregate", "Pathology.TLS-Mature", "Pathology.TLS-Immature", "Pathology.Lymph_Node")
tls_var <- intersect(names(merged_seurat_data@meta.data), tls_var)
for(tls_v in tls_var){
  idx <- is.na(merged_seurat_data@meta.data[[tls_v]])
  merged_seurat_data@meta.data[[tls_v]][idx] <- 0
}


# SCTransform the merged dataset
merged_seurat_data <- Seurat::SCTransform(merged_seurat_data, assay="Spatial", variable.features.n = n_var_features,
                                          return.only.var.genes = TRUE, verbose = FALSE, conserve.memory = FALSE)
cat("SCTransform complete.\n")

merged_seurat_data <- FindVariableFeatures(merged_seurat_data, nfeatures = n_var_features)
## Run PCA and Harmony ----
merged_seurat_data <- Seurat::RunPCA(merged_seurat_data, seed.use = seed)
  
merged_seurat_data <- RunHarmony(merged_seurat_data, group.by.vars = batch_vars, reduction = "pca",
                                 assay.use = "SCT", theta = theta, nclust = nclust_harmony, lambda = lambda,
                                 reduction.save = "harmony")
cat("Harmony Integration complete.\n")

## Run UMAP on new Harmony Embedding ----
merged_seurat_data <- RunUMAP(merged_seurat_data, reduction = "harmony", dims = 1:50)
merged_seurat_data <- FindNeighbors(merged_seurat_data, reduction = "harmony", dims = 1:50, force.recalc = TRUE)
cat("UMAP complete.\n")
## Plot different variables on integrated UMAP ----
plot_group_vars <- c("user.Sample_Name", "user.Block_ID", "Pathology.Group")
for(group in plot_group_vars){
  dim_plot_my(merged_seurat_data, reduction = "umap", group.by = group) + ggplot2::theme(aspect.ratio = 1)
  ggplot2::ggsave(paste0(output_path, "pdac_integrated_umap_harmony_", group, ".png"), height = 10, width = 15)
}

dim_plot_my(merged_seurat_data, reduction = "umap", group.by = "user.Block_ID", split.by = "user.Block_ID", ncolumns = 5) + ggplot2::theme(aspect.ratio = 1)
ggplot2::ggsave(paste0(output_path, "pdac_integrated_umap_harmony_block_facet.png"),  height = 10, width = 20)

dim_plot_my(merged_seurat_data, reduction = "umap", group.by = "Pathology.Group", split.by = "Pathology.Group", ncolumns = 5) + ggplot2::theme(aspect.ratio = 1)
ggplot2::ggsave(paste0(output_path, "pdac_integrated_umap_harmony_pathology_facet.png"), height = 10, width = 20)

dim_plot_my(merged_seurat_data, reduction = "umap", group.by = "Pathology.Group", split.by = "user.Block_ID", ncolumns = 5) + ggplot2::theme(aspect.ratio = 1)
ggplot2::ggsave(paste0(output_path, "pdac_integrated_umap_harmony_pathology_facet_blocks.png"),  height = 10, width = 20)

path_group_vars <- names(merged_seurat_data@meta.data)[grepl("Pathology.*.*",names(merged_seurat_data@meta.data))]
path_group_vars <- path_group_vars[!grepl("percent|Total|Focus|Group", path_group_vars)]
for(path_annotation in path_group_vars){
  feature_plot_my(merged_seurat_data, var = path_annotation) + ggplot2::theme(aspect.ratio = 1)
  ggplot2::ggsave(paste0(output_path, "integrated_umap_harmony", path_annotation, ".png"),  height = 10, width = 20)

}

## Plot Signatures ----
# TODO: Get this function to work on a merged object
# gene_list <- seurat_select_abundant_genes_my(merged_seurat_data, sct_threshold=sct_threshold, spot_threshold=spot_threshold)
gene_list <- merged_seurat_data@assays$SCT@var.features
# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select,
                                            sig_length_max=1000, ratio_threshold=0)
# Add signature scores to a seurat object
merged_seurat_data <- add_signature_scores_my(merged_seurat_data, signature_list, prefix="sig.", assay='SCT')

# Seurat renames column names, this step finds new signature names
sig_names_upd <- colnames(merged_seurat_data@meta.data[grep("sig.", colnames(merged_seurat_data@meta.data))])

for(i in c(1:length(sig_names_upd))){
  feature_plot_my(merged_seurat_data, var = sig_names_upd[i], min.cutoff = "q10", max.cutoff = "q90") %>% print()
  ggplot2::ggsave(paste0(output_path, sig_names_upd[i], "_integrated.png"), width = 7, height = 7)
}

# merged_seurat_data <- readRDS(file = paste0(output_path, save_integrated_name, ".rds"))
## Clustering ----
# Set parameters for the clustering analysis
params <- cluster_params_my()

params[["latent_vars"]] <- c("user.Sample_Name")
params[["resolution_list"]] <- seq(.1, 1.3, .1)
params[["seed"]] <- seed
merged_seurat_data <- harmony_cluster_analysis_my(merged_seurat_data, params, paste0(output_path, "spatial/"))

## Write Integrated Object ----
if(save_integrated_rds == TRUE){
  saveRDS(merged_seurat_data, file = save_integrated_name, compress = FALSE)
  cat(paste0("File Created:", save_integrated_name, "\n"))
}
#### 

# Visualize clustering behavior with clustree

meta_data <- merged_seurat_data@meta.data[,grepl("SCT_snn_res.", names(merged_seurat_data@meta.data))]
filename <- paste0(output_path,"clustree.png")
p4 <- clustree::clustree(meta_data, prefix="SCT_snn_res.")
ggplot2::ggsave(filename, plot = p4, width = 10, height = 8)
cat(paste0("File Created: ", filename, "\n"))

## Other Plots
# merged_seurat_data <- readRDS(save_integrated_name)

### All H&E Plots 
# p4 <- Seurat::SpatialDimPlot(merged_seurat_data, pt.size = 0, combine = FALSE)
# ggsave(plot = patchwork::wrap_plots(p4, nrow = 9, ncol = 10, guides = "collect"), file = paste0(output_path, "H&Es.png"), width = 30, height = 30)

## Clusters on combined UMAP
cluster_vars <- names(merged_seurat_data@meta.data)[grep("snn_res", names(merged_seurat_data@meta.data))]
for(var in cluster_vars){
  p <- dim_plot_my(merged_seurat_data, group.by = var, reduction = "umap") + ggplot2::theme(aspect.ratio = 1, legend.position = "none") 
  write_plot2file_my(p, filename = paste0(output_path, var, "_umap"), width =8, height = 8)
}

