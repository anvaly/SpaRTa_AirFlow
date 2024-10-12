# Add various data to the Seurat objects
# 1. Pick the number of PCs to use based on variance explained
# 2. Cluster at a few resolutions and visualize the result via clustree
# 3. Perform FindMarkers for each resolution that produces a unique number of clusters


## SETUP ENVIRONMENT ----


# The library doesn't work without loading because it has dependency on ggraph
library(clustree)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_dea.R")
source("code/R/Utils/utils_10X_matrix.R")

source("code/R/b_Cluster_Biomarkers/Seurat_cluster_utils.R")

## PARAMETERS ----

# Setup parameters
params <- cluster_params_my()

# PCA, minimum variance explained
pca_var_threshold <- 0.2
# List of resolutions to test in clustering
resolution_list <- c(0.8)
# All cluster resolutions tested
cluster_names <- c(paste0("SCT_snn_res.", resolution_list), "CRISPRGeneCategory")
de_reference <- c("NON-TARGETING", "background")

# Minimum cluster size
cluster_num_min <- 3
# Maximun cluster size
cluster_num_max <- 200
# Minimum number of spots in a cluster
cluster_size_min <- 5

# Minimum % in spots
pct_min <- 0.5
# Minimum FC difference to test
logfc_threshold <- 0.25

# FC signficance threshold
sign_fc_threshold <- 0.5
# Signficance p-value threshold
sign_pval_adj_neglog10_threshold <- 2

# Color scheme
col_type <- "jet"

# If no covariates present, default to Wilcoxon for DEA
params[["test_use"]] <- "wilcox"  # wilcox; "LR", "negbinom", "poisson", or "MAST"
params[["assay"]] <- "SCT"  # Spatial or SCT
params[["resolution_list"]] <- 0.2
set.seed(981)

## PATHS ----


# Input folder
input_path <- "XXXX"

# Output folder
output_path <- "XXXX"
output_figs <- "XXXX"

# Seurat RDS files are tagged as follows
data_regex_files <- "_all.rds"

# Create output folders
dir.create(output_figs, showWarnings = FALSE)

# Cell Cycle Genes

cc_genes <- readr::read_csv("data/import/Signatures/mouse_cell_cycle_genes.txt")
s_genes <- dplyr::filter(cc_genes, phase == "S") %>% dplyr::select(geneID) %>% unlist %>% unname
g2m_genes <- dplyr::filter(cc_genes, phase =="G2/M") %>% dplyr::select(geneID) %>% unlist %>% unname
## INGEST DATA ----


# Find all files following the regex
file_list <- dir(input_path, pattern=data_regex_files, full.names=TRUE)

## WRANGLE DATA ----

for(rds_file in file_list){
  
  name <- gsub("_all", "", gsub("^.+/|\\.[^\\.]+$", "", rds_file))
  filename_rds <- paste0(output_path, name, "_annotated.rds")
  
  print(name)
  
  if(!file.exists(filename_rds)){
    
    
    ## Load data ----
    
    data_seurat <- readRDS(rds_file)
    sample_name <- data_seurat@misc$user.Sample_Name
    
    data_seurat@meta.data$CRISPRGeneCategory <- gsub("\\-\\d+$", "", data_seurat@meta.data$CRISPRGeneCategory)
    
    # Create a sample-specific output folder
    output_figs_sample <- paste0(output_figs, sample_name)
    dir.create(output_figs_sample, showWarnings = FALSE)
    filename_rds <- paste0(output_path, sample_name, "_annotated.rds")
    # 
    # ## PCA, UMAP, FindNeighbors, FindClusters ----
    # 
    # # Perform PCA, visualize variance explained by each PC, and report the number of useful PCs
    # filename <- paste0(output_figs_sample, "/pca_barplot_", sample_name, "")
    # seurat_pca <- seurat_pca_my(data_seurat, var_threshold=pca_var_threshold, output_file=filename)
    # 
    # # Perform UMAP and find neighbors on previously PCA Seurat data
    # data_seurat <- seurat_umap_nb_my(seurat_pca$data, num_dimensions=seurat_pca$num_pcs)
    # 
    data_seurat <- Seurat::CellCycleScoring(data_seurat, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)
    # 
    # # dim_plot_my(data_seurat, group.by = "CRISPRGuideCategory", pt.size = 1.8, split.by = "CRISPRGeneCategory", ncolumns = 4)
    # 
    # p1 <-dim_plot_my(data_seurat, group.by = "Phase", pt.size = .8)  + theme(aspect.ratio = 1) + ggtitle("Normal SCT")
    # 
    regressed_data_seurat <- Seurat::SCTransform(object = data_seurat, assay = "RNA", vars.to.regress = c("ribo_percent", "mito_percent", "S.Score", "G2M.Score"), vst.flavor = "v2", verbose = FALSE,
                                                 return.only.var.genes = FALSE)
    # seurat_pca <- seurat_pca_my(regressed_data_seurat, var_threshold=pca_var_threshold, output_file=filename)
    # 
    # # Perform UMAP and find neighbors on previously PCA Seurat data
    # regressed_data_seurat <- seurat_umap_nb_my(seurat_pca$data, num_dimensions=seurat_pca$num_pcs)
    # p2 <- dim_plot_my(regressed_data_seurat, group.by = "Phase", pt.size = .8) + ggtitle("CC regressed out") + theme(aspect.ratio = 1)
    # filename <- paste0(output_figs_sample, "/GeX_UMAP_cell_cycle")
    # write_plot2file_my(patchwork::wrap_plots(p1, p2, ncol = 2), num_col = 2, filename = filename)
    # 
    # regressed_data_seurat_guides <- subset(regressed_data_seurat, subset = CRISPRGuideCategory != "")
    # p3 <- dim_plot_my(regressed_data_seurat_guides,group.by = "CRISPRGuideCategory", pt.size = 1.5, split.by = "CRISPRGuideCategory", ncolumns = 5)
    # filename <- paste0(output_figs_sample, "/GeX_UMAP_by_guide")
    # write_plot2file_my(filename = filename, in_plot = patchwork::wrap_plots(p3), num_col = 5, num_row = 6)
    
    
    ## Unsupervised clustering
    
    # Perform PCA, UMAP, FindNeighbors, FindClusters, and FindMarkers
    regressed_data_seurat <- cluster_analysis_my(regressed_data_seurat, params, sample_name, output_figs_sample)
    
    ## Save annotated Seurat object ----
    
    # Save the full dataset to RDS
    saveRDS(regressed_data_seurat, file = filename_rds)
    
  }
}

anno_filename <- "XXXX"
data_seurat <- readRDS(anno_filename)

Seurat::VlnPlot(data_seurat, features = c("ENSMUSG00000016496"), assay = "SCT", group.by = "CRISPRGuideCategory", slot = "scale.data")
good_cells <- subset(data_seurat, subset = CRISPRGuideCategory != "")
Seurat::VlnPlot(good_cells, features = c("ENSMUSG00000095180"), assay = "SCT", group.by = "CRISPRGuideCategory", slot = "scale.data") +labs(x = "sgRNA Guide") + theme(plot.title = element_blank(), aspect.ratio = .7, legend.position = "none")

regressed_data_seurat_subset <- subset(regressed_data_seurat, subset = CRISPRGeneCategory != "")

dim_plot_my(regressed_data_seurat_subset, group.by = "CRISPRGeneCategory", pt.size = 1.4)
Seurat::VlnPlot(regressed_data_seurat, features = c("ENSMUSG00000074151"), assay = "SCT", group.by = "CRISPRGuideCategory", slot = "scale.data") # nlrc5


feature_plot_my(regressed_data_seurat_subset, var = "ENSMUSG00000095180", max.cutoff = 2)


good_cells <- subset(data_seurat, subset = CRISPRGuideCategory != "")
npcs = 50
data_seurat <- seurat_pca_my(good_cells, npcs = npcs, var_threshold=pca_var_threshold)
data_seurat <- seurat_umap_nb_my(data_seurat$data, num_dimensions=data_seurat$num_pcs)

dim_plot_my(data_seurat, group.by = "CRISPRGeneCategory", split.by = "CRISPRGeneCategory", ncolumns = 4)
