# Author: Andy Kavran 23 February 2022
# Combining a lot of plots (H&E, features, clusters, pathology, signatures) into 
# one figure.

## ENVIRONMENT ----

# Allow piping throughout the package
`%>%` <- magrittr::`%>%`

source("code/utils/utils_signatures.R")
source("code/utils/utils_ggplot.R")
source("code/utils/utils_tibble.R")
source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_pathology_halo.R")
source("code/R/Utils/utils_10X_image.R")
source("code/R/Utils/utils_colormap.R")
source("code/R/Utils/utils_10X_signatures.R")

source("code/R/b_Cluster_Biomarkers/Seurat_cluster_utils.R")
## PARAMETERS ----
#exclude_samples <- c()
cohort_name <- "Pt5_ROI3"
sig_select <-c("PDAC.collisson.classical", "PDAC.moffitt.basal",  "PDAC.P19.Fibroblast", "PDAC.U.Exocrine.Acinar", "PDAC.U.Endocrine.Beta",
               "PDAC.U.Immune.Tcell.CD8", "PDAC.U.Immune.Bcell", "PDAC.U.Immune.Macrophage",
               "BMS.Pathway.IFNg", "Syng.U.State.Hypoxia.metabolism","PDAC.moffitt.normalstroma")

# Pathology params
reference_class <- c("Tissue")
remove_class_completely <- c(NULL) 
others <- c("learning islets", "learning tumor", "learning normal ducts", "learning tumor", "learning necrosis",
            "learning stroma", "learning tumor", "learning exocrine", "learning - tumor", "learning - exocrine", 
            "learning - benign", "learning - iselts", "learning - blood vessel", "learning - muscle-like from adjacent tissue (muscle)", 
            "learning - benign (ducts)", "learning - benign (exocrine)", "learning iselts")
remove_class_from_percent <- c("Tissue", "Viable_Tissue", "Blur_Tissue", others)

# Signature Params
# Gene abundance filters
sct_threshold <- .5
spot_threshold <- 5
# Assay and slot to use to calculate signature score
assay <- "SCT"
slot <- "data"


## PATHS ----
#input_seurat_path <-"XXXX"
input_seurat_path <-"XXXX"
input_raw_path <- "XXXX"
# Input file with sample meta data
meta_file <- paste0(input_raw_path, "meta_data.txt")

input_annotations_metadata <- "XXXX"

# Path to a signature file
sig_path <- "data/import/Signatures/signatures_220210.txt"  # "data/import/Signatures/processed_rshiny_20201127.txt"


output_path <- "XXXX"
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----
#sample_name <- "PDAC_Pt5_ROI1"

# Read sample metadata file
cohort_meta_df <- readr::read_delim(file=meta_file, delim="\t")  %>% 
  dplyr::filter(Quality_Status == TRUE)

# Read pathology metadata file
path_meta_df <- readr::read_delim(input_annotations_metadata, delim = "\t") %>% 
  dplyr::filter(Quality_Status == TRUE)

sample_list <- dir(input_seurat_path, pattern=paste0(cohort_name, ".*_annotated.rds"))

for(sample_file in sample_list){
  
  seurat_filename <- paste0(input_seurat_path, sample_file)
  data_seurat <- readRDS(seurat_filename)
  sample_name <- data_seurat@misc$user.Sample_Name
  out_filename <- paste0(output_path, sample_name, "_summary_plots")
  # if(file.exists(paste0(out_filename, ".png")) == TRUE){
  #   cat("Skipping because file exists: ", paste0(out_filename, ".png\n"))
  #   next()
  # }
  
  raw_filename <- cohort_meta_df %>% dplyr::filter(Sample_Name == sample_name) %>% dplyr::pull(FullPath)
  # Ingest raw (unfiltered) data
  
  data_seurat_raw <- read_10X_spatial_folder_my(raw_filename, filename="raw_feature_bc_matrix.h5")
  tissue_position_list <- readr::read_delim(paste0(raw_filename, "/spatial/tissue_positions_list.csv"), 
                                            delim=",", col_names = c("tissue", "id", "row", "col", "imagerow", "imagecol"))
  all_coord <- tissue_position_list$tissue
  real_coord <- data_seurat_raw %>% colnames
  missing_coord <- setdiff(all_coord, real_coord)
  if(length(missing_coord)>0){
    `%!in%` <- Negate(`%in%`)
    tissue_position_list <- tissue_position_list %>% dplyr::filter(tissue %!in% missing_coord)
  }
  # Put a full image into the raw Seurat data object instead of the filtered one
  data_seurat_raw@images$slice1@coordinates <- tissue_position_list %>%
    tibble::column_to_rownames("tissue") %>%
    dplyr::rename(tissue = id)
  
  
   #path_file <- path_meta_df %>% dplyr::filter(Sample_Name == sample_name) %>% dplyr::pull("Pathology_Full_Path") 
  
  # store plots in list
  p_list <- list()
  ## H&E PLOT ----
  p_list[[1]] <- Seurat::SpatialDimPlot(data_seurat, pt.size.factor = 0, crop = TRUE) +ggtitle(sample_name)+ theme(legend.position = "none", title = element_text(size = 8))
  ## FEATURE PLOT ----
  p_list[[2]] <- spatial_feature_plot_my(data_seurat_raw, feature="nFeature_Spatial", min.cutoff="100", max.cutoff="q98",
                                         name=paste0("Features per spot"), crop=TRUE, color = Turbo()) + theme(legend.position = "right")
  ## CLUSTER PLOT ----
  
  # Get cluster membership
  clust_res <- cohort_meta_df %>% dplyr::filter(Sample_Name == sample_name) %>% dplyr::pull(Clustering)
  #clust_res <- data_seurat@misc$user.Clustering
  
  # cluster colors
  n <- nrow(unique(data_seurat[[clust_res]]))
  cols <- Turbo(out.colors = n)
  # Variable levels
  levels_list <- sort(as.character(unlist(unname(unique(data_seurat[[clust_res]])))))
  names(cols) <- levels_list
  
  # Visualize cluster distribution
  p_list[[3]] <- spatial_dim_plot_my(data_seurat, group.by="SCT_snn_res.0.4", col=cols) + theme(legend.position = "right")
  
  ## PATHOLOGY PLOT ----
  # Seurat image object
  # image_structure <- data_seurat@images$slice1
  #
  # # Seurat spot center coordinates in full resolution: Y = imagerow, X = imagecol
  # X <- image_structure@coordinates[["imagecol"]]
  # Y <- image_structure@coordinates[["imagerow"]]
  # # Spot diameter at full resolution
  # spot_radius <- image_structure@scale.factors$spot_diameter_fullres / 2
  #
  # # Load the annotations as a list of owin-s
  # polys <- halo_annot_to_poly_my(path_file)
  # pathology_classes_init <- sort(names(polys))
  # # Remove certain classes from consideration entirely: legacy classes or areas of tissue we want to filter out
  # pathology_classes <- setdiff(pathology_classes_init, c(reference_class, remove_class_completely, remove_class_from_percent))
  #
  # # Subset pathology input to classes of interest for subsequent plotting
  # polys <- polys[c("Tissue", pathology_classes)]
  #
  # # Define colors for pathology plot
  # num_classes <- length(pathology_classes)
  # cols <- RColorBrewer::brewer.pal(num_classes, name = "Accent")
  # names(cols) <- sort(pathology_classes)
  # cols[["Tissue"]] <- "black"
  # cols[["None"]] <- "white"
  #
  # p_list[[4]] <- plot_pathology_ingestion_check_my(polys, X, Y, spot_radius, filename = NULL, main_class="Tissue",
  #                                                title=NULL, cols=cols, class_list=pathology_classes)

  path_major_group <- data_seurat@meta.data$Pathology.Group
  # Define colors for pathology plot
  num_classes <- path_major_group %>% unique %>% length
  cols <- Turbo(out.colors = num_classes)
  p_list[[4]] <- spatial_dim_plot_my(data_seurat, group.by = "Pathology.Group", cols = cols)
  
  # p_list[[4]] <- ggplot() + theme_void()

  ## SIGNATURE PLOTS ----
  # Find genes abundant in this sample
  gene_list <- seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold,
                                               assay=assay, slot=slot)
  
  # Load signatures and filter them down to only well represented genes
  signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select,
                                              sig_length_min=1, sig_length_max=1000, ratio_threshold=0)
  # Add signature scores to a seurat object
  data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.", assay=assay)
  
  # Seurat renames column names, this step finds new signature names
  sig_names_upd <- colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))])
  
  # Set a min and a max value to a Seurat object for a spatial plot
  data_seurat <- set_spatial_min_max_my(data_seurat, sig_names_upd)
  sig_names_upd <- colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))])
  sig_num <- sig_names_upd %>% length
  
  for(ii in (1:sig_num)){ 
    # there are 4 plots before the first plot here.
    p_list[[ii+4]] <- spatial_feature_plot_my(data_seurat, feature=sig_names_upd[ii], name=sig_names_upd[ii],min.cutoff = 0,
                                              max.cutoff = "q99", slot="data", crop = TRUE) + theme(legend.position = "right")
  }
  
  ## WRITE PLOTS ----
  out_plot <- patchwork::wrap_plots(p_list, ncol = 4)
  #out_plot <- ggpubr::ggarrange(plotlist = p_list, nrow = 5)
  #out_filename <- paste0(output_path, sample_name, "_summary_plots")
  write_plot2file_my(in_plot = out_plot, filename = out_filename, width = 8, height = 11.5)
}

