# Author: Anna Lyubetskaya. Date: 20-08-06
# Pathology v signature


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_signatures.R")

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_signatures.R")
source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_10X_vis.R")

source("code/utils/utils_ggplot.R")

library(Seurat)


## PARAMETERS ----


# The cohort of interest: regex and name
# HumanPanc_FF_cca_sct; HumanPanc_ROI1_FFPE_cca_sct; HumanPanc_S001890_FFPE_C_Apr21_annotated_pathology; HumanPanc_S001891_FFPE_D_Apr21_annotated_pathology
# Syng_B16_cca_sct, Syng_MC38_cca_sct, Syng_B16-MC38_cca_sct
cohort_name <- "Syng_MC38_cca_sct"

# Gene abundance filters
sct_threshold <- 0
spot_threshold <- 0

# Specific signatures to use
# c("bailey.pdac.adex", "bailey.pdac.immunogenic", "bailey.pdac.progenitor", "bailey.pdac.squamous", "collisson.pdac.classical", "collisson.pdac.exocrine", "collisson.pdac.quasimesenchymal", "moffitt.pdac.basal", "moffitt.pdac.classical")
sig_select <- c("U.State.Hypoxia.metabolism", "U.Immune.Tcell", "U.Colon.CRC_Fibroblast", "U.Stroma.CAF")

# Pathology class to compare to
# Pathology.Epithelium.percent; Pathology.Melanin; Pathology.Necrosis; Pathology.Distance.Tissue.filled
pathology_class <- "Pathology.Necrosis"
# Normalize pathology distances by sample
normalize_dist_by_sample <- TRUE

# Or single genes to plot
# c("TRAF2", "TRAF3", "BIRC2", "BIRC3", "MAP3K14", "CHUK", "RELB", "NFKB2", "PTPN2", "ADAR")
# c("TYRP1", "PMEL", "MLANA")
gene_interest_list <- NULL

# Other parameters to plot
# c("nCount_Spatial", "nFeature_Spatial", "Pathology.Necrosis")
param_list <- c("nCount_Spatial", "nFeature_Spatial", "Pathology.Necrosis")

# Order of clusters to plot
cluster_order <- NULL

# User defined color schema - when full control of vis is desired
# HumanPanc_FF_cca_sct_0.2
#cols <- list("0"="#db5a7c", "1"="#96257d", "2"="#fae573", "3"="#4cb0e1", 
#             "4"="#d3d3d3", "5"="#a0d5b5", "6"="#e5e5aa")
# HumanPanc_ROI1_FFPE_cca_sct_0.2
#cols <- list("0"="#fcb352", "1"="#e5e5aa", "2"="#4cb0e1", "3"="#056db5", 
#             "4"="#96257d", "5"="#aee0ea", "6"="#153c65")
# HumanPanc_S001890_FFPE_C_Apr21_0.3
#cols <- list("0"="#d3d3d3", "1"="#dd8231", "2"="#96257d", "3"="#4cb0e1", 
#             "4"="#a0d5b5", "5"="#aee0ea", "6"="#fcb352")
# HumanPanc_S001891_FFPE_D_Apr21_0.6
#cols <- list("0"="#aee0ea", "1"="#96257d", "2"="#e5e5aa", "3"="#fcb352", 
#             "4"="#dd8231", "5"="#4cb0e1", "6"="#d3d3d3")
# Syng_B16_cca_sct_0.2
#cols <- list("0"="#056db5", "1"="#fae573", "2"="#96257d", "3"="#dd8231", 
#             "4"="#e5e5aa", "5"="#4cb0e1", "6"="#aee0ea")
# Syng_MC38_cca_sct_0.2
cols <- list("0"="#db5b7c", "1"="#dd8231", "2"="#96257d", "3"="#fae573", 
             "4"="#fcb352", "5"="#e090b1")
# Syng_B16-MC38_cca_sct_0.2
#cols <- list("0"="#056db5", "1"="#db5b7c", "2"="#4cb0e1", "3"="#dd8231", 
#             "4"="#d3d3d3", "5"="#96257d")


## PATHS ----


# Input folder
input_path <- paste0("XXXX")

# Output path
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_name, "/")

# Create output folder
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)

# Path to a signature file
#sig_path <- NULL
#sig_path <- "data/import/Signatures/signatures_pdac_t1000_May21.txt"
sig_path <- "data/import/Signatures/signatures_syngeneics_t100_Aug21.txt"


## INGEST DATA ----


# Ingest an RDS Seurat object
data_seurat <- readRDS(input_path)


## CALCULATE SIGNATURE SCORE ----


# Find genes abundant in this sample
gene_list <- seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold,
                                             split_by="user.Sample_Name", assay="SCT", slot="data")

# Calculate signature scores and add them to Seurat meta data
if(!is.null(sig_path)){
  # Load signatures and filter them down to only well represented genes
  signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select,
                                              sig_length_min=1, sig_length_max=1000, ratio_threshold=0)
  
  # Add signature scores to a seurat object
  data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.")
  
}

# Add individual signatures to meta data as though they are signatures
if(!is.null(gene_interest_list)){
  signature_list <- intersect(gene_interest_list, gene_list)
  names(signature_list) <- signature_list
  
  # Add signature scores to a seurat object
  data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.")
  
}

# Include parameters already present in the meta data
if(!is.null(param_list)){
  signature_list <- param_list
  names(signature_list) <- signature_list
}

# Seurat meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

# Re-normalize pathology distance by sample (% of max)
if(normalize_dist_by_sample == TRUE){
  dist_norm_df <- meta_df[c(pathology_class, "user.Sample_Name")] %>%
    dplyr::group_by(user.Sample_Name) %>%
    dplyr::mutate(MaxDist = max(!!rlang::sym(pathology_class)),
                  DISTANCE = !!rlang::sym(pathology_class) / MaxDist * 100)
  
  meta_df[pathology_class] <- dist_norm_df$DISTANCE
}

# Preferred clustering resolution
clust_res <- gsub(".0$", "", data_seurat@misc[["user.Clustering"]])

# Factorize the meta data to plot in a specific order
if(!is.null(cluster_order)){
  meta_order <- match(cluster_order, unique(meta_df[[clust_res]]))
  meta_df[[clust_res]] <- factor(meta_df[[clust_res]], 
                                 levels=unique(meta_df[[clust_res]])[meta_order])
}


## VISUALIZE DATA ----


# Define colors
if(is.null(cols)){
  cols <- define_cols_for_var_my(data_seurat, clust_res)
}

# Seurat renames column names, this step finds new signature names
sig_names_upd <- c(colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))]), param_list)

# Iterate through signatures
for(sig in sig_names_upd){
  
  # Find mean signature and pathology score by cluster
  clust_df <- meta_df %>%
    dplyr::group_by(!!rlang::sym(clust_res)) %>%
    dplyr::summarise(Signature = mean(!!rlang::sym(sig)),
                     Pathology = mean(!!rlang::sym(pathology_class)),
                     Size = round(dplyr::n_distinct(Coordinate)/500)+1)

  # Plot cluster mean signature score v cluster mean pathology score  
  p <- create_scatter_plot_my(clust_df, x_label="Pathology", y_label="Signature", 
                              fill_label=clust_res, filename=NULL, size = "Size",
                              labels=c(pathology_class, sig, cohort_name), stroke=0.5) + 
    ggplot2::scale_fill_manual(values=cols)
  
  filename <- paste0(output_path, "scatter_clust_", cohort_name, "_", sig, "_", pathology_class)
  write_plot2file_my(p, filename)

    
  # Plot spot signature score v spot pathology score  
  p <- create_scatter_plot_my(meta_df, x_label=pathology_class, y_label=sig, 
                              fill_label=clust_res, filename=NULL, size = 0.5,
                              labels=c(pathology_class, sig, cohort_name), stroke=0) + 
    ggplot2::scale_fill_manual(values=cols)
  
  filename <- paste0(output_path, "scatter_spot_", cohort_name, "_", sig, "_", pathology_class)
  write_plot2file_my(p, filename)

  # Plot spot signature score v spot pathology score  
  p <- create_scatter_plot_my(meta_df, x_label=pathology_class, y_label=sig, 
                              fill_label=clust_res, filename=NULL, size = 0.5,
                              labels=c(pathology_class, sig, cohort_name), stroke=0) + 
    ggplot2::scale_fill_manual(values=cols)
  
  filename <- paste0(output_path, "scatter_spot_", cohort_name, "_", sig, "_", pathology_class)
  write_plot2file_my(p, filename)
  
  # Find maximum value of the pathology variable
  path_max <- max(meta_df[[pathology_class]])

  # Calculate deciles of the pathology feature  
  meta_loc_df <- meta_df %>% 
    dplyr::mutate(Part = round(!!rlang::sym(pathology_class) / path_max * 20)) %>%
    dplyr::arrange(Part) %>%
    dplyr::mutate(Part = factor(as.character(Part), levels=as.character(0:20)))
  
  # Plot spot signature score v spot pathology score  
  filename <- paste0(output_path, "box_spot_", cohort_name, "_", sig, "_", pathology_class)
  create_box_plot_my(meta_loc_df, x_label="Part", y_label=sig, 
                     fill_label="orig.ident", filename=filename, labels=c(pathology_class, sig, cohort_name))
  
  # Plot spot signature score v spot pathology score  
  filename <- paste0(output_path, "vio_spot_", cohort_name, "_", sig, "_", pathology_class)
  create_violin_plot_my(meta_loc_df, x_label="Part", y_label=sig, 
                     fill_label="orig.ident", filename=filename, labels=c(pathology_class, sig, cohort_name))
  
}
