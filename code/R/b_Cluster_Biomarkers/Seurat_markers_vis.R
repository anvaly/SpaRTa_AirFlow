# Author: Anna Lyubetskaya. Date: 21-01-26
# Create heatmaps of specific genes and clusters


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")
source("code/utils/utils_signatures.R")
source("code/utils/utils_heatmap.R")
source("code/R/Utils/utils_10X_matrix.R")
source("code/R/b_Cluster_Biomarkers/Seurat_cluster_utils.R")
source("code/R/c_Signatures/ST_analysis_signatures_corr_utils.R")

library(Seurat)


## PARAMETERS ----


## Parameters to find the right file

# The name of the input sample / cohort
cohort_name <- "PDAC108_path14_5K_harmony"
# Stat method
stat_method <- "wilcox"
# Resolution column and file name
resolution <- NULL  #  "Pathology_Group", "integrated_snn_res.0.9", "SCT_snn_res.0.3"
# Values to use in the heatmap
file_type <- "top10"  # unique_top10 or top10


## Parameters to select genes

# Filter genes by signature annotation: Signature_Name field
filter_sig <- "in-sig"  # "in-sig", "non-sig", "all"
# Select only up-regulated genes
filter_direction <- "UP"  # UP, DN, NULL


## Parameters to plot

# Expression data to extract
assay <- "SCT"
slot <- "scale.data"
cell_value <- paste0(assay, "_", slot)

# Use the following column to annotate clusters
annotation_col <- "Sig_Name" # "Sig_Name", "MSigDB", "KEGG

# Number of most frequent annotation terms to use
annotation_topn <- 1

# Additional spot annotations to add to the heatmap
# c("user.Sample_Name", "Pathology.Group")
row_annotation_user <- c("Pathology.Group")

# PCT threshold for plotting
pct_min <- 0

# Select maximum spot size
max_spot_size <- FALSE

# Clustering approach
cluster_method <- "pearson"


## Plotting cluster order controls


# Set defaults
pt.size <- 0.02
cluster_order <- NULL
cols <- NULL


## Harmony full integration, res 0.1


if(cohort_name == "PDAC108_path14_5K_harmony"){
  
  # Pt size for PCA / UMAP
  pt.size <- 0.02
  
  # User defined cluster order
  cluster_order <- c(3, 6, 4, 1, 8, 0, 5, 2, 7, 10, 11, 9)
  
  # User defined color schema - when full control of vis is desired
  cols <- c("#056DB5", "#153C65", "#4CB0E1", "#9A2626", "#96257D", "#E5E5AA",
            "#A0D5B5", "#FF9F2C", "#DF7126", "#AEE0EA", "#D4EEF5", "#939393")
  
  resolution <- "integrated_snn_res.0.1"
}


## Harmony-defined epi niche integrated by RPCA, res 0.4


if(cohort_name == "PDAC108_path14_harmonyepi_rpca_sct"){
  
  # Pt size for PCA / UMAP
  pt.size <- 0.1
  
  # User defined cluster order
  # cluster_order <- c(12,11,3,10,0,4,6,14,9,1,8,5,2,7,13)
  cluster_order <- c(3, 10, 0, 4, 6, 9)
  
  # User defined color schema - when full control of vis is desired
  # cols <- c("#056DB5","#153C65",
  #           "#4CB0E1","#D08EB3","#F16666","#9A2626","#96257D","#360F2E","#681D59",
  #           "#E5E5AA","#5E5E39","#358E5B","#FF9F2C","#D66100","#939393")
  cols <- c("#4CB0E1","#D08EB3","#F16666","#9A2626","#96257D","#681D59")
  
  resolution <- "integrated_snn_res.0.4"
}


## Harmony-defined stroma niche integrated by RPCA, res 0.2


if(cohort_name == "PDAC108_path14_harmonystr_rpca_sct"){
  
  # Pt size for PCA / UMAP
  pt.size <- 0.1
  
  # User defined cluster order
  cluster_order <- c(11,7,0,3,1,4,8,5,2,6,10,9)
  
  # User defined color schema - when full control of vis is desired
  cols <- c("#056DB5","#9A2626","#E5E5AA","#B5B572","#5E5E39","#184928",
            "#A0D5B5","#358E5B","#FF9F2C","#D66100","#49392E","#939393")
  
  resolution <- "integrated_snn_res.0.2"
}


## Harmony-defined immune niche integrated by RPCA, res 0.2


if(cohort_name == "PDAC108_path14_harmonyimm_rpca_sct"){
  
  # Pt size for PCA / UMAP
  pt.size <- 0.1
  
  # User defined cluster order
  cluster_order <- c(10,16,9,0,1,4,3,7,14,2,11,12,13,8,5,6,15)
  
  # User defined color schema - when full control of vis is desired
  cols <- c("#056DB5","#153C65","#9A2626","#E5E5AA","#B5B572",
            "#5E5E39","#878545","#358E5B","#184928","#FBE99E","#FEDC4B",
            "#FAE573","#FFCE07","#FF9F2C","#D66100","#A45024","#939393")
  
  resolution <- "integrated_snn_res.0.4"
}


names(cols) <- cluster_order

# User defined gene list to plot
# default = NULL
# user_gene_list_file <- "XXXX"
# user_gene_list_file <- "XXXX"
user_gene_list_file <- "XXXX"


## PATHS ----


# Composite file name
filename_composite <- paste0("markers_significant_annotated_", file_type, "_", cohort_name, "_", 
                             gsub(".+_snn_res.", "", resolution), "_", stat_method, ".txt")

# Input folder
input_path_object <- paste0("XXXX")
input_path_table <- paste0("XXXX", filename_composite)


# Output subfolder
if(is.null(user_gene_list_file)){
  output_subfolder <- paste0(file_type, "_", cohort_name, "_", 
                             gsub(".+_snn_res.", "", resolution), "_", stat_method)
} else{
  output_subfolder <- paste0(file_type, "_", cohort_name, "_", 
                             gsub(".+_snn_res.", "", resolution), "_", stat_method, "_", 
                             gsub(".+/|.txt", "", user_gene_list_file))
}

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, output_subfolder, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Read in a user defined gene list
if(!is.null(user_gene_list_file)){
  user_gene_list <- readr::read_delim(user_gene_list_file, delim="\t")$Symbol
}

# Seurat data
data_seurat <- readRDS(input_path_object)

# Subset the cohort by cluster if necessary
if(!is.null(cluster_order)){
  barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data[[resolution]] %in% cluster_order),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# Preferred clustering resolution
if(is.null(resolution)){
  resolution <- data_seurat@misc$user.Clustering
}

# Adjust point size if necessary
if(max_spot_size == TRUE){
  data_seurat@misc$user.pt.size.factor <- max(data_seurat@meta.data$user.pt.size.factor)
}

# Update output name
output_filename <- paste0(cohort_name, "_", stat_method, "_", file_type, "_pct", pct_min, "_", gsub(".+_snn_res.", "", resolution), "_", filter_sig, "_", filter_direction)

# Annotated DE genes
dea_df <- readr::read_delim(input_path_table, delim="\t") %>%
  dplyr::filter(pct_1 >= pct_min)


## WRANGLE DEA DATA ----


# Filter DE genes by signature
if(filter_sig == "in-sig"){
  dea_df <- dea_df %>%
    dplyr::filter(!is.na(Sig_Name))
} else if(filter_sig == "non-sig"){
  dea_df <- dea_df %>%
    dplyr::filter(is.na(Sig_Name))
}

# Filter DE genes by direction
if(!is.null(filter_direction)){
  dea_df <- dea_df %>%
    dplyr::filter(direction == filter_direction)
}

# Analyze dea_df stats
dea_stat_df <- dea_df %>%
  dplyr::select(dplyr::all_of(c("Symbol", "cluster", annotation_col))) %>%
  tidyr::separate_rows(!!rlang::sym(annotation_col), sep=";") %>%
  dplyr::group_by(cluster, !!rlang::sym(annotation_col)) %>%
  dplyr::summarise(Count = dplyr::n_distinct(Symbol),
                   Symbol = sort(unique(paste(Symbol, collapse="; ")))) %>%
  dplyr::arrange(cluster, desc(Count))

# Write signature summary by cluster to file
filename <- paste0(output_path, "/", output_filename, "_summary.txt")
readr::write_delim(dea_stat_df, filename, delim="\t")

# Select most abundant annotation
dea_stat_df <- dea_stat_df %>%
  dplyr::slice_max(order_by=Count, n=annotation_topn)

# Select genes of interest
gene_list <- dea_stat_df %>%
  tidyr::separate_rows(Symbol, sep="; ") %>%
  dplyr::pull(Symbol) %>%
  unique()

# Refine gene annotation based on abundant categories
gene_stat_df <- dea_stat_df %>%
  tidyr::separate_rows(Symbol, sep=";") %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(cluster = as.character(cluster)) %>%
  dplyr::summarise(Annotation = paste(sort(unique(unlist(strsplit(paste(!!rlang::sym(annotation_col), collapse=";"), "\\;")))), collapse=";"))

# Make sure that column names are compatible with meta data
colnames(gene_stat_df) <- gsub("cluster", resolution, colnames(gene_stat_df))


## WRANGLE SEURAT DATA ----


# Extract SCT expression data as a wide tibble
data_wide_df <- tibble::as_tibble(t(as.matrix(Seurat::GetAssayData(data_seurat, assay=assay, slot=slot))), 
                                  rownames="Coordinate")

# This is due to a number trim I have in marker filtering step
gene_list <- intersect(gene_list, colnames(data_wide_df))

# Extract SCT expression data as a long tibble; calculate z-scores
data_df <- seurat_expression_to_long_tibble_my(data_seurat, assay=assay, slot=slot, gene_list=gene_list)


## WRANGLE SEURAT META DATA ----


# Extract coordinate meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")
meta_df[[resolution]] <- as.character(meta_df[[resolution]])

# Add signature annotations to each cluster
meta_df <- meta_df %>%
  dplyr::left_join(gene_stat_df, by=resolution) %>%
  dplyr::mutate(Annotated_Cluster = paste(!!rlang::sym(resolution), Annotation))


# Perform clustering for the dot plot or not
if(is.null(cluster_order)){
  cluster.idents <- TRUE
} else{
  cluster.idents <- FALSE
}

# If user didn't define the cluster order, just use the default cluster order
if(is.null(cluster_order)){
  cluster_order <- unique(sort(meta_df[[resolution]]))
}

# Factorize the meta data to plot in a specific order
meta_order <- match(cluster_order, unique(meta_df[[resolution]]))
meta_df[[resolution]] <- factor(meta_df[[resolution]], levels=cluster_order)
meta_df[["Annotated_Cluster"]] <- factor(meta_df[["Annotated_Cluster"]], 
                                         levels=unique(meta_df[["Annotated_Cluster"]])[meta_order])

# Add meta data back to the Seurat object for Seurat plotting
data_seurat@meta.data <- meta_df %>%
  tibble::column_to_rownames("Coordinate")


## PLOT CLUSTER UMAP ----


# Create a PCA, UMAP, and spatial plots
output_name <- paste0(output_path, "/cluster_", output_filename)
Seurat_pca_umap_spatial_my(data_seurat, resolution, output_name, cols=cols, pt.size=pt.size)


## PLOT DOTPLOT ----


# Swap out the gene_list for user-defined if available
if(!is.null(user_gene_list)){
  gene_list <- user_gene_list
}

# Create a dot plot
p <- Seurat::DotPlot(data_seurat, features = gene_list, group.by = resolution,
                     assay = assay, col.min = 0, col.max = 1, cluster.idents=cluster.idents) +
  Seurat::RotatedAxis() +
  ggplot2::scale_y_discrete(limits=rev)

# Define the plot width
num_col <- round(length(gene_list)/50)+1
if(num_col < 1){
  num_col <- 2
}

# Write to file
filename <- paste0(output_path, "/", output_filename, "_dot")
write_plot2file_my(p, filename, num_row=1, num_col=num_col)


## PLOT HEATMAP ----


# Create a heatmap of DE biomarkers
p <- Seurat::DoHeatmap(object = data_seurat[gene_list, sample(colnames(data_seurat), size=30000, replace=F)],
                       features = gene_list,
                       group.by = resolution,
                       slot = slot,
                       assay = assay,
                       group.bar = TRUE,
                       group.colors = cols,
                       disp.min = -2,
                       disp.max = 2
) + ggplot2::theme(text = element_text(size=20))

# Write to file
filename <- paste0(output_path, "/", output_filename, "_hm")
write_plot2file_my(p, filename, 
                   num_row=ceiling(length(gene_list)/50)+1, 
                   num_col=ceiling(length(cluster_order)/10))


## PAIRWISE CORRELATION ANALYSIS ----


# Subset the wide tibble to genes of interest
data_wide_df <- data_wide_df %>%
  dplyr::select(c("Coordinate", dplyr::all_of(gene_list)))

# Calculate gene pairwise correlations
corrr_wide_df <- corrr::correlate(data_wide_df %>%
                                    tibble::column_to_rownames("Coordinate"), quiet = FALSE) %>%
  dplyr::mutate_if(is.numeric, round, 3)

# Create a long tibble of correlations
corr_df <- corrr_wide_df %>%
  df_wide2long_my(key="Gene2", val="R2") %>%
  dplyr::rename(Gene1 = term) %>%
  tidyr::replace_na(list(R2 = 0)) %>%
  dplyr::arrange(Gene1, Gene2)

# Write gene correlation matrix to file
filename <- paste0(output_path, "/", output_filename, "_corr.txt")
readr::write_delim(corr_df, filename, delim = "\t", append=FALSE, col_names=TRUE, ".txt")


# Create a correlation heatmap of gene correlations
filename <- paste0(output_path, "/", output_filename, "_corr_hm")
p <- correlation_plot_my(corrr_wide_df, scale=c(-round(max(abs(corr_df$R2)), 1),
                                                round(max(abs(corr_df$R2)), 1),
                                                0), filename=filename)
