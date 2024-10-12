# Author: Anna Lyubetskaya. Date: 20-11-20


##_ PARAMETERS _##



# Variables: data_seurat, signature_wide_df, output_folders, run_min

# Author: Anna Lyubetskaya. Date: 20-04-22
# Useful link: https://drsimonj.svbtle.com/exploring-correlations-in-r-with-corrr


## ENVIRONMENT ----


source("code/utils/utils_in_out.R")
source("code/utils/utils_signatures.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_ggplot.R")

source("code/R/Signatures/ST_analysis_signatures_group_utils.R")
#source("code/R/Utils/utils_10X_vis.R")
#source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_dea.R")


## PARAMETERS ----


# Path to processed Seurat data
sample_name <- "P-20200612-0001_ST_FFPE_HumanPanc_A"
#sample_name <- "P-20200625-0002_ST_FF_HumanPanc_D"

# Two signature of interest, one is the focus, the other is a condition
sig <- "CT_pdac_collisson_classical"


## LOAD SEURAT DATA ----


# Location of pre-processed data
input_path <- "C:/Users/lyubetsa/Documents/Projects/10X_ST/Reports/"

# All data - RDS
filename <- paste0(input_path, sample_name, "_all.Sobj.rds")
data_seurat <- readRDS(filename)

# Output information
output_folder <- paste0("C:/Users/lyubetsa/Documents/Projects/10X_ST/CompareSigGroups/")


## LOAD AND FILTER SIGNATURE DATA ----


# Path to signatures
sig_path <- "C:/Users/lyubetsa/Documents/Data/Signatures/processed_signatures_standard_Aug2020.txt"

# Load and filter signatures to only those genes that are well represented 
# Note: These signatures are filtered down to contain only well represented genes as defined above
signature_list <- read_filter_signatures_my(sig_path, rownames(data_seurat))

# Add "sig." prefix to signature names
names(signature_list) <- paste0("sig.", names(signature_list))

# Find a signatures list for every gene
sig_inverted_list <- invert_list_my(signature_list)


## EXTRACT SIGNATURE SCORE DATA ----


# Wide tibble of signature scores + nCount_Spatial and seurat_clusters
sig_wide_df <- tibble::as_tibble(data_seurat@meta.data, rownames="id") %>%
  dplyr::rename(Coordinate = id) %>%
  dplyr::select(dplyr::all_of(c("Coordinate", "nCount_RNA", "seurat_clusters", sig1, sig2)))


##_ MEASURE SIGNATURE STRENGTH IN NEIGHBORHOODS _##


# Visualize a pair of signatures being compared
filename <- paste0(output_folder, "/spatial_signature_pair")
spatial_visualization_gene_list_my(data_seurat, c(sig1, sig2), filename)

sig_select_df <- sig_wide_df %>%
  dplyr::select(dplyr::all_of(c("Coordinate", sig1, sig2)))

# For each coordinate, find its proximal coordinates and calculate a mean of a signature in that selection of spots
sig_select_df <- find_neighborhood_mean_signal_my(data_seurat, sig_wide_df, sig1, sig_select_df, radius=1) %>%
  dplyr::rename(Sig_area1 = Sig_area)

# For each coordinate, find its proximal coordinates and calculate a mean of a signature in that selection of spots
sig_select_df <- find_neighborhood_mean_signal_my(data_seurat, sig_wide_df, sig2, sig_select_df, radius=1) %>%
  dplyr::rename(Sig_area2 = Sig_area)


##_ BIN SPOTS BY SIGNATURE AND NEIGHBORHOOD STRENGTH _##


# Partition spots into hi-lo quadrants based on scores of two signatures
sig_binned_df <- sig_select_df %>%
  vector_categories_hi_lo_my(sig1, coef=1) %>%  # Bin first signature
  dplyr::rename(Category_sig1 = Category) %>%
  vector_categories_hi_lo_my(sig2, coef=1) %>%  # Bin second signature
  dplyr::rename(Category_sig2 = Category) %>%
  vector_categories_hi_lo_my("Sig_area1", coef=0.5) %>%  # Bin first signature neighborhoods
  dplyr::rename(Category_area1 = Category) %>%
  vector_categories_hi_lo_my("Sig_area2", coef=0.5) %>%  # Bin second signature neighborhoods
  dplyr::rename(Category_area2 = Category)

# Count number of categorical labels
bin_stat_df <- sig_binned_df %>% 
  dplyr::group_by(Category_sig1, Category_sig2, Category_area1, Category_area2) %>% 
  dplyr::summarise(n = dplyr::n_distinct(Coordinate))


##_ DEFINE COMPARISON _##


# Spots high in signature 1, low in signature 2, but high in signature 2 neighborhood
proximal <- sig_binned_df %>%
  dplyr::filter(Category_sig1 == "Hi" & Category_sig2 != "Lo" & 
                  Category_area1 == "Hi" & Category_area2 != "Lo")

# Spots high in signature 1, low in signature 2, and low in signature 2 neighborhood
distal <- sig_binned_df %>%
  dplyr::filter(Category_sig1 == "Hi" & Category_sig2 == "Lo" & 
                  Category_area1 == "Hi" & Category_area2 == "Lo")

comparison_groups <- list("proximal" = proximal$Coordinate,
                          "distal" = distal$Coordinate)


# Create a comparison column in meta data
data_seurat@meta.data[["Comparison"]] <- "other"

# Add groups to meta data
for(name in names(comparison_groups)){
  index <- which(rownames(data_seurat@meta.data) %in% comparison_groups[[name]])
  data_seurat@meta.data[index, "Comparison"] <- name
}

data_seurat@meta.data[["Comparison"]] <- as.factor(data_seurat@meta.data[["Comparison"]])

# Visualize groups to be compared
p <- Seurat::SpatialDimPlot(data_seurat, group.by="Comparison", cols=c("blue", "grey", "red"))

filename <- paste0(output_folder, "/spatial_signature_interest_groups")
write_plot2file_my(p, filename)


##_ COMPARE GROUPS _##


if(nrow(proximal) >= 3 && nrow(distal) >= 3){
  # Create a user-defined set of clusters for DEA and add them to active.ident of seurat object
  #data_seurat <- redefine_active_ident_seurat_my(data_seurat, )
  
  # Define genes to perform DEA on
  gene_dea_list <- rownames(data_seurat)[which(data_seurat@assays$SCT@meta.features$Good_data == TRUE)]
  
  # Perform DEA analysis on selected coordinates
  de_markers_df <- compare_clusters_my(data_seurat, find_all = FALSE, group_by_var="Comparison") %>%
    dplyr::filter(cluster_pair == "distal-proximal")
  
  if(run_min == FALSE){
    
    de_markers_plot_df <- de_markers_df %>% 
      dplyr::filter(p_val_adj_neg_log10 >= 2) %>%
      dplyr::inner_join(sig_inverted_list, by="Symbol") %>%
      dplyr::mutate(Expected = grepl(paste0(sig1, "|", sig2), Sig_Name))
    
    # a <- as.matrix(sort(table(gsub("sig.|brownea_|panc_|kumarn_|siemersn_|prostate_|drokhle_|powlesr_prostate_|", "", unlist(unname(sapply(de_markers_plot_df$Sig_Name, function(x) strsplit(x, ";"))))))))
    
    seurat_values_df <- seurat_expression_to_long_tibble_my(data_seurat, assay="SCT")
    
    spot_filt <- sig_binned_df %>%
      dplyr::filter(Category_sig1 == "Hi") %>%
      dplyr::pull(Coordinate)
    
    gene_filt <- de_markers_plot_df$Symbol
    
    params <- list(cell_value = "SCT_zscore",
                   row_label = "Symbol", 
                   col_label = "Coordinate", 
                   distance = "euclidean",
                   row_annotation = c("InSignature", "Expected"),
                   col_annotation = c("Category_area2"),  # "Category_sig1", "Category_sig2", "Category_area1", 
                   range = c(-2, 0, 2),
                   colors = c("red3", "white", "royalblue4"))
    
    create_heatmap_my(seurat_values_df, params, 
                      row_list=gene_filt, col_list=spot_filt,
                      row_meta_df = de_markers_plot_df, col_meta_df = sig_binned_df, 
                      filename=paste0(output_folder, "heatmap_groups.png"))
  }
}
