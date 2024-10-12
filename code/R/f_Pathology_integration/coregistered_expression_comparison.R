# Author: Anna Lyubetskaya. Date: 21-04-22
# For a pair of two experiments, perform spatially-aware comparison of expression patterns to demonstrate reproducibility
# Co-registration data is one of the inputs


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

library(foreach)
library(Seurat)

source("code/utils/utils_tibble.R")
source("code/utils/utils_ggplot.R")
source("code/utils/utils_signatures.R")

source("code/R/Utils/utils_10X_image.R")
source("code/R/Utils/utils_10X_vis.R")


## PARAMETERS ----


# Radius of the spatial neighborhood in spot number away from center
neighborhood_size <- 2
# Seurat assay to use Spatial / SCT
assay <- "SCT"
# Seurat assay to use data / scale.data
slot <- "data"

# Output name prefix
cohort_name <- "HumanPanc_ROI1_v2_vs_v1_SCT_data_N2"

# Pair of coregistered samples to compare
#sample1_name <- "HumanPanc_ROI2_FFPE_D_Dec20_T"
#sample2_name <- "HumanPanc_ROI2_FFPE_B_Apr21"
sample1_name <- "HumanPanc_ROI1_FFPE_B_Dec20_T"
sample2_name <- "HumanPanc_ROI1_FFPE_A_Apr21"


## PATHS ----


# File containing coregistration between sample1 and sample2
# PDAC FFPE ROI2
#coregistration_path <- "XXXX"
# PDAC FFPE ROI1
coregistration_path <- "XXXX"

# Path to a signature file
sig_path <- "data/import/Signatures/signatures_pdac_t1000_May21.txt"

# Input folder
input_path <- "XXXX"

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_name, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Read in first sample data
data_seurat_list <- list(readRDS(paste0(input_path, sample1_name, "_annotated.rds")),
                         readRDS(paste0(input_path, sample2_name, "_annotated.rds")))
names(data_seurat_list) <- c(sample1_name, sample2_name)

# List of genes detectable in both experiments
gene_list <- intersect(rownames(data_seurat_list[[1]]), rownames(data_seurat_list[[2]]))

if(!is.null(sig_path)){
  # Load signatures and filter them down to only well represented genes
  signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_length_min=1, sig_length_max=1000, ratio_threshold=0)
  gene_sig_df <- invert_list_my(signature_list)
  
  # List of genes present in the signatures and both Seurat objects
  gene_list <- intersect(gene_list, gene_sig_df$Symbol)
}

# XML of coordinates matching Sample1 centers to coordinates on Sample2 image
coregistration_xml <- xml2::read_xml(coregistration_path)


## CALCULATE EXPRESSION MEAN ----


neighborhood_expr_mean_my <- function(data_seurat, expr_matrix, sample, spot, gene, neighborhood_size){
  ## Cycle through samples, spots, and genes to calculate mean expression in a neighborhood of fixed size
  
  # Find the coordinates of the center spot
  center_coordinates <- data_seurat@images$slice1@coordinates[spot, c("row", "col")]
  
  # Find all spot names in the neighborhood
  coordinate_list <- data_seurat@images$slice1@coordinates %>%
    dplyr::filter(row >= center_coordinates[[1]] - neighborhood_size &
                    row <= center_coordinates[[1]] + neighborhood_size &
                    col >= center_coordinates[[2]] - neighborhood_size &
                    col <= center_coordinates[[2]] + neighborhood_size) %>%
    rownames()
  
  return( c(sample, spot, gene, round(median(expr_matrix[gene, coordinate_list]), 2)) )
}

# Store mean expression here
expr_mean_list <- list()

# Cycle through sample list
for(sample in names(data_seurat_list)){
  filename <- paste0(output_path, "table_", sample, "_", assay, "_", slot, "_", neighborhood_size, ".txt")
  
  if(!file.exists(filename)){
    
    # Define sample of interest
    data_seurat <- data_seurat_list[[sample]]
    # Extract appropriate expression matrix
    expr_matrix <- as.matrix(Seurat::GetAssayData(data_seurat, assay=assay, slot=slot))
    
    # Parallelize the process for speed
    cl <- parallel::makeCluster(parallel::detectCores())
    doParallel::registerDoParallel(cl)
    
    # Cycle through samples, spots, and genes to calculate mean expression in a neighborhood of fixed size    
    expr_mean <- foreach(spot = colnames(data_seurat), .combine=c) %:%
      foreach (gene = gene_list, .combine=data.frame) %dopar% 
      {
        neighborhood_expr_mean_my(data_seurat, expr_matrix, sample, spot, gene, neighborhood_size)
      }
    
    parallel::stopCluster(cl)
    
    expr_mean_list[[sample]] <- tibble::as_tibble(t(data.frame(expr_mean)))
    colnames(expr_mean_list[[sample]]) <- c("Sample_Name", "Coordinate", "Symbol", "ExpressionMedian")
    
    # Write data to file
    readr::write_delim(expr_mean_list[[sample]], filename, delim="\t")
    
  } else{
    
    expr_mean_list[[sample]] <- readr::read_delim(filename, delim="\t")
    
  }
}


## MATCH COORDINATES ----


# Sample1 Seurat barcode names
sample1_coordinate_names <- rownames(data_seurat_list[[sample1_name]]@images$slice1@coordinates)
sample1_coordinate_temp <- c(sapply(sample1_coordinate_names, function(x) c(paste(x, "UL"), paste(x, "LR"))))

# Extract nodes for reference
#coregistration_xml_nodes <- XML::getNodeSet(XML::xmlParse(coregistration_xml), "//Annotation")

# Grab all annotation nodes
coregistration_xml_all <- xml2::xml_find_all(coregistration_xml, ".//Annotation")

# Grab the first sample annotation of the two co-registered images
coregistration_xml_first <- coregistration_xml_all[[which(sapply(coregistration_xml_all, function(x) grepl(paste0("\\b", gsub("_T$", "", sample1_name), "\\b"), x)))]]

# XML coordinates to tibble
coregistration_df <- xml2::xml_find_all(coregistration_xml_first, ".//Vertices/V") %>% 
  purrr::map_df(function(x) { list(X = xml2::xml_attr(x, "X"), Y = xml2::xml_attr(x, "Y")) }) %>%
  dplyr::mutate(Barcode_Sample1 = sample1_coordinate_temp) %>%
  dplyr::relocate(Barcode_Sample1) %>%
  df_wide2long_my(key="Corner", val="Value") %>%
  dplyr::mutate(Corner = paste0(gsub(".+ ", "", Barcode_Sample1), Corner),
                Barcode_Sample1 = gsub(" .+", "", Barcode_Sample1),
                Value = as.numeric(Value)) %>%
  df_long2wide_my(rows="Barcode_Sample1", cols="Corner", value="Value") %>%
  dplyr::group_by(Barcode_Sample1) %>%
  dplyr::mutate(X = round(min(LRX, ULX) + abs(ULX - LRX)/2),
                Y = round(min(ULY, LRY) + abs(ULY - LRY)/2))

# Find spot coordinates for the second sample
sample2_coordinates <- tibble::as_tibble(data_seurat_list[[sample2_name]]@images$slice1@coordinates[, c("imagerow", "imagecol")], 
                                         rownames="Barcode_Sample2") %>%
  dplyr::rename(X = imagecol, Y = imagerow) %>%
  dplyr::arrange(X, Y)

# Match HALO co-registration coordinate with the nearest 10X spot
sample1_to_sample2 <- nearest_neighbor_df_my(sample2_coordinates, coregistration_df, value="Barcode_Sample1") %>%
  dplyr::select(Barcode_Sample1, Barcode_Sample2, Distance)

# Enumerate the barcodes from both samples in the same way
sample1_to_sample2[["Order"]] <- 1:nrow(sample1_to_sample2)


## VISUALIZE COREGISTRATION BY SPOT ----


# Add the barcode order to both samples
data_seurat_list[[sample1_name]]@meta.data["Order"] <- sample1_to_sample2[match(rownames(data_seurat_list[[sample1_name]]@meta.data), sample1_to_sample2$Barcode_Sample1), "Order"]
data_seurat_list[[sample1_name]]@meta.data[which(is.na(data_seurat_list[[sample1_name]]@meta.data["Order"])), "Order"] <- -2000
data_seurat_list[[sample2_name]]@meta.data["Order"] <- sample1_to_sample2[match(rownames(data_seurat_list[[sample2_name]]@meta.data), sample1_to_sample2$Barcode_Sample2), "Order"]

# Filter badly matched barcodes
sample1_to_sample2 <- sample1_to_sample2 %>%
  dplyr::filter(Distance <= data_seurat_list[[sample2_name]]@images$slice1@scale.factors$spot_diameter_fullres)

cat("Number of barcodes in the first sample =", length(sample1_coordinate_names), "\n",
    "Number of barcodes in the co-registration with first sample as reference =", nrow(coregistration_df), "\n",
    "Number of barcodes in the second sample =", nrow(sample2_coordinates), "\n",
    "Number of successfully co-registered items =", nrow(sample1_to_sample2), "\n",
    "Number of barcodes in the first sample preserved =", length(intersect(sample1_to_sample2$Barcode_Sample1, expr_mean_list[[1]]$Coordinate)), "\n",
    "Number of barcodes in the second sample preserved =", length(intersect(sample1_to_sample2$Barcode_Sample2, expr_mean_list[[2]]$Coordinate)), "\n")

# Visualize co-registered barcodes
p <- spatial_feature_plot_my(data_seurat_list[[sample1_name]], "Order", min.cutoff="q0", max.cutoff="q100", name=sample1_name)
filename <- paste0(output_path, "coreg_", sample1_name)
write_plot2file_my(p, filename)

p <- spatial_feature_plot_my(data_seurat_list[[sample2_name]], "Order", min.cutoff="q0", max.cutoff="q100", name=sample1_name)
filename <- paste0(output_path, "coreg_", sample2_name)
write_plot2file_my(p, filename)


## CALCULATE CORRELATION FOR EACH GENE ----


# Join sample1 and sample2 expression values using the co-registration map
expr_mean_df <- expr_mean_list[[1]] %>%
  dplyr::select(Sample_Name, Coordinate, Symbol, ExpressionMedian) %>%
  dplyr::rename(Barcode_Sample1 = Coordinate) %>%
  dplyr::inner_join(sample1_to_sample2, by="Barcode_Sample1") %>%
  dplyr::inner_join(expr_mean_list[[2]] %>%
                      dplyr::select(Sample_Name, Coordinate, Symbol, ExpressionMedian) %>%
                      dplyr::rename(Barcode_Sample2 = Coordinate), by=c("Barcode_Sample2", "Symbol")) %>%
  dplyr::mutate(ExpressionMedian.x = as.numeric(ExpressionMedian.x),
                ExpressionMedian.y = as.numeric(ExpressionMedian.y))

# Calculate correlation between two samples for every gene
cor_list <- list()
for(gene in gene_list){
  expr_mean_loc_df <- expr_mean_df %>%
    dplyr::filter(Symbol == gene)
  
  if(max(expr_mean_loc_df$ExpressionMedian.x) > 0 && max(expr_mean_loc_df$ExpressionMedian.y) > 0){
    cor_list[[gene]] <- round(cor(expr_mean_loc_df$ExpressionMedian.x, expr_mean_loc_df$ExpressionMedian.y), 3)
  }
}

# Create correlation tibble
cor_df <- tibble::tibble(Symbol = names(cor_list), R = unlist(unname(cor_list))) %>%
  dplyr::mutate(Correlated = R >= 1)

# Annotate gene information
expr_mean_df <- expr_mean_df %>%
  dplyr::inner_join(gene_sig_df, by="Symbol")


## VISUALIZE CORRELATIONS IN HISTOGRAMS ----


filename <- paste0(output_path, "cor_hist_", cohort_name, "_", assay, "_", slot, "_", neighborhood_size)
create_hist_plot_my(cor_df, x_label="R", fill_label="Correlated", intercept=c(-1, 0, 1), binwidth=0.1, 
                    filename=filename, labels=NULL)

# Select genes expressed well in both samples
genes_select <- expr_mean_df %>% 
  dplyr::group_by(Symbol) %>% 
  dplyr::summarise(Sample1_Count = sum(ExpressionMedian.x >= 1),
                   Sample2_Count = sum(ExpressionMedian.y >= 1)) %>%
  dplyr::filter(Sample1_Count >= 10 & Sample2_Count >= 10) %>%
  dplyr::pull(Symbol)

# Filter down correlations to only genes that are well expressed in both samples
cor_df <- cor_df %>%
  dplyr::filter(Symbol %in% genes_select)

filename <- paste0(output_path, "cor_hist_select_", cohort_name, "_", assay, "_", slot, "_", neighborhood_size)
create_hist_plot_my(cor_df, x_label="R", fill_label="Correlated", intercept=c(-1, 0, 1), binwidth=0.1, 
                    filename=filename, labels=NULL)


## VISUALIZE CORRELATIONS IN SCATTER PLOTS ----


# 10 most positively correlated genes
gene_pos_list <- cor_df %>%
  dplyr::arrange(desc(R)) %>%
  dplyr::slice_head(n=20)

# 10 genes with correlation closest to zero
gene_neu_list <- cor_df %>%
  dplyr::filter(!Symbol %in% gene_pos_list$Symbol) %>%
  dplyr::mutate(R = abs(R)) %>%
  dplyr::arrange(R) %>%
  dplyr::slice_head(n=20)

# 10 most negatively correlated genes
gene_neg_list <- cor_df %>%
  dplyr::filter(!Symbol %in% c(gene_pos_list$Symbol, gene_neu_list$Symbol)) %>%
  dplyr::arrange(R) %>%
  dplyr::slice_head(n=20)

# Create scatter plot for 10 most positively correlated genes
filename <- paste0(output_path, "scatter_pos_", cohort_name, "_", assay, "_", slot, "_", neighborhood_size)
create_scatter_plot_my(expr_mean_df %>%
                         dplyr::filter(Symbol %in% gene_pos_list$Symbol), x_label="ExpressionMedian.x", y_label="ExpressionMedian.y", 
                       fill_label="Sig_Name", facet_var=c("Symbol", "fixed"), filename=filename, 
                       labels=c(sample1_name, sample2_name, paste(cohort_name, assay, slot)), do_fit="log")

# Create scatter plot for 10 most negatively correlated genes
filename <- paste0(output_path, "scatter_neg_", cohort_name, "_", assay, "_", slot, "_", neighborhood_size)
create_scatter_plot_my(expr_mean_df %>%
                         dplyr::filter(Symbol %in% gene_neg_list$Symbol), x_label="ExpressionMedian.x", y_label="ExpressionMedian.y", 
                       fill_label="Sig_Name", facet_var=c("Symbol", "fixed"), filename=filename, 
                       labels=c(sample1_name, sample2_name, paste(cohort_name, assay, slot)), do_fit="log")

# Create scatter plot for 10 genes with correlation closest to zero
filename <- paste0(output_path, "scatter_neu_", cohort_name, "_", assay, "_", slot, "_", neighborhood_size)
create_scatter_plot_my(expr_mean_df %>%
                         dplyr::filter(Symbol %in% gene_neu_list$Symbol), x_label="ExpressionMedian.x", y_label="ExpressionMedian.y", 
                       fill_label="Sig_Name", facet_var=c("Symbol", "fixed"), filename=filename, 
                       labels=c(sample1_name, sample2_name, paste(cohort_name, assay, slot)), do_fit="log")


## VISUALIZE CORRELATIONS IN SPATIAL PLOTS ----


# Create spatial plots for genes above
for(gene in gene_pos_list$Symbol){
  
  filename <- paste0(output_path, "spatial_", cohort_name, "_", assay, "_", slot, "_", neighborhood_size, "_pos_", gene)
  batch_spatial_feature_plot_my(data_seurat_list, gene, output_file=filename)
}

# Create spatial plots for genes above
for(gene in gene_neg_list$Symbol){
  
  filename <- paste0(output_path, "spatial_", cohort_name, "_", assay, "_", slot, "_", neighborhood_size, "_neg_", gene)
  batch_spatial_feature_plot_my(data_seurat_list, gene, output_file=filename)
}

# Create spatial plots for genes above
for(gene in gene_neu_list$Symbol){
  
  filename <- paste0(output_path, "spatial_", cohort_name, "_", assay, "_", slot, "_", neighborhood_size, "_neu_", gene)
  batch_spatial_feature_plot_my(data_seurat_list, gene, output_file=filename)
}
