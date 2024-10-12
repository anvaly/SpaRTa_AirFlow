# Author: Anna Lyubetskaya. Date: 21-02-20

# Evaluate data trends within ST datasets:
# - Create a UMI and number of gene spatial plots including the spots outside of the tissue segment
# - Number of UMIs and genes detected in each spot under, adjacent, and outside of the tissue segment
# - Mean gene expression vs number of spots in which the gene is detected under, adjacent, and outside the tissue segment
# - Define a list of promiscuious genes that are expressed in >= X spots outside of the tissue segment
# This script runs from the original 10X output and is an attempt at evaluating ambient signal in an intuitive way


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
source("code/R/Utils/utils_colormap.R")

source("code/R/b_Cluster_Biomarkers/Seurat_cluster_utils.R")


## PARAMETERS ----


# Cohort name
cohort_name <- NULL


## PATHS ----


# Input file with sample meta data
meta_file <- "XXXX"

# Path to gene signatures
sig_path <- "data/import/Signatures/signatures_pdac_t1000_May21.txt"

# Output folder
output_path_init <- "XXXX"
# output_path_init <- "XXXX"

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)

# If probes were removed in the 10x pipeline, and you don't want to go back and re-run you can restore those probes from the raw data with this flag
return_probes <- FALSE
# File with probe offtarget information
probe_file <- "XXXX"


## INGEST DATA ----


# Read sample file
cohort_meta_df <- readr::read_delim(file=meta_file, delim="\t") # %>% 
  # dplyr::filter(grepl(cohort_name, Sample_Name) & !grepl("FFPE-probes", Protocol))


if(!is.null(sig_path)){
  # Load signatures and filter them down to only well represented genes
  signature_list <- read_filter_signatures_my(sig_path, NULL, sig_length_min=1, sig_length_max=1000, ratio_threshold=0)
  gene_sig_df <- invert_list_my(signature_list)
} else{
  gene_sig_df <- tibble::tibble(Symbol = NA)
}


## RAW DATA CHECK ----


for(i in c(33:nrow(cohort_meta_df))){
  
  ## INGEST 10X DATA ----
  
  
  # Input data path
  filename_in <- cohort_meta_df[[i, "FullPath"]]
  # Sample name
  sample_name <- cohort_meta_df[[i, "Sample_Name"]]
  
  # Create a sample-specific output folder
  output_path <- paste0(output_path_init, sample_name, "/")
  dir.create(output_path, showWarnings = FALSE)
  
  # Ingest raw (unfiltered) data
  data_seurat_raw <- read_10X_spatial_folder_my(filename_in, filename="raw_feature_bc_matrix.h5")
  
  # Ingest filtered data
  data_seurat <- read_10X_spatial_folder_my(filename_in)
  
  # Update filtered data if probes were filtered out in the pipeline
  if(return_probes == TRUE){
    data_seurat <- restore_filtered_probes(data_seurat_raw, data_seurat, probe_file)
  }

  
  ## REMOVE UNEXPRESSED GENES ----
  
  
  # Non-zero genes in the unfiltered data
  gene_list <- names(which(Matrix::rowSums(data_seurat_raw@assays$Spatial@counts) > 0))
  
  # Subset raw and filtered data
  data_seurat_raw <- subset(data_seurat_raw, features=gene_list)
  data_seurat <- subset(data_seurat, features=gene_list)
  

  ## EXCLUDE TISSUE ADJACENT SPOTS ----


  # 10X filtered spot list
  spot_under_tissue <- colnames(data_seurat)

  # Extract all Seurat coordinates
  xy <- tissue_position_list[, c("row", "col")]
  barcode_list <- tissue_position_list$tissue
  xy <- mapply(xy, FUN=as.numeric)

  # Adjacency distance
  distance <- 2

  # Find the matrix with the indices of points belonging to the set of the k nearest neighbours of each other
  spot_nb <- spdep::dnearneigh(xy, 0, distance, row.names=barcode_list)
  names(spot_nb) <- barcode_list

  # Find all spots within or adjacent to tissue
  spot_under_adjacent_tissue <- unique(unlist(unname(sapply(spot_nb[spot_under_tissue], function(x) barcode_list[x]))))
  # Check that all spots under tissue got included in the neighborhood analysis
  setdiff(spot_under_tissue, spot_under_adjacent_tissue)

  # Spots outside of tissue
  spot_outside_tissue <- setdiff(colnames(data_seurat_raw), spot_under_adjacent_tissue)
  # Spots adjacent to tissue
  spot_adjacent_tissue <- setdiff(spot_under_adjacent_tissue, spot_under_tissue)


  ## ANALYZE COORDINATE PROPERTIES ----


  # Expression long tibble
  data_meta_df <- tibble::as_tibble(data_seurat_raw@meta.data, rownames = "Coordinate")

  # UMI sum per spot, including filtered spots
  coordinate_df <- data_meta_df %>%
    dplyr::mutate(UnderTissue = Coordinate,
                  nCount_Spatial_log10 = log10(nCount_Spatial))

  # Assign spots into 3 buckets - under / adjacent / outside of tissue
  coordinate_df[which(coordinate_df$Coordinate %in% spot_under_tissue), "UnderTissue"] <- "UnderTissue"
  coordinate_df[which(coordinate_df$Coordinate %in% spot_adjacent_tissue), "UnderTissue"] <- "AdjacentTissue"
  coordinate_df[which(coordinate_df$Coordinate %in% spot_outside_tissue), "UnderTissue"] <- "OutsideTissue"

  # Factorize the tissue localization label
  coordinate_df["UnderTissue"] <- factor(coordinate_df[["UnderTissue"]], levels=c("UnderTissue", "AdjacentTissue", "OutsideTissue"))

  # Table of spot frequences by type
  spot_type_df <- tibble::as_tibble(data.frame(table(coordinate_df$UnderTissue))) %>%
    dplyr::rename(UnderTissue = Var1)


  # Violin plot of library sizes by spot type
  p1 <- create_violin_plot_my(coordinate_df, x_label="UnderTissue", y_label="nCount_Spatial_log10",
                              fill_label="orig.ident", filename=NULL,
                              labels=c("Spot", "Library Size (Sum of UMIs, log10)", "Library size of tissue-occupied and tissue-free spots"))

  # Violin plot of feature number by spot type
  p2 <- create_violin_plot_my(coordinate_df, x_label="UnderTissue", y_label="nFeature_Spatial",
                              fill_label="orig.ident", filename=NULL,
                              labels=c("Spot", "Number of features", "Library size of tissue-occupied and tissue-free spots"))

  # Write figure to file
  filename <- paste0(output_path, "violin_by_spot_type_", sample_name)
  write_plot2file_my(patchwork::wrap_plots(list(p1, p2), ncol=1), filename, num_col=1, num_row=2)


  ## ANALYZE EXPRESSION MEAN V SPOT-OCCUPANCY ----


  # Normalize all data - use LogNorm for simplicity of interpretation
  data_seurat_raw <- Seurat::NormalizeData(data_seurat_raw, normalization.method = "LogNormalize")

  # Extract raw UMIs for every spot in and out of tissue
  expr_all_umi_df <- seurat_expression_to_long_tibble_my(data_seurat_raw, assay="Spatial", slot="counts")
  table(expr_all_umi_df$Spatial_counts > 0)
  
  expr_all_umi_df <-  expr_all_umi_df %>%
    dplyr::filter(Spatial_counts > 0) %>%
    dplyr::rename(UMI = Spatial_counts)

  # Extract normalized UMIs for every spot in and out of tissue
  expr_all_norm_df <- seurat_expression_to_long_tibble_my(data_seurat_raw, assay="Spatial", slot="data")
  table(expr_all_norm_df$Spatial_data > 0)
  
  expr_all_norm_df <-  expr_all_norm_df %>%
    dplyr::filter(Spatial_data > 0) %>%
    dplyr::rename(LogNormUMI = Spatial_data)

  # Mean and spot-occupancy
  gene_expr_df <- expr_all_umi_df %>%
    dplyr::inner_join(expr_all_norm_df, by=c("Symbol", "Coordinate")) %>%
    dplyr::inner_join(coordinate_df, by=c("Coordinate")) %>%
    dplyr::group_by(Symbol, UnderTissue) %>%
    dplyr::summarise(CountPositive = dplyr::n_distinct(Coordinate),
                     ExpressionMeanLogUMI = round(mean(log2(UMI)), 3),
                     ExpressionMeanLogNormUMI = round(mean(LogNormUMI), 3)) %>%
    dplyr::inner_join(spot_type_df, by="UnderTissue") %>%
    dplyr::mutate(InSignature = Symbol %in% gene_sig_df$Symbol,
                  PercentPositive = round(CountPositive / Freq * 100, 2),
                  IsMitoRibo = grepl("^MT-|^RP[SL]", Symbol)) %>%
    dplyr::arrange(IsMitoRibo, UnderTissue, desc(PercentPositive), desc(ExpressionMeanLogUMI))

  # Write data to file
  filename <- paste0(output_path, "table_mean_v_occupancy_", sample_name, ".txt")
  readr::write_delim(gene_expr_df, filename, delim="\t")


  # Mean expression (in UMIs) vs spot occupancy
  p1 <- create_scatter_plot_my(gene_expr_df %>%
                                 dplyr::filter(IsMitoRibo == FALSE),
                               x_label="PercentPositive", y_label="ExpressionMeanLogUMI",
                               fill_label="InSignature", shape=19, size=0.5, facet_var=c("UnderTissue", "fixed"),
                               labels=c("Percent Spot with Detectable Gene Levels", "Expression, mean log2 UMI", sample_name),
                               filename=NULL, do_fit=FALSE)

  # Mean expression (in log norm UMIs) vs spot occupancy
  p2 <- create_scatter_plot_my(gene_expr_df %>%
                                 dplyr::filter(IsMitoRibo == FALSE),
                               x_label="PercentPositive", y_label="ExpressionMeanLogNormUMI",
                               fill_label="InSignature", shape=19, size=0.5, facet_var=c("UnderTissue", "fixed"),
                               labels=c("Percent Spot with Detectable Gene Levels", "Expression, mean log2 library size normalized UMI", sample_name),
                               filename=NULL, do_fit=FALSE)

  # Write figure to file
  filename <- paste0(output_path, "/scatter_mean_v_occupancy_", sample_name)
  write_plot2file_my(patchwork::wrap_plots(list(p1, p2), ncol=2), filename, num_col=2, num_row=3)


  ## COMPARE EXPRESSION UNDER AND OUTSIDE TISSUE ----


  # Create a tibble comparing various parameters for each gene under and outside the tissue
  gene_expr_wide_list <- list()
  for(m in c("ExpressionMeanLogNormUMI", "ExpressionMeanLogUMI", "PercentPositive")){

    gene_expr_wide_list[[m]] <- gene_expr_df %>%
      df_long2wide_my(rows="Symbol", cols="UnderTissue", value=m) %>%
      tidyr::replace_na(list(OutsideTissue = 0, UnderTissue = 0, AdjacentTissue = 0)) %>%
      dplyr::inner_join(gene_expr_df %>%
                          dplyr::select(Symbol, InSignature) %>%
                          unique, by="Symbol") %>%
      dplyr::mutate(Measurement = m) %>%
      dplyr::arrange(desc(OutsideTissue), desc(UnderTissue))

    # Write data to file
    filename <- paste0(output_path, "table_measure_", m, "_", sample_name, ".txt")
    readr::write_delim(gene_expr_wide_list[[m]], filename, delim="\t")
  }

  # Mean expression (in log norm UMIs) vs spot occupancy
  filename <- paste0(output_path, "scatter_under_outside_tissue_", sample_name)
  create_scatter_plot_my(dplyr::bind_rows(gene_expr_wide_list),
                         x_label="UnderTissue", y_label="OutsideTissue",
                         fill_label="InSignature", shape=19, size=0.5, facet_var=c("Measurement", "free"),
                         labels=c("Under Tissue", "Outside Tissue", sample_name),
                         filename=filename, do_fit=FALSE)
}
