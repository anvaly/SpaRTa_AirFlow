# Author: Anna Lyubetskaya. Date: 21-01-03
# Annotate biomarkers derived from clustering by signatures and gene over-representation analysis


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_heatmap.R")
source("code/utils/utils_signatures.R")
source("code/R/Utils/utils_10X_gsea.R")
source("code/R/Utils/utils_10X_matrix.R")


## PARAMETERS ----


cohort_expr <- "PDAC108_path14_harmonyepi_rpca_sct"
file_suffix <- "markers_significant"

# Parameters
reference_db <- "org.Hs.eg.db"  # org.Rn.eg.db, org.Hs.eg.db, org.Mm.eg.db
ref_short <- "human"   # rat, human, mouse
ref_twoletter <- "Hs"  # Hs, Mm, Rn
gsea_qvalue_threshold <- 0.01
gsea_refs <- c("MSigDB")  # "KEGG"

# Threshold for ambiguous genes
ambiguous_threshold <- 5
# Number of top genes to gather
topX <- 30


# Enforce specific constraints on key cohorts
resolution <- NULL
if(cohort_expr == "PDAC108_path14_harmonyepi_rpca_sct"){
  
  # User defined cluster order
  cluster_order <- c(3, 10, 0, 4, 6, 9)
  resolution <- "integrated_snn_res.0.4"
}


## PATHS ----


# Path to a signature file
sig_path <- "data/import/Signatures/signatures_pdac_20230720.txt"

# Input folder
input_path <- "XXXX"

# Input folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_expr, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, NULL, sig_length_min=1, sig_length_max=1000, ratio_threshold=0)


## WRANGLE DATA ----


# Invert signatures to get an annotated gene list
sig_gene_df <- invert_list_my(signature_list)

# Find the file with the DE gene corresponding to the resolution of interest
file_list <- dir(input_path, pattern=paste0(file_suffix, ".*", cohort_expr, ".*", ".txt"), full.names=TRUE, recursive=TRUE)


## DEFINE BACKGROUND GENE LIST ----


# Reference file - assuming its name can be derived from the first file in the list
ref_name <- dir(input_path, full.names = TRUE)[grep("ref_gene_list", dir(input_path, full.names = TRUE))]

# Read the reference file
if(length(ref_name) > 0){
  ref_df <- readr::read_delim(ref_name, delim="\t")
} else{
  ref_df <- readr::read_delim(paste0("data/import/References/", ref_twoletter, ".txt"), delim="\t")
}

# MT/ribo gene list
exclude_list <- unique(c(ref_df$Symbol[identify_mt_ribo_genes_my(ref_df$Symbol)],
                         identify_unannotated_genes_my(ref_df$Symbol)))


## ANNOTATE GENE LISTS ----


# Join Seurat markers with signature annotation
# Perform gene over-representation analysis on the significant markers
for(f in file_list){
  
  # Left join signature annotations from references to biomarker information from clustering
  markers_df <- readr::read_delim(f, delim="\t") %>%
    dplyr::left_join(sig_gene_df, by="Symbol") %>%
    tidyr::replace_na(list(InSignature = FALSE, Num_Sigs = 0))
  
  print("These biomarkers are not found in the reference:")
  warning(setdiff(markers_df$Symbol, ref_df$Symbol))
  
  
  # Redirect outputs to a new folder
  f <- gsub(".+/", output_path, f)
  
  
  # Make sure all marker genes are in the reference
  markers_df <- markers_df %>%
    dplyr::filter(Symbol %in% ref_df$Symbol)
  
  if(!is.null(gsea_refs)){
    # Add gene over-representation data
    output_file <- gsub(".txt", "_gsea", f)
    
    enrichment_df <- markers_cluster_profile_my(markers_df, ref_df$Symbol, reference_db, ref_short, 
                                                gsea_qvalue_threshold, gsea_refs, output_path=output_file) %>%
      dplyr::mutate(ID = gsub("HALLMARK_", "", ID)) %>%
      dplyr::group_by(ID) %>%
      dplyr::mutate(IDMaxScore = max(nlog10qvalue)) %>%
      dplyr::ungroup()
    
    # Select top enriched patywas for vis
    enrichment_filt_df <- enrichment_df %>%
      dplyr::group_by(cluster) %>%
      dplyr::top_n(5, IDMaxScore) %>%
      dplyr::ungroup()
    
    # Setup heatmap parameters
    params <- list(
      cell_value = "nlog10qvalue",
      row_label = "ID", 
      col_label = "cluster", 
      distance = "pearson",
      row_annotation = NULL,
      col_annotation = NULL,
      range = c(0, round(max(enrichment_filt_df$nlog10qvalue), -1) / 4, round(max(enrichment_filt_df$nlog10qvalue), -1)),
      colors = c("white", "blue3", "navy")
    )
    
    
    # Create a heatmap of top enriched pathways
    filename <- paste0(output_file, "_hm.png")
    hm <- create_heatmap_my(enrichment_filt_df, params, filename=filename, height=round(length(unique(enrichment_filt_df$ID))/2), 
                            width=length(unique(enrichment_filt_df$cluster))+1)
    
    
    # Set up bespoke vis parameters
    if(!is.null(resolution) && grepl(gsub("integrated_snn_res.", "", resolution), f)){
      enrichment_filt_df <- enrichment_filt_df %>%
        dplyr::filter(cluster %in% cluster_order) %>%
        dplyr::mutate(ID = ifelse(nchar(ID) > 45, paste0(substr(ID, 1, 45), "..."), ID))
      
      params$column_order <- as.character(cluster_order)
      
      filename <- paste0(output_file, "_hm_bespoke.png")
      hm <- create_heatmap_my(enrichment_filt_df, params, filename=filename, height=round(length(unique(enrichment_filt_df$ID))/2), 
                              width=length(unique(enrichment_filt_df$cluster))+1)
    }
    
  }
  
  # Sort annotated markers before writing them to file
  markers_df <- markers_df %>%
    dplyr::filter(!Symbol %in% exclude_list) %>%
    dplyr::arrange(cluster, desc(InSignature), desc(direction), desc(p_val_adj_neg_log10))
  
  # Find markers that appear in multiple clusters to remove them from top10
  markers_stat_df <- markers_df %>%
    dplyr::group_by(Symbol, direction) %>%
    dplyr::summarise(ClusterCount = dplyr::n_distinct(cluster)) %>%
    dplyr::ungroup()
  
  # Add number of clusters to which every gene belongs as an additional field
  markers_df <- markers_df %>%
    dplyr::inner_join(markers_stat_df, by=c("Symbol", "direction"))
  
  # Output file name
  output_file <- gsub(file_suffix, paste0(file_suffix, "_annotated"), f)
  readr::write_delim(markers_df, output_file, delim="\t")
  
  # Find only unique markers
  markers_unique_df <- markers_df %>%
    dplyr::filter(ClusterCount == 1)
  
  # Output file name
  output_file <- gsub(file_suffix, paste0(file_suffix, "_annotated_unique"), f)
  readr::write_delim(markers_unique_df, output_file, delim="\t")
  
  
  # Remove ambiguous genes
  # Select only upregulated genes for the top10 lists
  # Select top N genes for each cluster for easier manual verification
  markers_top10_df <- markers_df %>%
    dplyr::filter(direction == "UP") %>%
    dplyr::filter(ClusterCount <= ambiguous_threshold) %>%
    dplyr::mutate(pct_approximate = round(pct_1/3, -1)*3) %>%
    dplyr::group_by(cluster, InSignature, direction) %>%
    dplyr::arrange(desc(pct_approximate), desc(p_val_adj_neg_log10)) %>%
    dplyr::slice_head(n=topX) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(cluster, desc(InSignature), desc(direction), desc(p_val_adj_neg_log10))
  
  # Output file name
  output_file <- gsub(file_suffix, paste0(file_suffix, "_annotated_top10"), f)
  readr::write_delim(markers_top10_df, output_file, delim="\t")
  
  # Remove ambigious genes
  # Select only upregulated genes for the top10 lists
  # Select top N genes for each cluster for easier manual verification
  markers_top10_df <- markers_unique_df %>%
    dplyr::filter(direction == "UP") %>%
    dplyr::mutate(pct_approximate = round(pct_1/3, -1)*3) %>%
    dplyr::group_by(cluster, InSignature, direction) %>%
    dplyr::arrange(desc(pct_approximate), desc(p_val_adj_neg_log10)) %>%
    dplyr::slice_head(n=topX) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(cluster, desc(InSignature), desc(direction), desc(p_val_adj_neg_log10))
  
  # Output file name
  output_file <- gsub(file_suffix, paste0(file_suffix, "_annotated_unique_top10"), f)
  readr::write_delim(markers_top10_df, output_file, delim="\t")
}
