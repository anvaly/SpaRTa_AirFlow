# Author: Anna Lyubetskaya. Date: 22-02-10

# Calculate signature score and the corresponding empirical p-value
# This takes a while


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/R/Utils/utils_10X_signatures.R")
source("code/R/Utils/utils_10X_matrix.R")

library(Seurat)


## PARAMETERS ----


# The cohort of interest: regex or name
cohort_regex <- "PDAC"

# Specific signatures to use
sig_select <- c("PDAC.collisson.classical", "PDAC.collisson.quasimesenchymal",
                "PDAC.P19.Fibroblast", "PDAC.moffitt.activatedstroma", "PDAC.moffitt.normalstroma",
                "PDAC.P19.Macrophage", "BMS.MHCI_AP", "PDAC.U.Pathway.MHCII",
                "BMS.WZ.CD8.effector.CuiCell2021", "PDAC.U.Immune.Tcell.CD8",
                "BMS.Pathway.IFNa", "BMS.Pathway.IFNg", "BMS.Pathway.TGFB", "BMS.NFKB2_NIK_pathway1")


## PATHS ----


# Input folder
#input_path <- "XXXX"
input_path <- "XXXX"

# Output path
output_path <- "XXXX"

# Create output folder
dir.create(output_path, showWarnings = FALSE)

# Path to a signature file
sig_path <- "data/import/Signatures/signatures_220210.txt"


## INGEST DATA ----


# Find all files following the regex
file_list <- dir(input_path, pattern=paste0(cohort_regex, ".*.rds"), full.names=TRUE)


## CALCULATE SIGNATURE SCORE ----


# Loop through Seurat objects
for(f in file_list){
  
  # Ingest the Seurat object
  data_seurat <- readRDS(f)
  
  # Find genes abundant in this sample
  gene_list <- seurat_select_abundant_genes_my(data_seurat)
  
  # Sample name
  sample_name <- data_seurat@misc$user.Sample_Name
  
  # Load signatures and filter them down to only well represented genes
  signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select)
  
  # Add signature scores to a seurat object
  data_seurat <- add_signature_scores_my(data_seurat, signature_list)

  # Find column names in the Seurat object - Seurat can rename original signatures
  sig_names <- colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))])

  # Loop through signatures
  for(sig in sig_names){
    # Write/read random gene expression scores to file
    filename_random <- paste0(output_path, "/sig_random_", sample_name, "_", sig, ".txt")
    
    # Create a random signature background distribution and score the target signature against it
    data_seurat <- signature_empirical_pvalue_my(data_seurat, signature_list[[gsub("sig.", "", sig)]], gene_list, output_path, sample_name, sig, 
                                                 sig_random_filename=filename_random, n_simulations=1000, assay="SCT", output_col_name=paste0(sig, "_pv"))
  }

  # Find column names in the Seurat object - Seurat can rename original signatures
  sig_names <- colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))])
  
  # Extract coordinate meta data
  meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate") %>%
    dplyr::mutate(user.Sample_Name = sample_name) %>%
    dplyr::select(dplyr::all_of(c("Coordinate", "user.Sample_Name", sig_names))) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), round, 3))
  
  # Write stats to file
  filename <- paste0(output_path, "/sig_", sample_name, ".txt")
  readr::write_delim(meta_df, filename, delim="\t")
  
}
