# Author: Anna Lyubetskaya. Date: 23-07-27
# Parse WGCNA results


## ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`


## PARAMETERS ----


run_name <- "PDAC84_path12"


## PATHS ----


# Path to a signature file
input_path <- "XXXX"

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, run_name, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Seurat data
wgcna_data <- readRDS(input_path) %>%
  dplyr::mutate(OS_Hazard_Ratios = round(OS_Hazard_Ratios, 2),
                PFI_Hazard_Ratios = round(PFI_Hazard_Ratios, 2),
                OS_Sig = OS_Padj <= 0.05,
                PFI_Sig = PFI_Padj <= 0.05) 
# %>%
#   dplyr::select(OS_Hazard_Ratios, OS_Sig, PFI_Hazard_Ratios, PFI_Sig,
#                 best_H, best_C3, best_C6, best_C7, best_C8,
#                 cell_type, Enrichment_Signature, Overlap_with_targets)

# Write the base table to a file
filename <- paste0(output_path, "wgcna_result.txt")
readr::write_delim(wgcna_data, filename, delim="\t")
