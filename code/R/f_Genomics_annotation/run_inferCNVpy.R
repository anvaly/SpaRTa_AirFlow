# Author: Anna Lyubetskaya. Date: 22-01-10
# Run inferCNV

# Python version of inferCNV: https://icbi-lab.github.io/infercnvpy/index.html
# Manual: https://github.com/broadinstitute/inferCNV/wiki
# Vignette: https://bioconductor.org/packages/devel/bioc/vignettes/infercnv/inst/doc/inferCNV.html

# Some standard input files are here: data/import/inferCNV/hg38_gencode_v27.txt


## SETUP PYTHON ----


## cd /XXXX
## Download miniforge and put it in the folder: https://pypi.org/project/infercnvpy/
## bash Mambaforge-Linux-x86_64.sh
## miniforge/bin/pip install infercnv


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/R/Utils/utils_10X_in_out.R")

# Configure reticulate with my own Python on Stash
library(reticulate)

Sys.setenv(RETICULATE_PYTHON="XXXX")
use_condaenv("XXXX",
             required = TRUE)

py_sp <- import("scanpy")
py_icp <- import("infercnvpy")
py_ad <- import("anndata")
py_p <- import("pickle")


## PARAMETERS ----


# Cohort name for outputs
cohort_name <- "PDAC84_path12_merge"

# Establish the control pathology group
groups_ref <- c("Lymph_Node", "TLS_Mature", "TLS_Aggregate", "TLS_Immature", "Adipose",
                "Blood_Vessels", "Muscle", "Nerve", "Blood", "Exocrine_Endocrine")


## PATHS ----


# Input folder
input_path <- "XXXX"

# inferCNV data
inferCNV_file <- "data/import/inferCNV/hg38_gencode_v27.txt"

# The inferCNV reference
infercnv_ref_path <- "XXXX"

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, "/", cohort_name, "/")

dir.create(output_path_init, showWarnings=FALSE)
dir.create(output_path, showWarnings=FALSE)


## INGEST DATA ----


# Ingest a set of RDS Seurat objects
data_seurat <- readRDS(paste0(input_path, cohort_name, ".rds"))


## WRANGLE INPUTS ----


## Reference


# Ingest reference for inferCNV
ref_df <- readr::read_delim(infercnv_ref_path, delim="\t", col_names=FALSE)

colnames(ref_df) <- c("Symbol", "chromosome", "start", "end")

ref_df <- ref_df %>%
  dplyr::mutate(Symbol = toupper(Symbol)) %>%
  tibble::column_to_rownames("Symbol")

  
## Meta data


# Extract meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

# Find the max value across reference pathology classes
meta_df["inferCNVcontrol"] <- matrixStats::rowMaxs(as.matrix(meta_df[paste0("Pathology.", groups_ref, ".percent")]))

# Establishing the control and target groups
meta_df <- meta_df %>%
  dplyr::mutate(Group = ifelse(Pathology.Group %in% groups_ref & inferCNVcontrol >= 80, 
                               paste0(Pathology.Group, "_", user.Sample_Name), 
                               paste0("malignant_", Pathology.Group, "_", user.Sample_Name))) %>%
  dplyr::select(user.Block_ID, Pathology.Group, Coordinate, Group)

# Control groups
control_list <- unique(meta_df$Group[which(!grepl("malignant", meta_df$Group))])

table(meta_df$Group)


## INFERCNVPY CALL FUNCTION ----


run_infercnvpy_my <- function(data_seurat, meta_df, ref_df, control_list, block="data"){
  ## Wrangle data from R to Python and run inferncvpy through reticulate
  
  gc()
  
  table(meta_df$Coordinate == colnames(Seurat::GetAssayData(data_seurat, assay="SCT", slot="data")))
  
  # Create a Python AnnData object: X, obs, var
  py_adata <- py_ad$AnnData(Seurat::GetAssayData(data_seurat, assay="SCT", slot="data"), # X
                            ref_df[rownames(Seurat::GetAssayData(data_seurat, assay="SCT", slot="data")),],
                            meta_df %>% tibble::column_to_rownames("Coordinate")
  )
  
  # Transpose the object
  py_adata <- py_ad$AnnData$transpose(py_adata)
  
  # Run inferCNVpy - it augments the input object
  py_icp$tl$infercnv(py_adata, 
                     reference_key="Group", 
                     reference_cat=control_list, 
                     window_size=250)
  
  
  # Save infercnvpy result to an object
  filename <- paste0(output_path, block, "_res.RDS")
  py_adata$write(filename)

  # Write infercnvpy heatmap to file
  filename <- paste0(output_path, block, "_hm.png")
  png(file=filename, width=10, height=10, units="in", res=200)
  p <- py_icp$pl$chromosome_heatmap(py_adata, groupby="Pathology.Group")
  dev.off()
}


## SAMPLE ITERATOR ----


meta_df$user.Block_ID <- "all"

for(block in unique(meta_df$user.Block_ID)){

  print(block)
  
  # Subset meta data down to a single block
  meta_loc_df <- meta_df %>%
    dplyr::filter(user.Block_ID == block)  #  | Group %in% control_list

  # Subset data down to a single block 
  data_loc_df <- subset(data_seurat, cells=meta_loc_df$Coordinate)
  
  # Find controls available in the current object
  control_loc_list <- intersect(control_list, meta_loc_df$Group)
  
  # Wrangle data from R to Python and run inferncvpy through reticulate
  run_infercnvpy_my(data_loc_df, meta_loc_df, ref_df, control_loc_list, block)

}
