# Author: Anna Lyubetskaya. Date: 22-01-10
# Run inferCNV

# Manual: https://github.com/broadinstitute/inferCNV/wiki
# Vignette: https://bioconductor.org/packages/devel/bioc/vignettes/infercnv/inst/doc/inferCNV.html

# Some standard input files are here: data/import/inferCNV/hg38_gencode_v27.txt


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/R/Utils/utils_10X_in_out.R")
library(infercnv)


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

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, "/", cohort_name, "/")

dir.create(output_path_init, showWarnings=FALSE)
dir.create(output_path, showWarnings=FALSE)


## INGEST DATA ----


# Ingest a set of RDS Seurat objects
data_seurat <- readRDS(paste0(input_path, cohort_name, ".rds"))


## WRANGLE inferCNV INPUTS----


# Extract meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")
# Find the max value across reference pathology classes
meta_df["inferCNVcontrol"] <- matrixStats::rowMaxs(as.matrix(meta_df[paste0("Pathology.", groups_ref, ".percent")]))

# Establishing the control and target groups
meta_df <- meta_df %>%
  dplyr::mutate(Group = ifelse(Pathology.Group %in% groups_ref & inferCNVcontrol >= 80, 
                               user.Sample_Name,  # paste0(Pathology.Group, "_", user.Sample_Name), 
                               paste0("malignant_", user.Sample_Name))) %>%
  dplyr::select(user.Block_ID, Pathology.Group, Coordinate, Group)

# Control groups
control_list <- meta_df$Group[which(!grepl("malignant", meta_df$Group))]

table(meta_df$Group)


for(block in unique(meta_df$user.Block_ID)){
  
  gc()
  
  # Subset data down to manage the size
  meta_loc_df <- meta_df %>%
    dplyr::filter(user.Block_ID == block | Group %in% control_list)
  
  # Write sample - group assignments to file
  readr::write_delim(meta_loc_df %>%
                       dplyr::select(Coordinate, Group), 
                     paste0(output_path, block, "_samples.txt"), 
                     delim="\t", col_names=FALSE)
  
  
  # Extract the raw counts matrix from the Seurat object
  data_matrix <- as.data.frame(Seurat::GetAssayData(subset(data_seurat, cells=meta_loc_df$Coordinate), assay="Spatial"))
  write.table(data_matrix, paste0(output_path, block, "_matrix.txt"), sep="\t", quote=FALSE)
  
  table(colnames(data_matrix) == meta_loc_df$Coordinate)
  
  
  # Make sure that the input gene order file has upper case symbols
  order_df <- readr::read_delim(inferCNV_file, delim="\t", col_names=FALSE)
  readr::write_delim(order_df %>%
                       dplyr::mutate(X1 = toupper(X1)), inferCNV_file, delim="\t", col_names=FALSE)
  
  
  ## CREATE inferCNV OBJECT ----
  
  
  # Create an inferCNV object from previously generated data
  infercnv_obj <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = paste0(output_path, block, "_matrix.txt"),
    annotations_file = paste0(output_path, block, "_samples.txt"),
    delim = "\t",
    gene_order_file = inferCNV_file,
    ref_group_names = control_list)
  
  
  ## RUN inferCNV ----
  
  
  infercnv_obj_default <- infercnv::run(
    infercnv_obj,
    out_dir=paste0(output_path, block),
    cutoff=0.1,  # from manual: cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    plot_steps=TRUE,
    denoise=TRUE,
    HMM=FALSE,
    #HMM_type="i3",
    cluster_by_groups=TRUE, 
    no_prelim_plot=FALSE,
    analysis_mode="subclusters",
    png_res=300,
    num_threads=8  # parallel::detectCores() - 1
  )
  
}
