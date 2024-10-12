# Author: Anna Lyubetskaya. Date: 21-06-15
# Run decontX using raw and filtered ST data to understand contamination levels (ambient signal)

# celda Bioconductor - http://bioconductor.org/packages/release/bioc/html/celda.html
# celda Git - https://github.com/campbio/celda
# decontX - https://rdrr.io/bioc/celda/f/vignettes/decontX.Rmd

# Installation
# Sys.getenv('R_HOME')
# cat /etc/issue
# deb https://cloud.r-project.org/bin/linux/ubuntu hirsute-cran40
# cd /opt/tbio/domino_202102/binaries/
# sudo apt-get update
# sudo apt-get install r-base


## ENVIRONMENT ----


# Allow piping throughout the package
`%>%` <- magrittr::`%>%`

if(packageVersion("celda") < 1.8){
  install.packages("code/R/Utils/celda_1.8.1.tar.gz", repos = NULL, type="source")
}

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_vis.R")


## PARAMETERS ----


cohort_name <- "HumanPanc"


## PATHS ----


# Input location containing 10X folders
input_path <- "XXXX"

# Input file with sample meta data
meta_file <- paste0(input_path, "meta_data.txt")

# Output folder
output_path_init <- "XXXX"
output_models <- paste0(output_path_init, "DecontX_models/")
output_figs <- paste0(output_path_init, "Figures/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_models, showWarnings = FALSE)
dir.create(output_figs, showWarnings = FALSE)


## INGEST DATA ----


# Read sample file
cohort_meta_df <- readr::read_delim(file=meta_file, delim="\t")  %>%
  dplyr::filter(grepl(cohort_name, Sample_Name))# %>% tidyr::drop_na()


## RAW DATA CHECK ----


for(i in 1:nrow(cohort_meta_df)){
  
  
  ## INGEST FILTERED 10X DATA ----
  
  
  # Input data path
  filename_in <- cohort_meta_df[[i, "FullPath"]]
  # Sample name
  sample_name <- cohort_meta_df[[i, "Sample_Name"]]
  
  # Create a sample-specific output folder
  #output_path <- paste0(output_path_init, sample_name, "/")
  #dir.create(output_path, showWarnings = FALSE)
  
  # Ingest filtered data
  data_seurat <- read_10X_spatial_folder_my(filename_in)


  ## INGEST RAW 10X DATA AND RUN DECONTX ----
  
  
  rds_filename <- paste0(output_models, sample_name, "_decontx.RDS")
  
  if(!file.exists(rds_filename)){
    
    # Ingest raw (unfiltered) data
    data_seurat_raw <- read_10X_spatial_folder_my(filename_in, filename="raw_feature_bc_matrix.h5")
    
    # 10X spot list not under tissue
    spot_bad_list <- setdiff(colnames(data_seurat_raw), colnames(data_seurat))
    
    # Genes appearing in both raw and filtered datasets
    gene_list <- intersect(rownames(data_seurat), rownames(data_seurat_raw))
    
    # Create data matrices for spots under and outside of the tissue
    data_matrix <- as.matrix(Seurat::GetAssayData(data_seurat, assay="Spatial", slot="data"))[gene_list, ]
    bkgr_matrix <- as.matrix(Seurat::GetAssayData(data_seurat_raw, assay="Spatial", slot="data"))[gene_list, spot_bad_list]

    print(sample_name)
    print(dim(data_matrix))
    print(dim(bkgr_matrix))
    
    # Use the filtered data as the primary dataset and raw data as the backgound - decontX will automatically find the spots present in one vs the other
    logfile <- paste0(output_models, sample_name, "_decontx.log")
    decontx_res <- celda::decontX(x = data_matrix, background = bkgr_matrix, varGenes = 10000, #z = rep(1, length(colnames(data_matrix))),
                                  logfile = logfile, verbose = TRUE)
    
    saveRDS(decontx_res, rds_filename)
  } else{
    decontx_res <- readRDS(rds_filename)
  }
  
  
  ## INVESTIGATE DECONTX OUTPUT ----
  
  
  # Add contaimination level to meta data
  data_seurat@meta.data["Contamination"] <- round(decontx_res$contamination * 100, 1)
  # Add contaimination level to meta data
  data_seurat@meta.data["nCount_Spatial_decontX"] <- round(unname(colSums(as.matrix(decontx_res$decontXcounts)), 1))
  # Add contaimination level to meta data
  data_seurat@meta.data["nFeature_Spatial_decontX"] <- round(unname(colSums(as.matrix(decontx_res$decontXcounts) > 0), 1))
  
  # Contamination spatial plot
  p1 <- spatial_feature_plot_my(data_seurat, feature="Contamination", min.cutoff="q0", max.cutoff="q100", 
                                name=paste0(sample_name, "\nContamination percent"))
  
  # UMI count spatial plot
  p2 <- spatial_feature_plot_my(data_seurat, feature="nCount_Spatial", min.cutoff="q1", max.cutoff="q99", 
                                name=paste0(sample_name, "\nUMIs per spot BEFORE decontX"))
  
  # Feature count spatial plot
  p3 <- spatial_feature_plot_my(data_seurat, feature="nFeature_Spatial", min.cutoff="q1", max.cutoff="q99", 
                                name=paste0(sample_name, "\nFeatures per spot BEFORE decontX"))

  # UMI count spatial plot
  p4 <- spatial_feature_plot_my(data_seurat, feature="nCount_Spatial_decontX", min.cutoff="q1", max.cutoff="q99", 
                                name=paste0(sample_name, "\nUMIs per spot AFTER decontX"))
  
  # Feature count spatial plot
  p5 <- spatial_feature_plot_my(data_seurat, feature="nFeature_Spatial_decontX", min.cutoff="q1", max.cutoff="q99", 
                                name=paste0(sample_name, "\nFeatures per spot AFTER decontX"))
  
  # Write figure to file
  filename <- paste0(output_figs, sample_name, "_contamination_spatial")
  write_plot2file_my(patchwork::wrap_plots(list(p1, p2, p3, p4, p5), nrow=1), filename, num_col=5, num_row=1)
}
