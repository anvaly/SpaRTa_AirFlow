# Author: Anna Lyubetskaya. Date: 23-03-14
# Run cell2location: https://github.com/BayraktarLab/cell2location


## SETUP PYTHON ----


# Full installation of cell2location

#### Setup environment
## cd /XXXX
## echo $HOME
## HOME=/XXXX
## echo $HOME

#### Install miniconda
## wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
## bash Miniconda3-latest-Linux-x86_64.sh
## export PYTHONNOUSERSITE="literallyanyletters"

#### Setup environment and install cell2location
## miniconda3/bin/conda create -y -n cell2loc_env python=3.9
## source ~/miniconda3/etc/profile.d/conda.sh
## conda activate cell2loc_env
## pip install tensorflow -t XXXX --upgrade
## pip install keras -t XXXX --upgrade
## pip install cell2location[tutorials] -t XXXX --upgrade
## conda activate cell2loc_env
## python -m ipykernel install --user --name=cell2loc_env --display-name='Environment (cell2loc_env)'

#### Add the path to sys
## python
## import sys
## sys.path.insert(0, '/XXXX')


## SETUP ENVIRONMENT ----


# Solution works with reticulate version 1.22
# See this: https://github.com/rstudio/reticulate/issues/1176
if(packageVersion("reticulate") != "1.22"){
  remotes::install_version("reticulate", "1.22")
}


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

Sys.setenv(RETICULATE_PYTHON="XXXX")
Sys.getenv("RETICULATE_PYTHON")

# Configure reticulate with my own Python on Stash
library(reticulate)
packageVersion("reticulate")

use_condaenv("cell2loc_env", 
             "XXXX",
             required=TRUE)

# # Configure tensorflow / keras - https://www.reddit.com/r/tensorflow/comments/mte98s/r_error_in_py_call_implcallable_dotsargs/
# if(FALSE){
# 
#   install.packages("tensorflow")
#   library(tensorflow)
#   tensorflow::install_tensorflow()
#   
#   install.packages("keras")
#   library(keras)
#   keras::install_keras()
#   
#   system("XXXX")
#   
#   reticulate::py_config()
# }
# 
# library(tensorflow)
# library(keras)

# Import python packages
py_sys <- import("sys")
py_sys$path[1] <- "XXXX"

py_tf <- import("tensorflow")
py_tf <- import("keras")
py_sc <- import("scanpy")
py_np <- import("numpy")
py_ad <- import("anndata")
py_c2l <- import("cell2location")

# import matplotlib.pyplot as plt
# import matplotlib as mpl


## PARAMETERS ----


# Cohort name for outputs
cohort_name <- "PDAC_DonorA_s6_annotated"


## PATHS ----


# Cell type single cell reference
sc_ref <- "XXXX"

# File with unspecific probes in the Symbol column
unspecific_probes_file <- "XXXX"

# Input folder
input_path <- "XXXX"

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, "/", cohort_name, "/")

dir.create(output_path_init, showWarnings=FALSE)
dir.create(output_path, showWarnings=FALSE)


## INGEST DATA ----


# Read unspecific probes if defined
if(!is.null(unspecific_probes_file)){
  unspecific_probes_list <- readr::read_delim(unspecific_probes_file, delim="\t")
}


## Ingest RDS objects

# Ingest a set of RDS Seurat objects: ST
data_seurat_st <- readRDS(paste0(input_path, cohort_name, ".rds"))

# Ingest a set of RDS Seurat objects: SC
data_seurat_sc <- readRDS(sc_ref)

# Gene list
gene_list <- intersect(rownames(data_seurat_st@assays$Spatial)[which(!grepl("^MT-|^RP[SL]", 
                                                                            rownames(data_seurat_st@assays$Spatial)))],
                       rownames(data_seurat_sc@assays$RNA))
# Exclude unspecific probes
gene_list <- setdiff(gene_list, unspecific_probes_list$Symbol)

# Subset object to the common list of genes
data_seurat_st <- subset(data_seurat_st, features=gene_list)
data_seurat_sc <- subset(data_seurat_sc, features=gene_list)

# Extract meta data: ST
meta_st_df <- tibble::as_tibble(data_seurat_st@meta.data, rownames="Coordinate")

# Extract meta data: SC
meta_sc_df <- tibble::as_tibble(data_seurat_sc@meta.data, rownames="Coordinate")


## Transform RDS objects to Python AnnData format


unique(meta_sc_df$cell_type)

# ST: Create a Python AnnData object: X, obs, var
py_adata_st <- py_ad$AnnData(Seurat::GetAssayData(data_seurat_st, assay="Spatial", slot="data"),
                             rownames(Seurat::GetAssayData(data_seurat_st, assay="Spatial", slot="data")),
                             meta_st_df %>% tibble::column_to_rownames("Coordinate")
)

# ST: Transpose the object
py_adata_st <- py_ad$AnnData$transpose(py_adata_st)

# SC: Create a Python AnnData object: X, obs, var
py_adata_sc <- py_ad$AnnData(Seurat::GetAssayData(data_seurat_sc, assay="RNA", slot="data"),
                             rownames(Seurat::GetAssayData(data_seurat_sc, assay="RNA", slot="data")),
                             meta_sc_df %>% tibble::column_to_rownames("Coordinate")
)

# SC: Transpose the object
py_adata_sc <- py_ad$AnnData$transpose(py_adata_sc)


## PREPARE REFERENCE DATA ----


# Check reference gene abundances
genes_selected <- py_c2l$utils$filtering$filter_genes(py_adata_sc, cell_count_cutoff=5, 
                                                      cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# Prepare anndata for the regression model
py_c2l$models$RegressionModel$setup_anndata(adata=py_adata_sc,
                                            # 10X reaction / sample / batch
                                            batch_key="orig.ident",
                                            # cell type, covariate used for constructing signatures
                                            labels_key="cell_type")

# Create the regression model
mod <- py_c2l$models$RegressionModel(py_adata_sc)

# View anndata_setup as a sanity check
mod$view_anndata_setup()

# Train the model to estimate the reference cell type signatures.
mod$train(max_epochs=250, use_gpu=FALSE)
