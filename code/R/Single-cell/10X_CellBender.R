# Author: Anna Lyubetskaya. Date: 21-02-20
# Re-run 10X filtering procedure using CellBender

# Run Cell Bender on 10X raw output using a shell script that invokes a Stash version of Python + CellBinder via command line
# Appropriate conda, dependencies, and CellBender are installed here /XXXX

########################################################################################################
#
# To install conda and CellBender, run this in TERMINAL
#
# Install minconda: https://waylonwalker.com/install-miniconda/
## cd /XXXX
## mkdir -p miniconda3
## wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3/miniconda.sh
## bash miniconda3/miniconda.sh -b -u -p miniconda3
## miniconda3/bin/conda init bash
#
# Setup CellBender: https://cellbender.readthedocs.io/en/latest/installation/index.html
## bash
## conda create -n CellBender python=3.7
## source activate CellBender
## conda install -c anaconda pytables
## conda install pytorch torchvision -c pytorch
## cd /XXXX
## git clone https://github.com/broadinstitute/CellBender.git
## pip install -e CellBender
#
########################################################################################################
#
# To test installation, run this in TERMINAL
## cd /XXXX
## source activate CellBender
## cellbender remove-background
#
########################################################################################################


## SETUP PYTHON-RETICULATE ENVIRONMENT ----


library(reticulate)
library(foreach)

Sys.setenv(RETICULATE_PYTHON="XXXX")

use_condaenv("CellBender", 
             "XXXX",
             required = TRUE)

cb <- import("cellbender")


## ENVIRONMENT ----


# Allow piping throughout the package
`%>%` <- magrittr::`%>%`


## PATHS ----


# Input location containing 10X folders
input_path <- "XXXX"

# Input file with sample meta data
meta_file <- paste0(input_path, "meta_data_snRNA.txt")



## INGEST DATA ----


# Read sample file
meta_df <- readr::read_delim(file=meta_file, delim="\t") %>%
  tidyr::drop_na()


## RUN CELLBENDER ----


# Parallelization didn't work for some reason
# cl <- parallel::makeCluster(2)
# doParallel::registerDoParallel(cl)
# foreach::foreach(i = 12:15) %dopar% {  # 1:nrow(meta_df)
# };
# parallel::stopCluster(cl)


for(i in 1:nrow(meta_df)){
  # Input data path
  filename_in <- paste0(meta_df[[i, "FullPath"]], "raw_feature_bc_matrix.h5")
  # Output data path
  filename_out <- paste0(meta_df[[i, "FullPath"]], "cellbender_feature_bc_matrix.h5")
  
  # CellBender parameters
  realcell_num <- meta_df[[i, "CellBender_RealCells"]]
  droplet_num <- meta_df[[i, "CellBender_Droplets"]]
  
  # Run CellBender
  command <- paste("cellbender remove-background", 
                   "--input", filename_in, 
                   "--output", filename_out, 
                   "--expected-cells", realcell_num, 
                   "--total-droplets-included", droplet_num, 
                   "--fpr 0.01 --epochs 150")
  
  if(!file.exists(filename_out)){
    system(command)
  }
}
