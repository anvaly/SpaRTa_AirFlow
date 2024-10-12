# Author: Anna Lyubetskaya. Date: 21-03-05
# Following: https://theislab.github.io/scanpy-in-R/


install.packages("renv")

renv::init()

renv::install("reticulate")

renv::use_python()

pkgs <- c(
  "renv",
  "reticulate",
  "png",
  "ggplot2",
  "BiocManager",
  "Seurat"
)

bioc_pkgs <- c(
  "SingleCellExperiment",
  "scater",
  "multtest"
)

# If you are using an {renv} environment
renv::install(pkgs)

# Otherwise do it the normal way
install.packages(pkgs)

# Install Bioconductor packages
BiocManager::install(bioc_pkgs, update = FALSE)

py_pkgs <- c(
  "scanpy",
  "python-igraph",
  "louvain"
)

reticulate::py_install(py_pkgs)

renv::snapshot()


## TEST ----


## Call this function from R via reticulate...
# https://scanpy.readthedocs.io/en/stable/api/scanpy.read_10x_h5.html
scanpy.read_10x_h5("XXXX")
