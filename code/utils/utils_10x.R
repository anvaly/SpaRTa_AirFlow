# Author: Anna Lyubetskaya. Date: 20-02-03
# Functions to do with analyzing 10X data.


environment_10x_my <- function(){
  ## Check necessary package versions
  ## If install.packages doesn't work, try grabbing them from the BRAN website: http://bran.pri.bms.com/
  
  # You need Seurat 3.2 or newer
  if(packageVersion("Seurat") < '3.2'){
    install.packages("Seurat")
    devtools::install_github("satijalab/seurat", ref = "spatial")
  }
}
