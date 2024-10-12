# Author: Anna Lyubetskaya. Date: 20-01-20


matrix_to_Seurat_my <- function(input_path){
  ## Read a 10X matrix and create a Seurat object
  
  # Read the counts matrix
  mat <- read.table(file = input_path, header=TRUE)
  
  # Create a Seurat object
  data_seurat <- Seurat::CreateSeuratObject(mat)
  
  data_seurat <- Seurat::RenameCells(data_seurat, new.names=gsub("\\.", "-", colnames(data_seurat)))
  
  return(data_seurat)
}


read_split_10X_matrix_my <- function(input_path, type="filtered"){
  ## Read the 10X information from TSV files to a matrix
  
  # Choose the filtered or raw matrix set
  matrix_path <- paste0(input_path, type, "_feature_bc_matrix/matrix.mtx.gz")
  feature_path <- paste0(input_path, type, "_feature_bc_matrix/features.tsv.gz")
  barcode_path <- paste0(input_path, type, "_feature_bc_matrix/barcodes.tsv.gz")

  # Read the counts matrix
  mat <- Matrix::readMM(file = matrix_path)
  
  # Read gene name information
  feature_names_df <- readr::read_delim(feature_path, delim="\t", col_names=FALSE)
  # Read barcode name information
  barcode_names_df <- readr::read_delim(barcode_path, delim="\t", col_names=FALSE)

  # Add column and rownames to the object
  colnames(mat) <- barcode_names_df$X1
  rownames(mat) <- toupper(feature_names_df$X1)
  
  # Split features into GEX and Guides
  feature_indx <- feature_split_by_type_my(feature_names_df$X3)
  
  # Split the matrix into GEX and CRISPR parts
  mat_gex <- mat[feature_indx$gex,]
  mat_crispr <- mat[feature_indx$crispr,]
  
  # Write Guide matrix
  readr::write_delim(tibble::as_tibble(t(as.matrix(mat_crispr)), rownames="Barcode"), gsub(".mtx.gz", "_crispr.txt", matrix_path), delim="\t")
  
  # Write GEX and CRISPR matrices to file
  write.table(mat_gex, file=gsub(".mtx.gz", "_gex.mtx", matrix_path), row.names=TRUE, col.names=TRUE)
  write.table(mat_crispr, file=gsub(".mtx.gz", "_crispr.mtx", matrix_path), row.names=TRUE, col.names=TRUE)
}


feature_split_by_type_my <- function(feature_list){
  ## Find the indexes of GEX and non-GEX features

  return(list("gex" = which(feature_list == "Gene Expression"),
              "crispr" = which(feature_list != "Gene Expression")))
}

