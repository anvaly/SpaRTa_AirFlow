# Author: Andy Kavran Date: 16-June-2022
# Funcitons focused on scanorama data integration
# 

prep_scanorama_integration_my <- function(seurat_list, num_features = 5000, sct = FALSE, assay = "SCT", slot = "data"){
  #given a list of seurat objects, prepare them for scanorama
  hvg <- Seurat::SelectIntegrationFeatures(object.list = seurat_list, nfeatures = num_features)
  if(sct == TRUE){
    seurat_list <- Seurat::PrepSCTIntegration(object.list = seurat_list, anchor.features = hvg, sct.clip.range = NULL) # TODO decide if we want to clip values
    expression <- lapply(seurat_list, get_SCT_residuals_my)
    #since we use scaled data, this already subsets to the hvgs
    names(expression) <- names(seurat_list)
  }else{
    seurat_list <- Seurat::PrepSCTIntegration(object.list = seurat_list, anchor.features = hvg, assay = assay)
    #subset genes
    expression <- lapply(seurat_list, GetAssayData, slot = slot, assay = assay)
    expression <- lapply(expression, Matrix::t) # scanorama expects cells by genes
    expression <- lapply(expression, subset_hvg_my, hvg = hvg)
    names(expression) <- names(seurat_list)
  }
  return(expression)
}

get_SCT_residuals_my <- function(seurat_data){
  exprs <- Seurat::GetAssayData(seurat_data, slot = "scale.data", assay = "SCT")
  exprs <- t(exprs) # transpose to cells by genes
  return(as.matrix(exprs))
}

get_gene_names_my <- function(expr_matrix){
    gene_names <- colnames(expr_matrix) # genes should be columns
    return(gene_names)
}

subset_hvg_my <- function(data_seurat_mat, hvg){
  # subsets the seurat data matrix to only include genes previously identified as highly variable.
  subset_data <- data_seurat_mat[,hvg]
  return(subset_data)
} 