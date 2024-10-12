
calculate_guide_entropy_my <- function(normalized_guide_counts){
  # normalized_guide_counts is the matrix data_seurat@assays$RNA@data
  # after normalizing the data by Seurat::NormalizeData(normalization.method = "RC", scale.factor = 1)
  
  # need to get rid of zeros in entropy calculation. TODO implement this cleaner
  temp_mat <- Matrix::t(normalized_guide_counts) %>% unname(force = TRUE)
  my_list <- apply(temp_mat, 1, as.list) # creates a list of lists
  my_list_2 <- lapply(my_list,unlist, recursive = FALSE) # flatten it to be a list of vectors
  guide_pi <- lapply(my_list_2, function(x) {x[x!=0]}) # remove the zero elements
  
  log2_pi <- lapply(guide_pi, log2)
  
  
  plp <- Map(`*`, guide_pi, log2_pi)
  entropy_vector <- lapply(plp, sum) %>% unlist
  return(entropy_vector*-1)
}
