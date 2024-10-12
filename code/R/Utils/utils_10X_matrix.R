# Author: Anna Lyubetskaya. Date: 20-08-05


source("code/utils/utils_tibble.R")


identify_mt_ribo_genes_my <- function(gene_list){
  ## Identify MT and ribo genes
  
  return(grep("^MT-|^RP[SL]", gene_list))
}


identify_unannotated_genes_my <- function(symbol_list){
  ## Find genes that have no GO annotation and standardizes symbol names
  
  if(!"mygene" %in% rownames(installed.packages())){
    BiocManager::install("mygene")
  }

  # Use mygene library to find GO annotations for all genes
  gene_annotation_res <- mygene::queryMany(symbol_list, scopes="symbol", fields=c("go"), species=ref_short)
  
  # Identify genes that have no GO annotation at all
  gene_unannotated_list <- gene_annotation_res$query[which(unlist(sapply(1:length(gene_annotation_res$go.BP), 
                                                                         function(x) is.null(gene_annotation_res$go.BP[[x]]) && 
                                                                           is.null(gene_annotation_res$go.CC[[x]]) && 
                                                                           is.null(gene_annotation_res$go.MF[[x]])
  )))]
  
  # Identify serial symbol names
  exclude_list <- gene_unannotated_list[which(grepl("RGD\\d\\d+\\.*\\d+|LOC\\d\\d+\\.*\\d+|AABR\\d\\d+\\.*\\d+|NEWGENE-\\d\\d+\\.*\\d+|LINC\\d\\d+\\.*\\d+|A[CFJLP]\\d\\d+\\.*\\d+", 
                                                    gene_unannotated_list))]
  
  return(exclude_list)
}


calculate_mt_ribo_my <- function(data_seurat){
  ## Calculate mitochondrial and ribosomal content
  
  # Calculate MT content
  data_seurat <- Seurat::PercentageFeatureSet(data_seurat, "^MT-", col.name = "mito_percent")
  
  # Calculate Ribo content
  data_seurat <- Seurat::PercentageFeatureSet(data_seurat, "^RP[SL]", col.name = "ribo_percent")
  
  
  return(data_seurat)
}


seurat_select_abundant_genes_my <- function(data_seurat, sct_threshold=0.5, spot_threshold=10, assay="SCT", slot="data", split_by=NULL){
  ## Find all genes that are well represented in the dataset  
  
  # Define subsets of barcodes to find abundant genes in each of them
  barcode_list <- list()
  if(!is.null(split_by)){
    # Meta data values by which to split the barcodes
    split_vals <- unique(data_seurat@meta.data[[split_by]])
    
    # Identify a list of barcodes for each meta data value
    for(v in split_vals){
      barcode_list[[v]] <- data_seurat@meta.data %>%
        dplyr::filter(!!rlang::sym(split_by) == v) %>%
        rownames()
    }
  } else{
    if(length(rownames(data_seurat@meta.data)) <= 100000){
      barcode_list[["All"]] <- rownames(data_seurat@meta.data)
    } else{
      barcode_list[["All"]] <- sample(rownames(data_seurat@meta.data), 100000, replace=FALSE)
    }
  }
  
  # Go through every set of barcodes if defined and find abundant genes
  gene_list <- list()
  for(v in names(barcode_list)){
    data_seurat_subset <- subset(data_seurat, cells=barcode_list[[v]])
    data_matrix <- Seurat::GetAssayData(data_seurat_subset, assay=assay, slot=slot)
    gene_list[[v]] <- rownames(data_matrix)[which(Matrix::rowSums(as.matrix(data_matrix) >= sct_threshold) >= spot_threshold)]
  }
  
  return(Reduce(intersect, gene_list))
}


seurat_gene_properties_my <- function(data_seurat, gene_list, assay="SCT", sct_threshold=1){
  ## Extract expression matrix, transform it into a long tibble, and calculate z-score for each gene across all barcodes
  
  # Extract appropriate data matrix, the list of genes, and the number of spots
  exprs <- Seurat::GetAssay(data_seurat, assay = assay)
  gene_symbols <- rownames(exprs)
  spot_num <- length(colnames(exprs))
  
  # Match genes in signature to rownames of sparse matrix
  keep_gene_idx <- match(gene_list, gene_symbols)
  exprs <- exprs[keep_gene_idx,]
  
  # Calculate stats
  expression_df <- tibble::as_tibble(as.matrix(exprs), rownames = "Symbol") %>%
    df_wide2long_my(key="Coordinate", val="Expression") %>%
    dplyr::group_by(Symbol) %>%
    dplyr::summarise(ExpressionMean = round(mean(Expression), 3),
                     ExpressionMedian = round(median(Expression), 3),
                     ExpressionSD = round(sd(Expression), 3),
                     ExpressionMax = round(max(Expression), 3),
                     CountSpots = sum(Expression >= sct_threshold),
                     PercentSpots = round(CountSpots / spot_num * 100)) %>%
    dplyr::arrange(ExpressionMean)

  return(expression_df)
}


seurat_expression_to_long_tibble_my <- function(data_seurat, assay="SCT", slot="scale.data", gene_list=NULL){
  ## Extract expression matrix, transform it into a long tibble, and calculate z-score for each gene across all barcodes
  
  if(is.null(gene_list)){
    gene_list <- rownames(data_seurat)
  }
  
  expression_df <- tibble::as_tibble(as.matrix(Seurat::GetAssayData(data_seurat, assay=assay, slot=slot)[gene_list,]), rownames="Symbol") %>%
    df_wide2long_my(key="Coordinate", val=paste0(assay, "_", slot))

  return(expression_df)
}


pseudo_bulk_counts_my <- function(data_seurat, log_transform=TRUE, assay="Spatial", slot="data"){
  # Calculate pseudo bulk counts a Seurat object
  # Counts represented as a log2 value with pseudo-count = 1

  # Extract the counts matrix
  counts_matrix <- as.matrix(Seurat::GetAssayData(data_seurat, assay=assay, slot=slot))
  
  # Make sure that there are no negative values in the matrix
  counts_matrix_shifted <- counts_matrix - min(counts_matrix)
  
  # Calculate row sums
  counts_vector <- rowSums(counts_matrix_shifted)
  
  # Perform log2 transformation with a pseudocount = 1
  if(log_transform == TRUE){
    counts_vector <- round(log2(counts_vector + 1), 5)
  }

  return(as.list(counts_vector))
}


occupancy_counts_my <- function(data_seurat, log_transform=TRUE, assay="RNA", slot="data"){
  # Calculate pseudo bulk counts a Seurat object
  # Counts represented as a log2 value with pseudo-count = 1
  
  # Extract the counts matrix and calculate row sums
  counts_vector <- rowSums(as.matrix(Seurat::GetAssayData(data_seurat, assay=assay, slot=slot)) != 0)
  
  # Perform log2 transformation with a pseudocount = 1
  if(log_transform == TRUE){
    counts_vector <- round(log2(counts_vector + 1), 4)
  }
  
  return(as.list(counts_vector))
}


pseudo_bulk_vectors_combine_my <- function(pseudo_count_list, sample_names){
  # Combine a named list of lists into an expression wide tibble
  
  # Create a tibble of pseudo counts from a list of named lists; substitute NA for 0
  data_df <- list2tibble_my(pseudo_count_list, rownames="Symbol", transpose=TRUE)
  
  for(i in 1:ncol(data_df)){
    if(length(which(is.na(data_df[[i]]))) > 0){
      data_df[[i]] <- as.numeric(data_df[[i]])
      data_df[which(is.na(data_df[[i]])), i] <- 0
    }
  }
  
  return(data_df)
}
