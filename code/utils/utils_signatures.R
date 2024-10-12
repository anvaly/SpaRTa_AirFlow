# Author: Anna Lyubetskaya. Date: 20-05-01


read_signatures_my <- function(sig_path){
  ## Read a standard set of signatures
  
  # Load all signatures
  signature_df <- readr::read_delim(file=sig_path, delim="\t")
  
  # Create a named list of signatures
  signature_list <- sapply(signature_df$Gene_list, function(x) stringr::str_split(toupper(x), ","))
  names(signature_list) <- signature_df$Signature_name
  
  return(signature_list)
}


read_filter_signatures_my <- function(sig_path, gene_list=NULL, sig_names=NULL, sig_length_min=1, sig_length_max=1000, ratio_threshold=0){
  ## Read a standard set of signatures and filter them
  
  # Load all signatures
  signature_init_df <- readr::read_delim(file=sig_path, delim="\t")
  
  # Create a named list of signatures
  signature_init_list <- sapply(signature_init_df$Gene_list, function(x) stringr::str_split(toupper(x), ","))
  names(signature_init_list) <- signature_init_df$Signature_name
  
  # Filter signatures
  signature_list <- filter_signature_set_my(signature_init_list, gene_list=gene_list, sig_names=sig_names, 
                                            sig_length_min=sig_length_min, sig_length_max=sig_length_max, ratio_threshold=ratio_threshold)
  
  return(signature_list)
}


read_filter_convert_signatures_my <- function(sig_path, gene_list=NULL, sig_names=NULL, sig_length_min=1, sig_length_max=1000, ratio_threshold=0, match_table_path){
  ## Read a standard set of signatures of human genes, convert to mouse homolog, and filter them
  ## Skip the signatures starting with "Syng", assuming they are already mouse signatures.
  ## Only works from human to mouse
  ## match_table_file <- "XXXX"
  
  # Load all signatures
  signature_init_df <- readr::read_delim(file=sig_path, delim="\t")
  
  # Create a named list of signatures
  signature_init_list <- sapply(signature_init_df$Gene_list, function(x) stringr::str_split(toupper(x), ","))
  names(signature_init_list) <- signature_init_df$Signature_name
  
  # Load match table
  match_table <- readr::read_delim(file=match_table_path, delim="\t")
  match_table <- match_table[!is.na(match_table$Symbol_human),]
  
  # Convert human gene symbols to mouse homolog
  hash_h2m <- hash::hash(keys = match_table$Symbol_human, values = match_table$Symbol_mouse)
  for(i in 1:length(signature_init_list)){
    if(!stringr::str_detect(names(signature_init_list)[i], pattern = "Syng")){
      signature_init_list[[i]] <- signature_init_list[[i]][signature_init_list[[i]] %in% hash::keys(hash_h2m)]
      signature_init_list[[i]] <- toupper(hash::values(hash_h2m, keys=signature_init_list[[i]]))
    }
    
  }
  
  # Filter signatures
  signature_list <- filter_signature_set_my(signature_init_list, gene_list=gene_list, sig_names=sig_names, 
                                            sig_length_min=sig_length_min, sig_length_max=sig_length_max, ratio_threshold=ratio_threshold)
  
  return(signature_list)
}


filter_signature_set_my <- function(signature_list, gene_list=NULL, sig_names=NULL, sig_length_min=1, sig_length_max=250, ratio_threshold=50){
  ## Remove all gene signatures that have no user-defined genes
  ## This function is specific to an expression dataset
  
  # Filter down to only signatures of interest by name
  if(!is.null(sig_names)){
    signature_list <- signature_list[sig_names]
  }
  
  # Intersect every signature with a user-defined list of genes
  if(!is.null(gene_list)){
    signature_filt_list <- lapply(signature_list, function(x) intersect(x, gene_list))
  } else{
    signature_filt_list <- signature_list
  }
  
  # Find lengths of filtered genes
  signature_lengths_before <- sapply(signature_list, function(x) length(x))
  # Find lengths of filtered genes
  signature_lengths_after <- sapply(signature_filt_list, function(x) length(x))
  
  # Calculate the difference in signature lengths before and after
  ratio_list <- round(signature_lengths_after / signature_lengths_before * 100)
  
  # Remove signatures that lost more than X% of their genes
  signature_filt_list <- signature_filt_list[ratio_list >= ratio_threshold &
                                               signature_lengths_before <= sig_length_max &
                                               signature_lengths_after >= sig_length_min]
  
  cat("Signatures filtered. Before = ", length(signature_list), ". After = ", length(signature_filt_list), "\n")
  
  return(signature_filt_list)
}


sig_lengths_df_my <- function(signature_list){
  ## Create a tibble of signature lenths
  
  # Signature length
  sig_lengths <- sapply(names(signature_list), function(x) length(signature_list[[x]]))
  
  # Create a tibble with signature lengths
  sig_len_df <- tibble::tibble("Signature_name" = names(sig_lengths),
                               "sig_len" = unname(sig_lengths)) %>%
    dplyr::arrange(Signature_name)

  return(sig_len_df)  
}


annotate_signature_pairs_my <- function(signature_list){
  ## For a list of signatures, find the intersection between each pair of signatures
  ## This function is NOT specific to an expression dataset
  ## Run: signatures_list <- collect_signatures_my(sig_names=NULL, in_path="C:/Users/lyubetsa/Documents/Data/Signatures/")
  
  # Signature length
  sig_lengths <- sapply(names(signature_list), function(x) length(signature_list[[x]]))
  
  # Find every 2 element combination of signature names
  sig_pairs <- combn(names(sig_lengths), 2)
  
  # Create a tibble with intersection values between two signatures
  sig_intersect_df <- tibble::tibble("Signature_name" = sig_pairs[1,],
                                     "Signature_other" = sig_pairs[2,],
                                     "Intersect_num" = 0,
                                     "Intersect_percent_smaller" = 0,
                                     "Sig_length" = 0,
                                     "Sig_length_other" = 0)
  
  # Calculate intersection between every pair of signatures
  for(j in 1:ncol(sig_pairs)){
    sig_n1 <- sig_pairs[1,j]
    sig_n2 <- sig_pairs[2,j]
    
    # Number of genes in both signatures
    value <- length(intersect(signature_list[[sig_n1]], signature_list[[sig_n2]]))
    sig_intersect_df[[j, "Intersect_num"]] <- value
    
    # Number of genes in both signatures divided by the length of the smaller signature
    percent <- round(value / min(sig_lengths[[sig_n1]], sig_lengths[[sig_n2]]) * 100)
    sig_intersect_df[[j, "Intersect_percent_smaller"]] <- percent
    
    # Signature lengths
    sig_intersect_df[[j, "Sig_length"]] <- sig_lengths[[sig_n1]]
    sig_intersect_df[[j, "Sig_length_other"]] <- sig_lengths[[sig_n2]]
  }

  sig_intersect_sym_df <- sig_intersect_df %>%
    dplyr::rename(Signature_temp = Signature_name) %>%
    dplyr::rename(Signature_name = Signature_other, Signature_other = Signature_temp)
  
  return(rbind(sig_intersect_df, sig_intersect_sym_df))
}


calculate_signature_scores_my <- function(signatures, eset_object){
  ## Calculate signatures scores from a named list of signatures and a normalized eSet object
  
  library(bmsgsa2)
  
  # Create a signature object
  signature_object <- bmsgsa2::makeGeneList(signatures)   # identifier_type=rep("ensembl_v86", length(signatures)
  
  # Perform gene signature analysis
  signature_scores <- bmsgsa2::scoreList(signature_object, eset_object, algorithm="median_z_score")
  
  # Make a signature score tibble
  signatures_df <- tibble::as_tibble(signature_scores@scores, rownames="id") %>%
    dplyr::rename(DepMap_ID = id)
  
  return(signatures_df)
}


collect_msig_my <- function(species="Homo sapiens"){
  ## Use msigdbr package to collect MSigDB gene lists

  m_init_df <- msigdbr::msigdbr(species=species)
    
  m_df <- m_init_df %>%
    dplyr::group_by(gs_name) %>%
    dplyr::summarise(Gene_list = paste0(sort(unique(gene_symbol)), collapse=",")) %>%
    dplyr::rename(Signature_name = gs_name)

  # in_path <- "C:/Users/lyubetsa/Documents/Data/Signatures/"
  # readr::write_delim(m_df, path = paste0(in_path, "/signatures_msigdb_20200802.txt"), delim = "\t")
  
  # Create a named list of signatures
  signature_list <- sapply(m_df$Gene_list, function(x) stringr::str_split(x, ","))
  names(signature_list) <- m_df$Signature_name
  
  return(signature_list)  
}


invert_list_my <- function(init_list){
  ## From a list of lists, invert names and values to a tibble
  
  value_list <- unique(unlist(unname(init_list)))
  
  name_list <- unlist(unname(sapply(value_list, 
                                    function(x) paste(sort(unique(names(init_list[which(sapply(unname(init_list), 
                                                                                               function(y) x %in% y))]))), collapse=";")
  )))
  
  df <- tibble::tibble("Symbol" = value_list,
                       "Sig_Name" = name_list,
                       "Num_Sigs" = stringr::str_count(name_list, ";") + 1) %>%
    dplyr::mutate(InSignature = Num_Sigs > 0)
  
  return(df)
}
