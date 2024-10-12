# Author: Anna Lyubetskaya. Date: 20-01-06


edit_group_name_my <- function(df, column){
  # Edit group name
  
  df <- df %>%
    dplyr::mutate(Column = gsub("^\\s+|\\s+$", "", !!rlang::sym(column)),
                  Column = gsub(" ", "_", Column))  # Column = tolower(Column)
  
  # If the input column name doesn't match the new column name, perform a substitution
  if(column != "Column"){
    df <- df %>%
      dplyr::select(-dplyr::all_of(column))
  
    colnames(df) <- gsub("^Column$", column, colnames(df))
  }
  
  return(df)
}


edit_symbol_my <- function(df, column){
  # Edit gene symbol
  
  df <- df %>%
    dplyr::mutate(Symbol = gsub("^\\s+|\\s+$", "", !!rlang::sym(column)),
                  Symbol = toupper(Symbol))
  
  # If the input column name doesn't match the new column name, perform a substitution
  if(column != "Symbol"){
    df <- df %>%
      dplyr::select(-dplyr::all_of(column))
    
    colnames(df) <- gsub("^Symbol$", column, colnames(df))
  }
  
  return(df)
}


check_reference_my <- function(gene_list, species="Rn"){
  ## Check a list of genes against a reference list
  ## Reference genome = "rat", "human", "mouse"
  
  # Folder with reference gene lists
  reference_path <- "data/import/References/"
  
  # Load reference gene list
  ref_df <- readr::read_delim(paste0(reference_path, species, ".txt"), delim="\t")
  
  # List of genes in the reference
  genes_included <- gene_list[which(gene_list %in% ref_df$Symbol)]
  genes_excluded <- setdiff(gene_list, genes_included)
  
  cat("Removed", length(genes_excluded), "genes of", length(gene_list), "genes total\n")
  # print("Excluding:")
  # print(sort(genes_excluded))
  
  return(list("genes_included" = genes_included,
              "genes_excluded" = genes_excluded))
}


check_limma_alias_my <- function(gene_list, species="Rn"){
  ## Attempt to find a primary gene symbol using limma alias2Symbol function
  ## species = "Hs", "Mm", "Rn"
  
  if(species == "Hs"){
    gene_list <- toupper(gene_list)
  } else{
    gene_list <- stringr::str_to_title(gene_list)
  }
  
  # Find a synonym for every gene in the list
  gene_primary_list <- sapply(gene_list, function(x) limma::alias2Symbol(x, species=species))
  # Index of genes that have exactly one gene synonym
  index <- which(sapply(names(gene_primary_list), 
                        function(x) length(gene_primary_list[[x]]) == 1 && x != gene_primary_list[[x]]))

  cat("Genes translated =", length(index), "\t")
  
  return(gene_primary_list[index]) 
}
