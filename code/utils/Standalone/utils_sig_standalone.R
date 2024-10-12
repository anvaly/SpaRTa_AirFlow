# Author: Anna Lyubetskaya. Date: 20-07-28


# Allow piping throughout the package
`%>%` <- magrittr::`%>%`
source("code/utils/utils_in_out.R")
source("code/utils/utils_signatures.R")


collect_signatures_my <- function(in_path=NULL){
  ## Collect a named list of signatures
  ## This function is NOT specific to an expression dataset
  ## Run: signatures_list <- collect_signatures_my(in_path="C:/Users/lyubetsa/Documents/Data/Signatures/")
  
  if(is.null(in_path)){
    in_path <- "XXXX"
  }
  
  # List of signature files
  paths <- paste0(in_path, c(#"signatures_bongene.txt",
                             #"signatures_celltypes_kumarn.txt",
                             #"signatures_discovery.txt",
                             #"signatures_emt_kumarn.txt",
                             #"signatures_general.txt",
                             #"signatures_listtracker_20200105.txt", 
                             #"signatures_myeloidGeneSet_suthrams_20200910.txt",
                             #"signatures_pancreas_brownea_20200723.txt",
                             #"signatures_pancreas_drokhle.txt",
                             #"signatures_prostate_drokhle.txt",
                             #"signatures_prostate_powlesr.txt",
                             #"signatures_siemersn.txt"
                             #"integrated_panc_20201011.txt"
                             "signatures_rshiny_20201127.txt"
                             ))
  
  # Load all signature files into a tibble
  signatures_df <- read_dir2file_my(paths)

  ## Validate a signature tibble
  signatures_df <- validate_signature_df_my(signatures_df) %>%
    dplyr::mutate(Gene_list = toupper(Gene_list))

  # Write standard signatures to file
  # Once this file is generated, it can be used an input for any analysis
  readr::write_delim(signatures_df, path = paste0(in_path, "/processed_signatures_standard.txt"), delim = "\t")
  
  # Create a named list of signatures
  signature_list <- sapply(signatures_df$Gene_list, function(x) stringr::str_split(x, ","))
  names(signature_list) <- signatures_df$Signature_name
  
  return(signature_list)
}


validate_signature_df_my <- function(signatures_df){
  ## Validate a signature tibble
  
  # Make sure that signatures don't have duplicate genes and spaces
  for(i in 1:nrow(signatures_df)){
    value <- gsub(" ", "", signatures_df[[i, "Gene_list"]])
    
    if(stringr::str_detect(value, ",")){
      value <- sort(unique(unlist(stringr::str_split(value, ","))))
    }
    signatures_df[[i, "Gene_list"]] <- paste0(value[value != ""], collapse=",")
  }
  
  # Check gene symbols within signatures against an Ensembl annotation
  signatures_df <- check_signature_genes(signatures_df)
  
  # Standardize signature names and gene lists
  signatures_df <- standardize_signatures_my(signatures_df) %>%
    dplyr::mutate(Gene_list = gsub(",,", ",", Gene_list)) %>%
    dplyr::mutate(Gene_list = gsub("^,|,$", "", Gene_list)) %>%
    dplyr::arrange(Signature_name)
  
  return(signatures_df)
}


standardize_signatures_my <- function(signatures_df){
  ## Standardize a set of signatures
  ## This function is NOT specific to an expression dataset
  
  # Remove duplicate signatures
  # Remove all punctuation except underscore from signature names
  signatures_df <- signatures_df %>%
    dplyr::group_by(Gene_list) %>%
    dplyr::summarise(Signature_name = paste0(unique(Signature_name), collapse="_")) %>%
    dplyr::mutate(Signature_name = gsub("_+", "_", gsub(" ", "_", gsub("[[:punct:]]", "_", Signature_name)))) %>%
    dplyr::ungroup()
  
  # Remove duplicate words from signature names
  signatures_df[["Signature_name"]] <- sapply(signatures_df[["Signature_name"]], function(x) 
    paste(unique(unlist(stringr::str_split(x, "_"))), collapse="_"))
  
  # Add an ID column
  signatures_df[["ID"]] <- 1:nrow(signatures_df)
  
  # Find duplicate signature names
  # Fix any signature names starting with a number
  signatures_df <- signatures_df %>%
    dplyr::group_by(Signature_name) %>%
    dplyr::mutate(Count = dplyr::n_distinct(Gene_list)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Signature_name = ifelse(Count == 1, Signature_name, paste0(Signature_name, "_", ID))) %>%
    dplyr::mutate(Signature_name = gsub("^\\d+_?", "", Signature_name))
  
  # Make sure that signature names are no more than X words
  signatures_df[["Signature_name"]] <- sapply(signatures_df[["Signature_name"]], function(x) 
    gsub("_NA", "", paste(unlist(stringr::str_split(x, "_"))[1:5], collapse="_")))
  
  return(signatures_df %>%
           dplyr::select(Signature_name, Gene_list))
}


check_signature_genes <- function(signatures_df, species="Rn"){
  ## Check gene symbols within signatures against an Ensembl annotation
  ## Species: Hs, Mm, Rn
  
  # sig_global_path <- "XXXX"
  sig_global_path <- "C:/Users/lyubetsa/Documents/Data/Signatures/"
  
  # Load gene symbol list
  genes_df <- readr::read_delim(file=paste0(sig_global_path, "hg38_symbol2ensembl.txt"), delim="\t")
  
  # Create a list of lists of signatures
  signature_list <- as.list(sapply(signatures_df$Gene_list, function(x) strsplit(x, ",")))
  names(signature_list) <- signatures_df$Signature_name

  # For each gene find a list of signatures where its mentioned
  # Check gene symbols against an Ensembl reference
  gene_sig_df <- invert_list_my(signature_list) %>%
    dplyr::filter(Symbol %in% unique(genes_df$Symbol) == FALSE)
  
  # Attempt to find a true gene symbol using limma
  gene_names_new <- sapply(gene_sig_df$Symbol, function(x) limma::alias2Symbol(x, species=species))
  gene_sig_df[["Symbol_new"]] <- unname(gene_names_new)
  
  # Select those genes that were uniquely identified and match that to Ensembl reference
  gene_sig_new_df <- gene_sig_df %>%
    dplyr::filter(Symbol_new %in% unique(genes_df$Symbol) == TRUE)
  
  # Substitute old gene symbols for new ones where possible
  gene_list <- as.list(gene_sig_new_df$Symbol_new)
  names(gene_list) <- gene_sig_new_df$Symbol
  
  # Substitute faulty gene symbols for good ones
  signatures_new_df <- signatures_df
  for(g in names(gene_list)){
    for(rnum in 1:nrow(signatures_new_df)){
      signatures_new_df[rnum, 3] <- gsub(paste0("\\b", g, "\\b"), gene_list[[g]], signatures_new_df[[rnum, 3]])
    }
  }

  # Remove unresolved genes
  #for(g in gene_sig_df$Symbol){
  #  for(rnum in 1:nrow(signatures_new_df)){
  #    signatures_new_df[rnum, 3] <- gsub(paste0("\\b", g, "\\b"), "", signatures_new_df[[rnum, 3]])
  #  }
  #}
  
  return(signatures_new_df)
}


parse_signatures_rds_my <- function(){
  ## Read Namit's RDS with signatures and output a two column table
  ## Output: signature name, comma-delimited gene list
  
  # sig_global_path <- "C:/Users/lyubetsa/Documents/Data/Signatures/"
  sig_global_path <- "XXXX"
  
  
  ## Read GeneTracker output
  
  # sig_path <- paste0(sig_global_path, "Originals/listtracker_20200105.rds")
  # gene_list <- readRDS(sig_path)
  # outfile <- paste0(sig_global_path, "/signatures_listtracker_20200105.txt")
  
  ## Read Silpa's RData (this is not an RDS!)
  
  # prefix <- "suthrams_myelod_"
  # sig_path <- paste0(sig_global_path, "Originals/suthrams_myeloidGeneSet_20200910.rds")
  # load(sig_path)
  # gene_list <- myeloidGeneSet
  # outfile <- paste0(sig_global_path, "signatures_myeloidGeneSet_suthrams_20200910.txt")
  
  ## Read Andy's GMT file
  ## BiocManager::install("qusage")
  
  prefix <- "brownea_pdac_"
  sig_path <- paste0(sig_global_path, "Originals/brownea_pdac_allSigs_072320.gmt")
  gene_list <- qusage::read.gmt(sig_path)
  outfile <- paste0(sig_global_path, "signatures_pancreas_brownea_20200723.txt")
  
  
  ## Parse and write signatures to a standard file
  
  sig_df <- tibble::tibble("Signature_name" = paste0(prefix, names(gene_list)), 
                           "Gene_list" = sapply(gene_list, function(x) paste0(sort(unique(x)), collapse=","))) %>%
    dplyr::arrange(Signature_name)

  readr::write_delim(sig_df, file=outfile, delim = "\t", append=FALSE, col_names = TRUE)
}


parse_signatures_table_my <- function(){
  ## Read a table where every row corresponds to one annotated gene
  ## Output: signature name, comma-delimited gene list

  # sig_global_path <- "XXXX"  
  input_table <- "C:/Users/lyubetsa/Documents/Projects_NewTech/References/RatColon/1-s2.0-S0092867420309946-mmc4.xlsx"
  output_file <- "C:/Users/lyubetsa/Documents/Data/Signatures/processed_signatures_rat_colon.txt"

  # Signature prefix
  prefix <- "RC_"
  
  # Annotation columns
  compartment_col <- "Broad category"
  cell_type_col <- "Cell type"
  cell_subtype_col <- "ident"
  # Gene symbol column
  gene_col <- "gene"
  
  # Load the table
  genes_df <- readxl::read_xlsx(input_table, sheet = "Mouse Colon All") 
  
  # Wrangle the table
  signatures_df <- genes_df %>%
    dplyr::mutate(Compartment = tolower(gsub("[[:punct:]]+", "-", !!rlang::sym(compartment_col))),
                  CellType = tolower(gsub("[[:punct:]]+", "-", !!rlang::sym(cell_type_col))),
                  CellSubtype = tolower(gsub("[[:punct:]]+", "-", !!rlang::sym(cell_subtype_col)))) %>%
    # dplyr::select(dplyr::all_of(c("Compartment", "CellType", "CellSubtype", gene_col))) %>%
    dplyr::mutate(Signature_name = paste(Compartment, CellType, sep="_")) %>%
    dplyr::group_by(Signature_name) %>%
    dplyr::summarise(Gene_list = paste0(sort(unique(!!rlang::sym(gene_col))), collapse=",")) %>%
    dplyr::select(Signature_name, Gene_list)

  # Check gene symbols within signatures against an Ensem,bl annotation
  signatures_df <- check_signature_genes(signatures_df, species="Rn")
  
  # Standardize signature names and gene lists
  signatures_df <- standardize_signatures_my(signatures_df)

  readr::write_delim(signatures_df, output_file, delim = "\t", append=FALSE, col_names = TRUE)
}
