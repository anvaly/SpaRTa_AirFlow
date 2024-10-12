# Author: Anna Lyubetskaya. Date: 21-02-28


markers_cluster_profile_my <- function(markers_df, gene_ref_list, reference_db, ref_short, 
                                       gsea_qvalue_threshold=0.1, gsea_refs=c("KEGG", "MSigDB"), 
                                       output_path=NULL, direction_filt=c("UP")){
  ## Perform gene over-representation analysis using clusterProfiler
  ## http://yulab-smu.top/clusterProfiler-book/index.html
  ## Tai Wang prefers the following references: KEGG, MSigDB C2, WikiPathways, Panther, ENCODE (TF), GO, Reactome
  
  
  # Define gene lists
  gsea_gene_lists_df <- markers_df %>% 
    dplyr::filter(direction %in% direction_filt) %>%
    dplyr::group_by(cluster, direction) %>% 
    dplyr::summarise(Gene_list = paste(sort(unique(Symbol)), collapse=";"),
                     Symbol_count = dplyr::n_distinct(Symbol)) %>%
    dplyr::mutate(cluster = as.character(cluster)) %>%
    dplyr::filter(Symbol_count >= 3)
  
  
  # Collect enrichment data here
  enrichment_list <- list()
  
  # Go through various references
  for(ref_name in gsea_refs){
    
    
    # Go through groups / gene lists
    for(i in 1:nrow(gsea_gene_lists_df)){
      
      # Target gene list
      gene_sym_list <- unname(unlist(strsplit(gsea_gene_lists_df[[i, "Gene_list"]], ";")))
      
      if(ref_name == "KEGG"){
        
        # KEGG over-representation analysis
        res <- clusterProfiler_kegg_my(gene_sym_list, gene_ref_list, ref_short, reference_db)
        
      } else if(ref_name == "MSigDB"){
        
        # MSigDB over-representation analysis - C2 as default reference
        res <- clusterProfiler_msigdb_my(gene_sym_list, gene_ref_list)
      }
      
      if(!is.null(res)){
        
        # Filter results of the over-representation alaysis
        filt_res <- res %>%
          clusterProfiler_filter_my(gsea_qvalue_threshold) %>%
          dplyr::mutate(cluster = gsea_gene_lists_df[[i, "cluster"]],
                        direction = gsea_gene_lists_df[[i, "direction"]],
                        DB = ref_name) %>%
          dplyr::relocate(cluster, direction) %>%
          tibble::remove_rownames()
        
        enrichment_list[[paste0(ref_name, i)]] <- filt_res
        
      }
    }
    
  }
  
  enrichment_df <- dplyr::bind_rows(enrichment_list)
  
  # Write GSEA result to a file
  if(nrow(enrichment_df) > 0){
    filename <- paste0(output_path, "_", ref_name, ".txt")
    readr::write_delim(enrichment_df, filename, delim="\t")
  }
  
  
  return(enrichment_df)
}


clusterProfiler_msigdb_my <- function(gene_sym_list, gene_ref_list){
  # Use clusterProfiler to perform gene over-representation analysis against MSigDB
  # Todo: add mouse reference from http://bioinf.wehi.edu.au/MSigDB/
  
  # MSigDB gene dataset
  msigdb_h <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, human_gene_symbol)
  
  msigdb_c2 <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2") %>% 
    dplyr::select(gs_name, human_gene_symbol)
  
  msigdb_c5 <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5") %>% 
    dplyr::select(gs_name, human_gene_symbol)
  
  # Hallmark over-representation analysis - always human gene Symbols as inputs
  if(length(intersect(msigdb_h$human_gene_symbol, gene_sym_list)) > 0){
    res <- clusterProfiler::enricher(gene = gene_sym_list, universe = gene_ref_list, 
                                     TERM2GENE=dplyr::bind_rows(msigdb_h, msigdb_c5))
    
    # Edit the gene list
    if(nrow(res@result) > 0){
      res@result[["Gene_list"]] <- gsub("/", ";", res@result[["geneID"]]) 
    }
  } else{
    res <- NULL
  }
  
  return(res)
}


clusterProfiler_kegg_my <- function(gene_sym_list, gene_ref_list, ref_short, reference_db){
  # Use clusterProfiler to perform gene over-representation analysis against KEGG
  
  
  ## Wrange the target gene list
  
  # Correct case for non-human gene list
  if(ref_short != "human" && ref_short != "hsa"){
    gene_sym_list <- stringr::str_to_title(gene_sym_list)
    gene_ref_list <- stringr::str_to_title(gene_ref_list)
  }
  
  # Match symbols to Entrez IDs
  eg <- clusterProfiler::bitr(gene_sym_list, fromType="SYMBOL", toType="ENTREZID", OrgDb=reference_db)
  
  
  ## Wrange the reference gene list
  
  # Match symbols to Entrez IDs
  eg_ref <- clusterProfiler::bitr(gene_ref_list, fromType="SYMBOL", toType="ENTREZID", OrgDb=reference_db)
  
  # Create a dictionary of entrez to symbol gene names
  gene_ref_dict <- eg_ref$SYMBOL
  names(gene_ref_dict) <- eg_ref$ENTREZID
  
  
  ## Run analysis
  
  # KEGG over-representation analysis - EntrezID specific to the organism
  res <- clusterProfiler::enrichKEGG(gene = eg$ENTREZID, universe = eg_ref$ENTREZID, organism = ref_short, keyType = "kegg")
  
  # Add symbol gene list
  if(!is.null(res)){
    res@result[["Gene_list"]] <- sapply(res@result$geneID, 
                                        function(x) paste(unname(gene_ref_dict[unlist(strsplit(x, "/"))]), collapse=";"))
  }
  
  return(res)
}


clusterProfiler_filter_my <- function(res, gsea_qvalue_threshold){
  ## Simplify and filter clusterProfiler data object to a tibble
  
  # Filtered results of the over-representation alaysis
  filt_res <- res@result %>%
    dplyr::filter(qvalue < gsea_qvalue_threshold) %>%
    dplyr::mutate(nlog10qvalue = round(-log10(qvalue), 1)) %>%
    dplyr::select(ID, Description, Count, GeneRatio, BgRatio, nlog10qvalue, Gene_list)
  
  return(filt_res) 
}


clusterProfiler_invert_my <- function(filt_res){
  ## Invert the filtered clusterProfiler tibble to create a gene-centric tibble
  
  if(nrow(filt_res) > 0){  
    # List of gene lists with KEGG term descriptions as names
    res_gene_list <- sapply(filt_res[["Gene_list"]], function(x) list(unname(unlist(strsplit(x, ";")))))
    names(res_gene_list) <- filt_res[["Description"]]
    
    # List of genes appearing in the over-representation analysis with the corresponding list of terms
    invert_filt_res <- invert_list_my(res_gene_list)
  } else{
    invert_filt_res <- NULL
  }
  
  return(invert_filt_res)
}
