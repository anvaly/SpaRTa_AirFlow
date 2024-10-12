# Author: Anna Lyubetskaya. Date: 19-11-21

# Packages uses:
# fgsea: https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
# Example with fgsea: https://www.pkimes.com/benchmark-fdr-html/additionalfile25_fgsea-human.html
# KEGG.db: http://bioconductor.org/packages/release/data/annotation/html/KEGG.db.html
# General resources: http://bioconductor.org/packages/release/BiocViews.html#___Pathways

# Packages to consider:
# msigdbr:https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html


# Important biomaRt helper functions:
# pages <- attributePages(ensembl)
# page_attributes <- listAttributes(ensembl, page="feature_page")

# Tool to try: Camera (pathway enrichment analysis): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3458527/


source("utils_biomart.R")


run_fgsea_my <- function(gsea_rank_list, annotation_df, filename){
  ## Perform GSEA via fgsea()
  ## todo: generalize id_go, name_go, and Hugo_Symbol to handle both GO and KEGG
  
  # Create a named list with names=term_id and values=gene_vector
  term_df <- annotation_df %>%
    dplyr::select(id_go, name_go, Hugo_Symbol) %>%
    dplyr::group_by(id_go, name_go) %>%
    dplyr::summarise(gene_list = paste(sort(unique(Hugo_Symbol)), collapse=";"))
  
  term_list <- setNames(as.list(strsplit(term_df$gene_list, ";")), term_df$id_go)
  
  # Perform GSEA test
  term_gsea <- fgsea::fgsea(term_list, gsea_rank_list, nperm=10000, maxSize=500, minSize=5) %>%
    tibble::as_tibble() %>%
    dplyr::inner_join(term_df, by=c("pathway" = "id_go")) %>%
    dplyr::select(-c(ES, NES, nMoreExtreme, leadingEdge)) %>%
    dplyr::arrange(pval)
  
  
}


annotate_symbol_go_kegg_my <- function(gene_ids){
  ## Extract GO and KEGG information for provided gene lists
  
  gene_ids <- gene_ids %>%
    dplyr::inner_join(ensembl2entrez_my(gene_ids$ensembl_gene_id), by="ensembl_gene_id")
  
  gene_go <- gene_ids %>%
    dplyr::inner_join(ensembl2go_my(gene_ids$ensembl_gene_id), by="ensembl_gene_id") %>%
    dplyr::filter(id_go != "")
  
  gene_kegg <- gene_ids %>%
    dplyr::inner_join(entrez2kegg_my(gene_ids$entrezgene_id), by="entrezgene_id") %>%
    dplyr::filter(id_kegg != "")
  
  return(list("go" = gene_go,
              "kegg" = gene_kegg))
}


entrez2kegg_my <- function(gene_list){
  ## Translate Entrez IDs to KEGG pathway IDs and names
  ## KEGG.db: KEGGPATHID2NAME, KEGGEXTID2PATHID
  ## Other resource: org.Hs.eg.db::org.Hs.egPATH
  ## Keep all values as characters for compatibility
  
  library(KEGG.db)
  
  gene_list <- unique(gene_list[!is.na(gene_list) & !is.null(gene_list)])
  kegg_ids <- sub("hsa", "", unlist(mget(gene_list, KEGGEXTID2PATHID, ifnotfound=NA)))
  
  kegg_list <- as.character(unique(unname(kegg_ids[!is.na(kegg_ids) & !is.null(kegg_ids)])))
  kegg_descriptions <- unlist(mget(kegg_list, KEGGPATHID2NAME, ifnotfound=NA))
  
  return(tibble::enframe(kegg_ids, name="entrezgene_id", value="id_kegg") %>%
           dplyr::inner_join(tibble::enframe(kegg_descriptions, name="id_kegg", value="name_kegg"), by="id_kegg"))
}
