# Author: Anna Lyubetskaya. Date: 20-01-09

# Important biomaRt helper functions:
# mart <- biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="useast.ensembl.org")
# pages <- biomaRt::attributePages(mart)
# page_attributes <- biomaRt::listAttributes(mart, page="feature_page")

# Getting archived versions
# archive_list <- biomaRt::listEnsemblArchives()
# name  date  url version current_release
# Ensembl GRCh37	Feb 2014	http://grch37.ensembl.org	GRCh37
# Ensembl 99	Jan 2020	http://jan2020.archive.ensembl.org	99
# Ensembl 98	Sep 2019	http://sep2019.archive.ensembl.org	98
# Ensembl 97	Jul 2019	http://jul2019.archive.ensembl.org	97
# Ensembl 96	Apr 2019	http://apr2019.archive.ensembl.org	96


access_biomart_my <- function(attribute_list, filter_field, filter_list, dataset="hsapiens_gene_ensembl", version=NULL){
  ## General function for querying Ensembl DB via biomaRt
  ## dataset = "hsapiens_gene_ensembl" or "mmusculus_gene_ensembl"
  ## host = uswest.ensembl.org or useast.ensembl.org

  # Connect Ensembl data via BioMart
  mart <- biomaRt::useMart(biomart="ensembl", dataset=dataset, host="useast.ensembl.org")

  # Grab selected fields using a filter
  biomart_data <- biomaRt::getBM(attributes=attribute_list, mart = mart, 
                                 filters = filter_field, values = filter_list)
  
  return(biomart_data)
}


ensembl2entrez_my <- function(gene_list, dataset="hsapiens_gene_ensembl"){
  ## Translate Ensembl IDs to Entrez IDs via BioMart
  
  attribute_list <- c("ensembl_gene_id", "entrezgene_id")
  filter_field <- "ensembl_gene_id"
  filter_list <- unique(gene_list)
  
  # Download Entrez IDs for a list of Ensembl IDs  
  biomart_data <- access_biomart_my(attribute_list, filter_field, filter_list, dataset=dataset) %>%
    dplyr::mutate(entrezgene_id = as.character(entrezgene_id))
  
  return(biomart_data)
}


ensembl2go_my <- function(gene_list, dataset="hsapiens_gene_ensembl"){
  ## Translate Ensembl IDs to GO IDs (and GO names) via BioMart
  
  attribute_list <- c("ensembl_gene_id", "go_id", "name_1006")
  filter_field <- "ensembl_gene_id"
  filter_list <- unique(gene_list)
  
  # Download GO IDs and their names for a list of Ensembl IDs
  biomart_data <- access_biomart_my(attribute_list, filter_field, filter_list, dataset=dataset) %>%
    dplyr::rename("name_go" = "name_1006", "id_go" = "go_id")

  return(biomart_data)
}


refseq2ensembl_my <- function(gene_list, dataset="mmusculus_gene_ensembl"){
  ## Find gene lengths from Ensembl IDs via BioMart
  
  attribute_list <- c("refseq_mrna", "ensembl_gene_id")
  filter_field <- "refseq_mrna"
  filter_list <- unique(gene_list)
  
  # Download Ensembl IDs for a list of RefSeq mRNA IDs
  biomart_data <- access_biomart_my(attribute_list, filter_field, filter_list, dataset=dataset) %>%

  return(biomart_data)
}


ensembl2length_my <- function(gene_list, dataset="mmusculus_gene_ensembl"){
  ## Find gene lengths from Ensembl IDs via BioMart
  
  attribute_list <- c("ensembl_gene_id", "transcript_length")
  filter_field <- "ensembl_gene_id"
  filter_list <- unique(gene_list)
  
  # Download gene lengts for a list of Ensembl IDs
  biomart_data <- access_biomart_my(attribute_list, filter_field, filter_list, dataset=dataset) %>%

  return(biomart_data)
}


ensembl2seq_my <- function(gene_list, dataset="mmusculus_gene_ensembl", version=NULL){
  ## Find protein sequences Ensembl IDs via BioMart
  ## Full coding sequence field: "coding"
  ## Exon coordinates are always provided on the + strand. Start < End
  ## Exon sequence is always provided according to the strand: forward or reverse
  
  # Exon sequences
  attribute_list <- c("external_gene_name", "ensembl_transcript_id", 
                      "transcript_length", "cds_length", "coding")
  filter_field <- "ensembl_transcript_id"
  filter_list <- unique(gene_list)
  
  # Download gene lengts for a list of Ensembl IDs
  biomart_data <- access_biomart_my(attribute_list, filter_field, filter_list, dataset=dataset, version=version)

  return(biomart_data)
}


symbol2ensembl_my <- function(gene_list, dataset="mmusculus_gene_ensembl"){
  ## Find ensembl gene ID for every gene symbol
  
  attribute_list <- c("external_gene_name", "ensembl_gene_id")
  filter_field <- "external_gene_name"
  filter_list <- unique(gene_list)
  
  biomart_data <- access_biomart_my(attribute_list, filter_field, filter_list, dataset=dataset)
  
  return(biomart_data)
}
  
  
symbol2coordinates_my <- function(gene_list, dataset="mmusculus_gene_ensembl"){
  ## Find CDS starts and ends, as well as transcript IDs for every gene symbol

  # Gather all transcripts for each exon
  attribute_list <- c("external_gene_name", "ensembl_transcript_id")  #"ensembl_exon_id"
  filter_field <- "external_gene_name"
  filter_list <- unique(gene_list)
  
  biomart_data_annot1 <- access_biomart_my(attribute_list, filter_field, filter_list, dataset=dataset)
  
  attribute_list <- c("ensembl_transcript_id", "ensembl_exon_id")#"genomic_coding_start", "genomic_coding_end")
  filter_field <- "ensembl_transcript_id"
  filter_list <- unique(biomart_data_annot1$ensembl_transcript_id)
  
  biomart_data_annot2 <- access_biomart_my(attribute_list, filter_field, filter_list, dataset=dataset)

  attribute_list <- c("ensembl_exon_id", "chromosome_name", "strand", "genomic_coding_start", "genomic_coding_end")
  filter_field <- "ensembl_exon_id"
  filter_list <- unique(biomart_data_annot2$ensembl_exon_id)
  
  biomart_data_annot3 <- access_biomart_my(attribute_list, filter_field, filter_list, dataset=dataset)

  return(biomart_data_annot1 %>%
           dplyr::inner_join(biomart_data_annot2, by="ensembl_transcript_id") %>%
           dplyr::inner_join(biomart_data_annot3, by="ensembl_exon_id") %>% 
           tidyr::drop_na() %>%
           unique) 
}


gene2coordinates_my <- function(gene_list, dataset="hsapiens_gene_ensembl"){
  ## Find gene chromosome, strand, start, and stop
  
  attribute_list <- c("ensembl_gene_id", "external_gene_name", "chromosome_name", "strand", 
                      "genomic_coding_start", "genomic_coding_end")
  filter_field <- "ensembl_gene_id"
  filter_list <- unique(gene_list)
  
  biomart_data_annot <- access_biomart_my(attribute_list, filter_field, filter_list, dataset=dataset)
  
  return(biomart_data_annot) 
}
