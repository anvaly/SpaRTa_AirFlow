# Author: Anna Lyubetskaya. Date: 21-05-19
# Parse various PDAC signatures mined by Brian Rabe, Eugene Drokhlyansky, and Anna Lyubetskaya
# /XXXX


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_tibble.R")
source("code/utils/utils_gene_lists.R")


## PARAMETERS ----


date <- "May21"
species <- "Hs"  ## rat, human, mouse

# Thresholds
alpha_thr <- 0.1
base_thr <- 1
pvalue_thr <- 0.01
mu_thr <- 0.5
fc_thr <- 1
topn <- 1000

# Files to read: File name, sheet name, column names, and the number at which data starts
file_list <- list(
  # Gene    p value    avg_logFC    pct.1    pct.2    Cell type
  list("pdac_peng2019_MOESM9.xls", "Sheet1", 2),  #1
  # EmptyColName    baseMean    baseMeanA    baseMeanB    foldChange    log2FoldChange    pval    padj
  list("panc_Muraro2016_CellSys_mmc3.xlsx", c("alpha", "beta", "delta", "pp", "epsilon", "duct", "acinar", "mesenchyme", "endothelial"), 1),  #4
  # Gene Symbol, Cell Type, Cell Sub Type, Compartment
  list("signatures_PDAC_combined.xlsx", "PDAC_biomarkers", 1)  #5
  # Gene    Weight
  # list("nicolle_table_s3", "Sheet 1", 3),  #2
  # Symbol    correlation    p_value    Entrez_Gene_ID    Definition    ILMN_Gene    Probe_Id
  #list("pek_table_s2.xlsx", "Sheet1 (2)", 2),  #3
  )


## PATHS ----


# Input folder
input_path <- "XXXX"

# Output folder
filename <- paste0("signatures_pdac_t", topn, "_", date, ".txt")
output_file1 <- paste0("data/import/Signatures/", filename)
output_file2 <- paste0("XXXX", filename)


## "pdac_peng2019_MOESM9.xlsx" ----


# Gene    p value    avg_logFC    pct.1    pct.2    Cell type
data_p19_df <- readxl::read_excel(paste0(input_path, file_list[[1]][[1]]), sheet=file_list[[1]][[2]], skip=file_list[[1]][[3]]-1) %>%
  dplyr::rename(Symbol = Gene, Group = `Cell type`) %>%
  dplyr::filter(avg_logFC >= fc_thr & `p value` <= pvalue_thr & pct.1 >= alpha_thr) %>%
  edit_group_name_my("Group") %>%
  edit_symbol_my("Symbol") %>%
  dplyr::mutate(Source = "P19",
                Group = paste(Source, Group, sep=".")) %>%
  dplyr::select(Group, Symbol)

print(sort(table(data_p19_df$Group)))


## "panc_Muraro2016_CellSys_mmc3.xlsx" ----


# EmptyColName    baseMean    baseMeanA    baseMeanB    foldChange    log2FoldChange    pval    padj
# Cycle through sheets
data_m16_list <- list()
for(compartment in file_list[[2]][[2]]){
  # Parse each sheet
  data_m16_list[[compartment]] <- readxl::read_excel(paste0(input_path, file_list[[2]][[1]]), sheet=compartment, skip=file_list[[2]][[3]]-1) %>%
    dplyr::filter(baseMean >= base_thr & log2FoldChange <= fc_thr*-1 & padj <= pvalue_thr) %>%
    dplyr::rename(Symbol = `...1`) %>%
    dplyr::mutate(Group = compartment) %>%
    dplyr::group_by(Group) %>%
    dplyr::top_n(topn, padj) %>%
    dplyr::ungroup()
}

data_m16_df <- dplyr::bind_rows(data_m16_list) %>%
  edit_group_name_my("Group") %>%
  edit_symbol_my("Symbol") %>%
  dplyr::mutate(Group = paste("M16", Group, sep=".")) %>%
  dplyr::select(Group, Symbol)

print(sort(table(data_m16_df$Group)))


## "signatures_PDAC_combined.xlsx" ----


data_manual_df <- readxl::read_excel(paste0(input_path, file_list[[3]][[1]]), sheet=file_list[[3]][[2]], skip=file_list[[3]][[3]]-1) %>%
  dplyr::rename(Symbol = `Gene Symbol`, CellType = `Cell Type`, CellSubType = `Cell Sub Type`) %>%
  edit_group_name_my("Compartment") %>%
  edit_group_name_my("CellType") %>%
  edit_group_name_my("CellSubType") %>%
  edit_symbol_my("Symbol") %>%
  dplyr::mutate(Group = paste(Compartment, CellType, CellSubType, sep="."),
                Group = gsub("\\.NA\\b", "", Group)) %>%
  edit_group_name_my("Group") %>%
  dplyr::mutate(Group = paste("U", Group, sep=".")) %>%
  dplyr::select(Group, Symbol)

print(sort(table(data_manual_df$Group)))


## WRITE SIGNATURES TO FILE ----


data_list <- list(data_manual_df, data_m16_df, data_p19_df)

write("Signature_name\tGene_list", output_file1)
write("Signature_name\tGene_list", output_file2)

for(d in data_list){
  # Compare the gene list against a reference extracted from 10X
  genes_list <- check_reference_my(d$Symbol, species=species)
  
  # Compare the gene list against a reference extracted from Ensembl via limma  
  if(length(genes_list$genes_excluded) > 0){
    gene_check <- check_limma_alias_my(genes_list$genes_excluded, species=species)
  }
  
  # Substitute genes with a reference synonym
  for(g in names(gene_check)){
    d$Symbol <- gsub(toupper(g), toupper(gene_check[[g]]), d$Symbol)
  }
  # Compare the gene list against a reference extracted from 10X
  genes_list <- check_reference_my(d$Symbol, species=species)
  
  df <- d %>%
    dplyr::filter(Symbol %in% genes_list$genes_included) %>%
    dplyr::rename(Signature_name = Group) %>%
    dplyr::group_by(Signature_name) %>%
    dplyr::summarise(Gene_list = paste(sort(unique(Symbol)), collapse=","))
  
  readr::write_delim(df, output_file1, append=TRUE, delim="\t")
  readr::write_delim(df, output_file2, append=TRUE, delim="\t")
}
