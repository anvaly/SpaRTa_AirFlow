# Author: Anna Lyubetskaya. Date: 20-01-27
# Parse various syngeneic signatures mined by Brian Rabe
# /XXXX


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_tibble.R")
source("code/utils/utils_gene_lists.R")


## PARAMETERS ----


date <- "Aug21"
species <- "Mm"  ## rat, human, mouse

# Thresholds
alpha_thr <- 0.1
pvalue_thr <- 0.01
mu_thr <- 0.5
fc_thr <- 0.5
topn <- 100

# Files to read: File name, sheet name, column names, and the number at which data starts
file_list <- list(
  # c("", "", "Tcell", "NK", "Bcell", "M", "CAF", "B16F10", "EMT6", "LL2", "CT26", "MC38", "RBC", "pDC"))
  list("Kumar_2018_table_s3.xlsx", "1-s2.0-S221112471831636X-mmc4 (", 2),
  # c("Symbol", "Gene symbol", "melanocyt", "pigment", "tyrosinas", "melanoma", "catenin", "Melanoma vs others") - pick UP for the last column
  list("Rambow_2015_table_s1.xlsx", "S1D", 2),
  # No first column, just a list of genes
  list("Sehgal_2020_table_s6.xlsx", "Table S6 E-M signature genes", 1),
  # c("id", "log2 Fold Change", "pct.1", "pct.2", "p value", "Adjusted p value")
  list("Sehgal_2020_table_s7.xlsx", "Table S7 E-M positive cell expr", 1),
  # list("Gene Symbol", "Cell Type", Cell Sub Type", "Compartment")
  list("syngeneic_signatures.xlsx", "syngeneic_biomarkers", 1))


## PATHS ----


# Input folder
input_path <- "XXXX"

# Output folder
filename <- paste0("signatures_syngeneics_t", topn, "_", date, ".txt")
output_file1 <- paste0("data/import/Signatures/", filename)
output_file2 <- paste0("XXXX", filename)


## "Kumar_2018_table_s3.xlsx" ----


data_k18_df <- readxl::read_excel(paste0(input_path, file_list[[1]][[1]]), sheet=file_list[[1]][[2]], skip=file_list[[1]][[3]]-1) %>%
  dplyr::select(-`...1`) %>%
  dplyr::rename(Symbol = `...2`) %>%
  df_wide2long_my(key="Group", val="Presence") %>%
  tidyr::drop_na() %>%
  edit_group_name_my("Group") %>%
  edit_symbol_my("Symbol") %>%
  dplyr::mutate(Source = "K18",
                Compartment = "",
                Subtype = "",
                Group = paste(Source, Compartment, Group, Subtype, sep=".")) %>%
  dplyr::select(Group, Symbol)

print(sort(table(data_k18_df$Group)))


## "Rambow_2015_table_s1.xlsx" ----


data_r15_df <- readxl::read_excel(paste0(input_path, file_list[[2]][[1]]), sheet=file_list[[2]][[2]], skip=file_list[[2]][[3]]-1) %>%
  dplyr::filter(`Melanoma vs others` == "UP") %>%
  dplyr::select(-Symbol, `Melanoma vs others`) %>%
  dplyr::rename(Symbol = `Gene symbol`) %>%
  df_wide2long_my(key="Subtype", val="Presence") %>%
  dplyr::filter(Presence == "yes") %>%
  edit_group_name_my("Subtype") %>%
  edit_symbol_my("Symbol") %>%
  dplyr::mutate(Source = "R15",
                Compartment = "",
                Group = "B16",
                Subtype = "",
                Group = paste(Source, Compartment, Group, Subtype, sep=".")) %>%
  dplyr::select(Group, Symbol)

print(sort(table(data_r15_df$Group)))


## "Sehgal_2020_table_s6.xlsx" ----


data_s20_6_df <- readxl::read_excel(paste0(input_path, file_list[[3]][[1]]), sheet=file_list[[3]][[2]], skip=file_list[[3]][[3]]-1)
data_s20_6_df <- data_s20_6_df %>%
  dplyr::add_row(Mmp2 = "Mmp2") %>%
  dplyr::rename(Symbol = Mmp2) %>%
  dplyr::mutate(Group = "MC38") %>%
  edit_group_name_my("Group") %>%
  edit_symbol_my("Symbol") %>%
  dplyr::mutate(Source = "S20_6",
                Compartment = "",
                Subtype = "postIO_intersect_EMTHallmark",
                Group = paste(Source, Compartment, Group, Subtype, sep="."))

print(sort(table(data_s20_6_df$Group)))


## "Sehgal_2020_table_s7.xlsx" ----


data_s20_7_df <- readxl::read_excel(paste0(input_path, file_list[[4]][[1]]), sheet=file_list[[4]][[2]], skip=file_list[[4]][[3]]-1) %>%
  dplyr::rename(Symbol = id) %>%
  dplyr::filter(`log2 Fold Change` >= fc_thr & `Adjusted p value` <= pvalue_thr & pct.1 >= alpha_thr) %>%
  dplyr::mutate(Group = "MC38") %>%
  edit_group_name_my("Group") %>%
  edit_symbol_my("Symbol") %>%
  dplyr::mutate(Source = "S20_7",
                Compartment = "",
                Subtype = "HighEMT_cells",
                Group = paste(Source, Compartment, Group, Subtype, sep=".")) %>%
  dplyr::select(Group, Symbol)

print(sort(table(data_s20_7_df$Group)))


## "syngeneic_signatures.xlsx" ----


data_manual_df <- readxl::read_excel(paste0(input_path, file_list[[5]][[1]]), sheet=file_list[[5]][[2]], skip=file_list[[5]][[3]]-1) %>%
  dplyr::rename(Symbol = `Gene Symbol`, CellType = `Cell Type`, CellSubType = `Cell Sub Type`) %>%
  edit_group_name_my("Compartment") %>%
  edit_group_name_my("CellType") %>%
  edit_group_name_my("CellSubType") %>%
  edit_symbol_my("Symbol") %>%
  dplyr::mutate(Group = paste(Compartment, CellType, CellSubType, sep="."),
                Group = gsub("\\/", "-", gsub("\\+", "", Group))) %>%
  edit_group_name_my("Group") %>%
  dplyr::mutate(Group = paste("U", Group, sep=".")) %>%
  dplyr::select(Group, Symbol)

print(sort(table(data_manual_df$Group)))


## WRITE SIGNATURES TO FILE ----


data_list <- list(data_manual_df, data_k18_df, data_r15_df, data_s20_6_df, data_s20_7_df)

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
    dplyr::summarise(Gene_list = paste(sort(unique(Symbol)), collapse=",")) %>%
    dplyr::mutate(Signature_name = gsub("\\bnone\b|\\bNone\\b", "", Signature_name),
                  Signature_name = gsub("\\.+", "\\.", Signature_name),
                  Signature_name = gsub("^\\.|\\.$", "", Signature_name))
  
  readr::write_delim(df, output_file1, append=TRUE, delim="\t")
  readr::write_delim(df, output_file2, append=TRUE, delim="\t")
}
