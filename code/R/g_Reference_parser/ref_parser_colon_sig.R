# Author: Anna Lyubetskaya. Date: 20-01-04
# Parse various colon signatures mined by Eugene Drokhlyansky
# /XXXX


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_tibble.R")
source("code/utils/utils_gene_lists.R")


## PARAMETERS ----


date <- "Mar21_short"
species <- "Rn"  ## rat, human, mouse

# Thresholds
alpha_thr <- 0.1
pvalue_thr <- 0.01
mu_thr <- 0.5
fc_thr <- 1
topn <- 100

# Files to read: File name, sheet name, column names, and the number at which data starts
file_list <- list(
  # c("ident", "gene", "alpha", "padjD", "padjC", "mastfc", "log2fc"))
  list("Drokhlyansky_2020_table_4.xlsx", "Mouse Colon All", 1),
  # Use as primary if performs well: list("Enteroendocrine", "Enterocyte", "Enterocyte progenitor (early)", "Enterocyte progenitor (late)", "Goblet", "Paneth", "Stem", "TA", "Tuft")
  list("Haber_2018_Table_3.xlsx", "Summary", 6),
  # Use as primary if Table 3 fails: Main table: list("Goblet", "Paneth", "Tuft", "Enteroendocrine", "Enterocyte (Proximal)", "Enterocyte (Distal)")
  list("Haber_2018_Table_4.xlsx", "Summary", 6),
  # list("Activated CD4 T", "Activated CD4 T adj-pval", "Activated CD4 T log fold change", ...)
  list("James_2020_Supp.xlsx", "Table S4", 1),
  # list("ident", "gene", "alpha", "padjD", "padjC", "mastfc", "log2fc"),
  list("Smillie_2019_colon_s2.xlsx", c("Epithelial", "Stromal", "Immune"), 1),  # "Subclusters"
  # list("ident", "gene", "alpha", "padjD", "padjC", "mastfc", "log2fc")
  list("Smillie_2019_colon_s3.xlsx", c("TFs", "GPCRs", "Transporters", "PRRs", "Cytokines"), 1),
  # list("Gene Symbol", "Cell Type", Cell Sub Type", "Compartment")
  list("signatures_colon.xlsx", "colon_biomarkers", 1))

file_supp_list <- list(
  # "Name", "Cell subsets"
  list("Smillie_2019_colon_s2.xlsx", "Lineage", 1)
)


## PATHS ----


# Input folder
input_path <- "XXXX"

# Output folder
filename <- paste0("signatures_colon_t", topn, "_", date, ".txt")
output_file1 <- paste0("data/import/Signatures/", filename)
output_file2 <- paste0("XXXX", filename)


## "Drokhlyansky_2020_table_4.xlsx" ----


data_d20_df <- readxl::read_excel(paste0(input_path, file_list[[1]][[1]]), sheet=file_list[[1]][[2]], skip=file_list[[1]][[3]]-1) %>%
  dplyr::filter(alpha >= alpha_thr & mastfc >= fc_thr & 
                  padjD <= pvalue_thr & padjC <= pvalue_thr & 
                  (mu >= mu_thr | ref_mu >= mu_thr)) %>%
  dplyr::rename(Group = ident, Symbol = gene) %>%
  edit_group_name_my("Group") %>%
  edit_symbol_my("Symbol") %>%
  dplyr::mutate(Source = "D20",
                Compartment = "C",
                Subtype = ifelse(grepl("\\d+", Group), gsub(".+_(\\d+)", "c\\1", Group), "na"),
                Group = gsub("_\\d+", "", Group),
                Group = paste(Source, Compartment, Group, Subtype, sep=".")) %>%
  dplyr::group_by(Group) %>%
  dplyr::top_n(topn, mastfc) %>%
  dplyr::ungroup() %>%
  dplyr::select(Group, Symbol)

print(sort(table(data_d20_df$Group)))



## "Haber_2018_Table_3.xlsx" ----


data_h18a_df <- readxl::read_excel(paste0(input_path, file_list[[2]][[1]]), sheet=file_list[[2]][[2]], skip=file_list[[2]][[3]]-1) %>%
  df_wide2long_my(key="Group", val="Symbol", start_col=1) %>%
  tidyr::drop_na() %>%
  edit_group_name_my("Group") %>%
  edit_symbol_my("Symbol") %>%
  dplyr::mutate(Source = "H18a",
                Compartment = "C",
                Subtype = ifelse(grepl("\\(", Group), gsub(".+\\((.+)\\)", "\\1", Group), "na"),
                Group = gsub("_\\(.+\\)", "", Group),
                Group = paste(Source, Compartment, Group, Subtype, sep="."))

data_h18a_df[["Order"]] <- 1:nrow(data_h18a_df)
data_h18a_df <- data_h18a_df %>%
  dplyr::group_by(Group) %>%
  dplyr::top_n(topn, Order) %>%
  dplyr::ungroup() %>%
  dplyr::select(Group, Symbol)

print(sort(table(data_h18a_df$Group)))


## "Haber_2018_Table_4.xlsx" ----


data_h18b_df <- readxl::read_excel(paste0(input_path, file_list[[3]][[1]]), sheet=file_list[[3]][[2]], skip=file_list[[3]][[3]]-1) %>%
  df_wide2long_my(key="Group", val="Symbol", start_col=1) %>%
  tidyr::drop_na() %>%
  edit_group_name_my("Group") %>%
  edit_symbol_my("Symbol") %>%
  dplyr::mutate(Source = "H18b",
                Compartment = "C",
                Subtype = ifelse(grepl("\\(", Group), gsub(".+\\((.+)\\)", "\\1", Group), "na"),
                Group = gsub("_\\(.+\\)", "", Group),
                Group = paste(Source, Compartment, Group, Subtype, sep="."))


data_h18b_df[["Order"]] <- 1:nrow(data_h18b_df)
data_h18b_df <- data_h18b_df %>%
  dplyr::group_by(Group) %>%
  dplyr::top_n(topn, Order) %>%
  dplyr::ungroup() %>%
  dplyr::select(Group, Symbol)

print(sort(table(data_h18b_df$Group)))


## "James_2020_Supp.xlsx" ----


# Read initial data
data_j20_init_df <- readxl::read_excel(paste0(input_path, file_list[[4]][[1]]), sheet=file_list[[4]][[2]], skip=file_list[[4]][[3]]-1)
# List of groups
group_list <- unique(gsub(" adj-pval| log fold change", "", colnames(data_j20_init_df)))

# A dictionary to rename the groups
group_dict <- list("Activated CD4 T" = "immune.t_cell.cd4_activated", 
                   "B cell IgA Plasma" = "immune.b_cell.iga_plasma", 
                   "B cell IgG Plasma" = "immune.b_cell.igg_plasma", 
                   "Follicular B cell" = "immune.b_cell.follicular", 
                   "B cell memory" = "immune.b_cell.memory", 
                   "CD8 T" = "immune.t_cell.cd8",
                   "ILC" = "immune.ilc.na", 
                   "Lymphoid DC" = "immune.dc.lymphoid", 
                   "Monocyte" = "immune.monocyte.na", 
                   "Mast" = "immune.mast.na", 
                   "Macrophage" = "immune.macrophage.na", 
                   "LYVE1 Macrophage" = "immune.macrophage.lyve1",
                   "NK" = "immune.nk_cell.na", 
                   "Tcm" = "immune.t_cell.tcm", 
                   "Tfh" = "immune.t_cell.tfh", 
                   "Th1" = "immune.t_cell.th1", 
                   "Th17" = "immune.t_cell.th17", 
                   "Treg" = "immune.t_cell.treg",
                   "cDC1" = "immune.dc.cdc1", 
                   "cDC2" = "immune.dc.cdc2", 
                   "cycling DCs" = "immune.dc.cycling", 
                   "pDC" = "immune.dc.pdc", 
                   "gd T" = "immune.t_cell.gd", 
                   "cycling gd T" = "immune.t_cell.gd_cycling")

# Break up the data, filter, and put back together
data_j20_list <- list()
for(g in group_list){

  # Select columns corresponding to the group
  data_loc_df <- data_j20_init_df[, grep(paste0("^", g), colnames(data_j20_init_df))][, 1:3]
  # Rename oclumns
  colnames(data_loc_df) <- c(g, "adj_pval", "log_fold_change")
  
  # Filter data using thresholds
  data_loc_df <- data_loc_df %>%
    dplyr::filter(adj_pval <= pvalue_thr & log_fold_change >= fc_thr) %>%
    dplyr::top_n(topn, log_fold_change)
  
  # Add data to the list
  if(nrow(data_loc_df) >= 1){
    data_j20_list[[g]] <- tibble::tibble(Group = group_dict[[g]],
                                      Symbol = unlist(unname(c(data_loc_df[, 1]))))
  }
}

# Stack lists into a tibble
data_j20_df <- dplyr::bind_rows(data_j20_list) %>%
  edit_group_name_my("Group") %>%
  edit_symbol_my("Symbol") %>%
  dplyr::mutate(Group = paste("J20", Group, sep=".")) %>%
  dplyr::select(Group, Symbol)

print(sort(table(data_j20_df$Group)))


## "Smillie_2019_colon_s2.xlsx" - lineage data ----


# References for Smillie cell types and subtypes
#data_s19_ref_df <- readxl::read_excel(paste0(input_path, file_supp_list[[1]][[1]]), sheet=file_supp_list[[1]][[2]], skip=file_supp_list[[1]][[3]]-1)

# Create a named list of cell types and subtypes
#data_s19_ref_list <- sapply(data_s19_ref_df[["Cell subsets"]], function(x) stringr::str_split(toupper(x), ", "))
#names(data_s19_ref_list) <- data_s19_ref_df$Name
# Invert the list of lists
#data_s19_ref_rev_df <- invert_list_my(data_s19_ref_list)


## "Smillie_2019_colon_s2.xlsx" ----


# Cycle through sheets
data_s19a_list <- list()
for(compartment in file_list[[5]][[2]]){
  # Parse each sheet
  data_s19a_list[[compartment]] <- readxl::read_excel(paste0(input_path, file_list[[5]][[1]]), sheet=compartment, skip=file_list[[5]][[3]]-1) %>%
    dplyr::filter(alpha >= alpha_thr & mastfc >= fc_thr & 
                    padjD <= pvalue_thr & padjC <= pvalue_thr & 
                    (mu >= mu_thr | ref_mu >= mu_thr)) %>%
    dplyr::rename(Group = ident, Symbol = gene) %>%
    dplyr::group_by(Group) %>%
    dplyr::top_n(topn, mastfc) %>%
    dplyr::ungroup()
  
  # Add compartment name
  data_s19a_list[[compartment]]$Compartment <- compartment
}

data_s19a_df <- dplyr::bind_rows(data_s19a_list) %>%
  dplyr::mutate(Group = gsub("\\+", "", Group)) %>%
  edit_group_name_my("Group") %>%
  edit_group_name_my("Compartment") %>%
  edit_symbol_my("Symbol") %>%
  dplyr::mutate(Group = paste("S19a", Compartment, Group, sep=".")) %>%
  dplyr::select(Group, Symbol)

print(sort(table(data_s19a_df$Group)))


## "Smillie_2019_colon_s3.xlsx" ----


# Cycle through sheets
data_s19b_list <- list()
for(compartment in file_list[[6]][[2]]){
  # Parse each sheet
  data_s19b_list[[compartment]] <- readxl::read_excel(paste0(input_path, file_list[[6]][[1]]), sheet=compartment, skip=file_list[[6]][[3]]-1) %>%
    dplyr::filter(alpha >= alpha_thr & mastfc >= fc_thr & 
                    padjD <= pvalue_thr & padjC <= pvalue_thr & 
                    (mu >= mu_thr | ref_mu >= mu_thr)) %>%
    dplyr::rename(Group = ident, Symbol = gene) %>%
    dplyr::group_by(Group) %>%
    dplyr::top_n(topn, mastfc) %>%
    dplyr::ungroup()
  
  # Add compartment name
  data_s19b_list[[compartment]]$Compartment <- compartment
}

data_s19b_df <- dplyr::bind_rows(data_s19b_list) %>%
  dplyr::mutate(Group = gsub("\\+", "", Group)) %>%
  edit_group_name_my("Group") %>%
  edit_group_name_my("Compartment") %>%
  edit_symbol_my("Symbol") %>%
  dplyr::mutate(Group = paste("S19b", Compartment, Group, sep=".")) %>%
  dplyr::select(Group, Symbol)

print(sort(table(data_s19b_df$Group)))


## "signatures_colon.xlsx" ----


data_manual_df <- readxl::read_excel(paste0(input_path, file_list[[7]][[1]]), sheet=file_list[[7]][[2]], skip=file_list[[7]][[3]]-1) %>%
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


data_list <- list(data_d20_df, data_j20_df, data_h18a_df)  # data_s19a_df, data_s19b_df, data_h18b_df, data_manual_df

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
    dplyr::filter(!grepl("^MT-|^RP[SL]", Symbol)) %>%
    dplyr::rename(Signature_name = Group) %>%
    dplyr::group_by(Signature_name) %>%
    dplyr::summarise(Gene_list = paste(sort(unique(Symbol)), collapse=","))
  
  readr::write_delim(df, output_file1, append=TRUE, delim="\t")
  readr::write_delim(df, output_file2, append=TRUE, delim="\t")
}
