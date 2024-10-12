# Author: Anna Lyubetskaya. Date: 23-07-20
# Parse various PDAC signatures from Chijimatsu


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_tibble.R")
source("code/utils/utils_gene_lists.R")


## PARAMETERS ----


date <- "20230720"
species <- "Hs"  ## rat, human, mouse

# Thresholds
alpha_thr <- 0.2
pvalue_thr <- 0
fc_thr <- 1
topn <- 1000

# Files to read: File name, sheet name, column names, and the number at which data starts
file_list <- list(
  # cluster	gene	avg_log2FC	pct.1	pct.2	p_val	p_val_adj
  list("1-s2.0-S2589004222009312-mmc2.xlsx", c("Acinar cell", "B cell", "Ductal cell type 1", "Ductal cell type 2", 
                                              "Endocrine cell", "Endothelial cell", "Fibroblast", "Macrophage", "Stellate cell", "T cell"), 1),
  # classical signature (original)	in scRNAseq data	Average Expression in malignant ductal cells	Average Expression in normal ductal cells	Correlation to classical score (r value)	Correlation to basal-like score (r value)	New classical signature
  # basal-like signature (original)	in scRNAseq data	Average Expression in malignant ductal cells	Average Expression in normal ductal cells	Correlation to classical score (r value)	Correlation to basal-like score (r value)	New classical signature
  list("1-s2.0-S2589004222009312-mmc4.xlsx", c("classical", "basal"), 1)
)


## PATHS ----


# Input folder
input_path <- "XXXX"

# Output folder
filename <- paste0("signatures_pdac_", date, ".txt")
output_file1 <- paste0("data/import/Signatures/", filename)
output_file2 <- paste0("XXXX", filename)


## Chijimatsu cell types ----
# https://www.sciencedirect.com/science/article/pii/S2589004222009312?via%3Dihub#appsec2


# cluster	gene	avg_log2FC	pct.1	pct.2	p_val	p_val_adj

# Cycle through sheets
data_chiji_list <- list()
for(compartment in file_list[[1]][[2]]){
  # Parse each sheet
  data_chiji_list[[compartment]] <- readxl::read_excel(paste0(input_path, file_list[[1]][[1]]), sheet=compartment, skip=file_list[[1]][[3]]-1) %>%
    dplyr::filter(pct.1 >= alpha_thr & avg_log2FC >= fc_thr & p_val_adj <= pvalue_thr) %>%
    dplyr::rename(Symbol = gene) %>%
    dplyr::mutate(Group = cluster) %>%
    dplyr::group_by(Group) %>%
    dplyr::top_n(topn, p_val_adj) %>%
    dplyr::ungroup()
}

# Collect all cell type data
data_chiji_df <- dplyr::bind_rows(data_chiji_list) %>%
  edit_group_name_my("Group") %>%
  edit_symbol_my("Symbol") %>%
  dplyr::mutate(Group = paste("PDAC.Chiji", Group, sep=".")) %>%
  dplyr::select(Group, Symbol) %>%
  dplyr::ungroup()

# Find repete biomarkers
data_chiji_sum_df <- data_chiji_df %>%
  dplyr::group_by(Symbol) %>%
  dplyr::mutate(Count = dplyr::n_distinct(Group)) %>%
  dplyr::select(-Group) %>%
  unique() %>%
  dplyr::arrange(desc(Count)) %>%
  dplyr::filter(Count > 1)

# Exclude any biomarker that hits more than 2 cell types
data_chiji_df <- data_chiji_df %>%
  dplyr::filter(!Symbol %in% data_chiji_sum_df$Symbol)

print(sort(table(data_chiji_df$Group)))


## Chijimatsu tumor ----


# classical signature (original)	in scRNAseq data	Average Expression in malignant ductal cells	Average Expression in normal ductal cells	Correlation to classical score (r value)	Correlation to basal-like score (r value)	New classical signature
data_classical_df <- readxl::read_excel(paste0(input_path, file_list[[2]][[1]]), sheet="classical", skip=file_list[[2]][[3]]-1) %>%
  dplyr::mutate(Group = "PDAC.Chiji.TumorClassical") %>%
  dplyr::rename(Symbol = `New classical signature`) %>%
  dplyr::filter(Symbol != "NA") %>%
  dplyr::select(Group, Symbol)

# basal-like signature (original)	in scRNAseq data	Average Expression in malignant ductal cells	Average Expression in normal ductal cells	Correlation to classical score (r value)	Correlation to basal-like score (r value)	New classical signature
data_basal_df <- readxl::read_excel(paste0(input_path, file_list[[2]][[1]]), sheet="basal", skip=file_list[[2]][[3]]-1) %>%
  dplyr::mutate(Group = "PDAC.Chiji.TumorBasal") %>%
  dplyr::rename(Symbol = `New classical signature`) %>%
  dplyr::filter(Symbol != "NA") %>%
  dplyr::select(Group, Symbol)


## WRITE SIGNATURES TO FILE ----


data_list <- list(data_chiji_df, data_classical_df, data_basal_df)

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
