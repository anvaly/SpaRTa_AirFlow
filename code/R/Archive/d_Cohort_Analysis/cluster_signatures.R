# Author: Anna Lyubetskaya. Date: 20-10-14


##_ SETUP ENVIRONMENT _##


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_in_out.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_heatmap.R")
source("code/utils/utils_signatures.R")

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_matrix.R")


## PARAMETERS ----


# The cohort of interest regex ID
cohort_name <- "FF_HumanPanc_B"

# DEA thresholds
fc_threshold <- 1
fdr_threshold <- 2  # in -log10(FDR) space

# Gene thresholds
spot_threshold <- 5
sct_threshold <- 1

# Correlation threshold for igraph analysis
corr_threshold <- 0.5


## PATHS ----


# Path to signatures
sig_path <- "C:/Users/lyubetsa/Documents/Data/Signatures/processed_integrated_panc_Oct2020.txt"

# Input folder
input_folder <- "C:/Users/lyubetsa/Documents/Projects_NewTech/10X_ST_Reports/"

# DEA files are tagged as follows
dea_regex_files <- "_dea.csv"
# Seurat RDS files are tagged as follows
data_regex_files <- "_all.Sobj.rds"

# Output folder
output_folder <- "C:/Users/lyubetsa/Documents/Projects_NewTech/10X_ST_Clust_Sig/"


## INGEST DATA ----


# Find all files following the regex
file_list <- dir(input_folder, pattern=paste0(cohort_name, ".*", data_regex_files), full.names=TRUE)

# Ingest clean Seurat data
seurat_list <- list()  
for(rds_file in file_list){
  seurat_list[gsub(paste0("^.+/|", data_regex_files), "", rds_file)] <- read_seurat_rds_my(rds_file, do_subset=TRUE, do_cluster=TRUE, marker_file=NULL)
}

# Read DEA data
dea_df <- read_dir2file_my(input_folder, in_regex=paste0(cohort_name, ".*", dea_regex_files), do_search=TRUE) %>%
  dplyr::mutate(Sample = gsub("^.+/|_dea.csv", "", File)) %>%
  dplyr::select(-File) %>%
  dplyr::filter(avg_logFC >= fc_threshold & p_val_adj_neg_log10 >= fdr_threshold)

# Load signatures
signature_list <- read_filter_signatures_my(sig_path, NULL, sig_length_max=1000, ratio_threshold=0)
signature_list <- signature_list[which(stringr::str_detect(names(signature_list), "CT_"))]
names(signature_list) <- gsub("CT_", "", names(signature_list))


## WRANGLE INPUT DATA ----


# Inverst signatures
gene_sig_df <- invert_list_my(signature_list) %>%
  dplyr::mutate(Sig_Name_Init = Sig_Name,
                Sig_Name = gsub("myeloid_[^;]+\\b", "myeloid", Sig_Name),
                Cat_Name = gsub("_[^;]+\\b", "", Sig_Name)) %>%
  dplyr::group_by(Cat_Name) %>%
  dplyr::mutate(Cat_Name = paste(unique(unlist(strsplit(Cat_Name, ";"))), collapse=";")) %>%
  dplyr::mutate(Num_Cat = stringr::str_count(Cat_Name, ";") + 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Sig_Name = gsub("pdac_|Immune_|Exocrine_|Endocrine_|Stroma_|_\\d+$", "", Sig_Name)) %>%
  dplyr::group_by(Symbol) %>%
  dplyr::mutate(Sig_Name = paste(unique(unlist(strsplit(Sig_Name, ";"))), collapse=";")) %>%
  dplyr::mutate(Num_Sig = stringr::str_count(Sig_Name, ";") + 1) %>%
  dplyr::select(Symbol, InSignature, Num_Cat, Num_Sig, Cat_Name, Sig_Name, Sig_Name_Init)

# Extract long tibbles of SCT values
data_list <- list()
meta_list <- list()

for(f in names(seurat_list)){
  data_list[[f]] <- seurat_expression_to_long_tibble_my(seurat_list[[f]], assay="SCT", 
                                                        spot_threshold=spot_threshold, sct_threshold=sct_threshold, do_stats=TRUE) %>%
    #spot_threshold=1, sct_threshold=-10, do_stats=TRUE) %>%
    dplyr::mutate(Sample = f)
  
  meta_list[[f]] <- tibble::as_tibble(seurat_list[[f]]@meta.data, rownames = "Coordinate") %>%
    dplyr::select(Coordinate, seurat_clusters) %>%
    dplyr::mutate(Sample = f)
}

# Long SCT tibble
data_df <- dplyr::bind_rows(data_list) %>%
  dplyr::mutate(Sample_Spot = paste0(Sample, "_", Coordinate))

# Meta data tibble
meta_df <- dplyr::bind_rows(meta_list)

# Gene properties of interest
gene_df <- data_df %>%
  dplyr::select(Symbol) %>%
  unique %>%
  dplyr::left_join(gene_sig_df %>%
                     dplyr::select(Symbol, InSignature, Num_Cat, Num_Sig, Cat_Name, Sig_Name), by="Symbol") %>%
  dplyr::left_join(dea_df %>% 
                     dplyr::group_by(Symbol) %>%
                     dplyr::summarise(FC = max(avg_logFC)), by="Symbol") %>%
  dplyr::mutate(InDE = abs(FC) > 0) %>%
  tidyr::replace_na(list(InSignature = FALSE, InDE = FALSE, Num_Sigs = 0))


## ANALYZE DATA ----


for(s in names(seurat_list)[1]){
  
  
  ## SUBSET DATA ----
  
  
  # Sample specific expression data
  data_loc_df <- data_df %>%
    dplyr::filter(Sample == s)
  
  # Sample specific meta data
  meta_loc_df <- meta_df %>%
    dplyr::filter(Sample == s)
  
  
  ## PAIRWISE CORRELATION AND STATS ----
  
  
  # Wide SCT tibble
  # Calculate gene pairwise correlations between all genes
  cor_matrix <- cor(data_loc_df %>%
                      df_long2wide_my(rows="Sample_Spot", cols="Symbol", value="SCT") %>%
                      dplyr::select(-Sample_Spot) %>% 
                      as.data.frame())
  
  # Create a long correlation tibble of select genes
  cor_df <- tibble::as_tibble(cor_matrix, rownames="Symbol") %>% 
    df_wide2long_my(key="Symbol2", val="R2") %>%
    dplyr::filter(Symbol != Symbol2)
  
  # Find gene most correlated to each other
  cor_stat_df <- cor_df %>%
    dplyr::group_by(Symbol) %>%
    dplyr::arrange(desc(R2), .by_group = TRUE) %>% 
    dplyr::summarise(gene_max_R2 = dplyr::first(R2), Symbol_closest = dplyr::first(Symbol2))
  
  
  ## FIND GENE CLUSTERS USING A HIERARCHICAL TREE ----
  
  
  # Perform hierarchical clustering
  hclust_obj <- hclust(as.dist(1-cor_matrix)) 
  
  # Cut the hierarchical tree
  hclust_cut_obj <- cutree(hclust_obj, h = mean(hclust_obj$height))
  
  # Size of gene groups detected
  hclust_group_sizes <- table(hclust_cut_obj)
  # Detect sizes with more than 2 genes
  hclust_group_ids <- hclust_group_sizes[which(unname(hclust_group_sizes) > 2)]
  
  # Extract genes that belong to groups with more than 2 genes
  hclust_gene_groups <- sort(hclust_cut_obj[which(unname(hclust_cut_obj) %in% names(hclust_group_ids))])
  
  # Gene groups by hierarchical clustering
  hclust_group_df <- tibble::tibble(Symbol = names(hclust_gene_groups), 
                                    Group = unname(hclust_gene_groups))
  
  # Find mean correlation for each gene pair in a hierarchical group
  hclust_group_stat_df <- cor_df %>% 
    dplyr::filter(Symbol %in% hclust_group_df$Symbol & Symbol2 %in% hclust_group_df$Symbol) %>%
    dplyr::inner_join(hclust_group_df, by="Symbol") %>%
    dplyr::inner_join(hclust_group_df %>%
                        dplyr::rename(Symbol2 = Symbol, Group2 = Group), by="Symbol2") %>%
    dplyr::filter(Group == Group2) %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(tree_group_mean_R2 = mean(R2))
  
  
  cat("Tree height\nMin =", min(hclust_obj$height), "; Max =", max(hclust_obj$height), "; Mean =", mean(hclust_obj$height), "\n")
  cat("Gene number =", length(hclust_gene_groups), "; Group number =", length(hclust_group_ids), "\n")
  cat("Tree cluster correlation: min =", min(hclust_group_stat_df$tree_group_mean_R2), "; max =", max(hclust_group_stat_df$tree_group_mean_R2), "\n")
  
  
  ## FIND GENE CLUSTERS USING A GRAPH ----
  
  
  # Find all positively correlated signatures with non-intersecting gene lists using igraph
  igraph_obj <- igraph::graph_from_data_frame(data.frame(cor_df %>%
                                                           dplyr::filter(R2 >= corr_threshold)), directed=FALSE, vertices=NULL)
  
  # igraph signature group lists
  igraph_list <- igraph::clusters(igraph_obj)$membership
  
  # igraph signature group tibble
  igraph_group_df <- tibble::tibble(Symbol = names(igraph_list), 
                                    Group_igraph = paste0("igraph", unname(igraph_list)))
  
  cat("igraph group sizes:", table(igraph_group_df$Group_igraph), "\n")
  
  
  ## INTEGRATE CORRELATION DATA ----
  
  
  # Combine pairwise correlation with annotation information
  integrated_properties_st_df <- hclust_group_df %>%
    dplyr::left_join(igraph_group_df, by="Symbol") %>%
    dplyr::inner_join(hclust_group_stat_df, by="Group") %>%
    dplyr::inner_join(cor_stat_df, by="Symbol") %>%
    dplyr::rename(Group_Tree = Group) %>%
    dplyr::mutate_if(is.numeric, round, 3) %>%
    dplyr::inner_join(gene_df %>%
                        dplyr::select(Symbol, InDE, Num_Cat, Num_Sig, Cat_Name, Sig_Name), by="Symbol") %>%
    dplyr::arrange(desc(tree_group_mean_R2), Group_Tree, Group_igraph, Symbol)
  
  
  # List of genes in signatures with their ST properties
  integrated_properties_sig_df <- gene_sig_df %>%
    dplyr::left_join(gene_df %>%
                       dplyr::mutate(InST = TRUE) %>%
                       dplyr::select(Symbol, InST), by="Symbol") %>%
    dplyr::left_join(integrated_properties_st_df, by=c("Symbol", "Num_Cat", "Num_Sig", "Cat_Name", "Sig_Name")) %>%
    dplyr::mutate(InSTCluster = Group_Tree >= 0) %>%
    dplyr::select(Symbol, Num_Cat, Num_Sig, Cat_Name, Sig_Name, InST, InSTCluster, InDE, Group_Tree, Group_igraph, Symbol_closest, gene_max_R2) %>%
    dplyr::arrange(Cat_Name, Sig_Name)
  
  
  ## WRITE TO FILE ----
  
  
  # Write to file
  filename <- paste0(output_folder, cohort_name, "_", s, "_st_summary.txt")
  readr::write_delim(integrated_properties_st_df, path=filename, delim = "\t", append=FALSE, col_names = TRUE)
  
  # Write to file
  filename <- paste0(output_folder, cohort_name, "_", s, "_sig_summary.txt")
  readr::write_delim(integrated_properties_sig_df, path=filename, delim = "\t", append=FALSE, col_names = TRUE)
  
  # Write to file
  filename <- paste0(output_folder, cohort_name, "_", s, "_cor.txt")
  readr::write_delim(cor_df %>%
                       dplyr::mutate_if(is.numeric, round, 3), path=filename, delim = "\t", append=FALSE, col_names = TRUE)
  
  
  ## CLUSTER DATA ----
  
  
  # Heatmap parameters
  params <- list(cell_value = "SCT_zscore",
                 row_label = "Symbol", 
                 col_label = "Coordinate", 
                 distance = "pearson",
                 row_annotation = c("InSignature", "InDE"),
                 col_annotation = c("seurat_clusters"),
                 range = c(-1, 0, 1),
                 colors = c("red4", "white", "blue4"))
  
  # Create heatmap of sample clusters using DE genes
  filename <- paste0(output_folder, cohort_name, "_", s, "_igraph_clusters.png")
  hm <- create_heatmap_my(data_loc_df, params, row_list = igraph_group_df$Symbol,
                          row_meta_df = gene_df, col_meta_df = meta_loc_df, filename=filename)
  
  # Create heatmap of sample clusters using DE genes
  filename <- paste0(output_folder, cohort_name, "_", s, "_tree_sig_clusters.png")
  hm <- create_heatmap_my(data_loc_df, params, row_list = integrated_properties_sig_df %>% dplyr::filter(InST == TRUE) %>% dplyr::pull(Symbol),
                          row_meta_df = gene_df, col_meta_df = meta_loc_df, filename=filename)
  
}
