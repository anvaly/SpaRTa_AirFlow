# Author: Anna Lyubetskaya. Date: 20-08-06

## Analyze mouse ARCHS4 data


## ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_archs4.R")
source("code/utils/utils_pca.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_rna_diff_expr.R")
source("code/utils/utils_ggplot.R")
source("code/utils/utils_heatmap.R")
source("code/utils/utils_signatures.R")


## PARAMETERS ----


# Signatures to select
sig_select <- c("PDAC.U.Immune.Macrophage", "BMS.WZ.CD8.effector.CuiCell2021",
                "BMS.Pathway.IFNg", "BMS.Pathway.IFNa", "BMS.MoCR.PTPN2", "BMS.MoCR.ADAR",
                "Syng.U.Tumor.MC38", "Syng.K18.B16F10", "Syng.K18.CT26", "Syng.K18.EMT6", "Syng.K18.LL2")

# Use the following signatures to filter samples
sig_filter <- NULL  # c("Syng.U.Tumor.MC38", "Syng.K18.B16F10")

# Rename signatures for visualizations
sig_rename <- c("Macrophages", "CD8", "IFNg sig", "IFNa sig", "PTPN2", "ADAR", "MC38", "B16F10", "CT26", "EMT6", "LL2")

# List of cell lines of interest
cell_line_list <- c("4T1", "CT26", "MC38", "B16", "LL2", "TRAMPC2", "KPC", "1956", "EL4", "A20", "EMT6", "Renca")


## PATHS ----


# Path to a signature file
sig_path <- "data/import/Signatures/signatures_220210.txt"

# Output folder
output_path <- "XXXX"


## INGEST ARCHS4 ----


# Setup and read ARCHS4 data
data_archs4 <- setup_read_archs4(cell_line_list=cell_line_list, group="all")


## NORMALIZE ARCHS4 ----


hist(log10(data_archs4$samples$readsaligned))

# Select samples with at least 10M aligned reads
samples_select <- data_archs4$samples$GEO_ID[which(data_archs4$samples$readsaligned >= 1000000 & 
                                                     data_archs4$samples$readsaligned <= 100000000)]

# Calculate CPM values
data_eset <- Biobase::ExpressionSet(as.matrix(data_archs4$expression[, c("Symbol", samples_select)] %>% 
                                                     tibble::column_to_rownames("Symbol")))

# Normalize expression data
data_norm_eset <- edgeR_normalize_my(data_eset)
data_norm_matrix <- Biobase::exprs(data_norm_eset)

# Find samples with suspiciously few well expressed genes
sample_list <- names(which(colSums((data_norm_matrix > 0)) >= 10000))
data_norm_matrix <- data_norm_matrix[, sample_list]

# Find genes with good expression levels
genes_represented <- names(which(rowSums(data_norm_matrix >= 3) >= 10))

# Normalized data wide tibble
data_norm_wide_df <- tibble::as_tibble(data_norm_matrix, rownames="Symbol") %>%
  dplyr::filter(Symbol %in% genes_represented)

# Normalized data long tibble
data_norm_df <- data_norm_wide_df %>%
  df_wide2long_my(key="Sample_ID", val="Expression")

# Meta data
meta_df <- data_archs4$samples %>%
  dplyr::filter(GEO_ID %in% sample_list) %>%
  dplyr::rename(Sample_ID = GEO_ID)

# Calculate z-scores for each gene
zscore_df <- df_zscore_my(data_norm_df, col_by="Symbol", value="Expression")

hist(log10(meta_df$readsaligned))


## INGEST AND WRANGLE SIGNATURES ----


# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, genes_represented, sig_names=sig_select)

# Calculate gene signature scores
sig_val_list <- list()
for(s in names(signature_list)){
  sig_val_list[[s]] <- zscore_df %>%
    dplyr::filter(Symbol %in% signature_list[[s]]) %>%
    dplyr::group_by(Sample_ID) %>%
    dplyr::summarise(Score = round(median(Expression_zscore), 3)) %>%
    dplyr::mutate(Signature = s)
}

# Create a single table of singature scores
sig_val_df <- dplyr::bind_rows(sig_val_list) %>%
  dplyr::inner_join(meta_df %>%
                      dplyr::select(Sample_ID, Cancer, cell_line), by="Sample_ID")

# Filter samples by select user-defined signatures
if(!is.null(sig_filter)){
  samples_sig_select <- sig_val_df %>% 
    dplyr::filter(Signature %in% sig_filter & Score >= 1) %>% 
    dplyr::pull(Sample_ID) %>% 
    unique()
} else{
  samples_sig_select <- unique(sig_val_df$Sample_ID)
}

# Parse the meta_data
meta_select <- meta_df %>%
  dplyr::filter(grepl("untreated|control|igg|pbs", tolower(characteristics_ch1))) %>%
  dplyr::pull(Sample_ID)

# Filter data down to samples with strong signature scores
sig_val_df <- sig_val_df %>%
  dplyr::filter(Sample_ID %in% intersect(samples_sig_select, meta_select))

# Long to wide tibble of signature scores with meta data
sig_val_wide_df <- sig_val_df %>%
  df_long2wide_my(rows="Sample_ID", cols="Signature", value="Score") %>%
  dplyr::inner_join(meta_df, by="Sample_ID")

# Print data to file
filename <- paste0(output_path, "signature_table.txt")
readr::write_delim(sig_val_wide_df, filename, delim="\t")


## VISUALIZE SIGNATURE SCORES ----


# Rename signatures for visualizations
if(!is.null(sig_rename)){
  names(sig_rename) <- sig_select
  print(sig_rename)
  
  sig_val_df[["Signature"]] <- sig_rename[sig_val_df$Signature]
}

# Create a box plot of gene signature levels by experiment group, fixed axes
filename <- paste0(output_path, "signature_boxplot_fixed")
create_violin_plot_my(sig_val_df, x_label="cell_line", y_label="Score", fill_label="Cancer", 
                   facet_var=c("Signature", "fixed"), filename=filename, labels=NULL)

# Create a box plot of gene signature levels by experiment group, free axes
filename <- paste0(output_path, "signature_boxplot_free")
create_violin_plot_my(sig_val_df, x_label="cell_line", y_label="Score", fill_label="Cancer", 
                   facet_var=c("Signature", "free"), filename=filename, labels=NULL)


# Create a heatmap of signature scores by sample
params <- list(cell_value = "Score",
               row_label = "Signature", 
               col_label = "Sample_ID", 
               distance = "pearson",
               row_annotation = NULL,
               col_annotation = c("cell_line"),
               range = c(-2, 0, 2),
               colors = c("royalblue4", "white", "red3"),
               row_order = NULL,
               col_order = NULL)

filename <- paste0(output_path, "signature_heatmap.png")
create_heatmap_my(sig_val_df, params, row_list=NULL, col_list=NULL, col_meta_df=meta_df, row_meta_df=NULL, filename=filename)


## PCA ----


params <- list(sample_value = "Sample_ID", sample_filter = data_norm_df$Sample_ID,
               feature_value = "Symbol", feature_filter = data_norm_df$Symbol,
               cell_value = "Expression", color_by = c("cell_line"))

filename <- paste0(output_folder, "pca_cell_line")
# Perform first round of PCA
pc1_df <- pca_from_long_tibble_my(data_norm_df, meta_df, params, filename=filename, table2file=FALSE, extra_plot=TRUE)

# Sample extreme outlier sample using first PCA
samples_select <- pc1_df %>% 
  dplyr::filter(.fittedPC2 < 1000) %>%
  dplyr::pull(Sample_ID)


## ARCHS4 DEA ----


# Library size mean and SD
reads_mean <- mean(data_archs4$samples$reads_aligned)
reads_sd <- sd(data_archs4$samples$reads_aligned)

# Select well represented genes for the PCA approved list of samples
genes_represented <- names(which(rowSums(data_norm_matrix[,samples_select] >= 5) >= 20))

# Select samples with comparable library sizes
#samples_select <- data_archs4$samples %>%
#  dplyr::filter(reads_aligned >= reads_mean - reads_sd &
#                  reads_aligned <= reads_mean + reads_sd) %>%
#  dplyr::pull(GEO_ID)

# DEA parameters
params_dea <- list(sample = "GEO_ID", sample_filter = samples_select,
                   feature = "Symbol", feature_filter = genes_represented,
                   batch = "Cancer", condition = "cell_line",
                   contrasts = list(c("cell_line", "B16", "MC38")),
                   num_cores = 1)

# DEA to select biomarkers
dea_list <- dea_from_wide_tibble_my(data_archs4$expression, data_archs4$samples, params_dea, 
                                  name_out=paste0(output_folder, "/dea"))

top_n <- c(25, 50, 75, 100)

# Select top n for every contrast
for(contrast in names(dea_list$contrasts)){
  gene_up_list <- dea_list$contrasts[[contrast]] %>% 
    dplyr::filter(IsSignificant == TRUE & log2FoldChange >= 2) %>% 
    dplyr::filter(!grepl("^GM\\d+$", Symbol)) %>%
    dplyr::arrange(desc(IsSignificant), desc(abs(log2FoldChange)))
  
  gene_dn_list <- dea_list$contrasts[[contrast]] %>% 
    dplyr::filter(IsSignificant == TRUE & log2FoldChange <= -2) %>% 
    dplyr::filter(!grepl("^GM\\d+$", Symbol)) %>%
    dplyr::arrange(desc(IsSignificant), desc(abs(log2FoldChange)))
  
  for(n in top_n){    
    filename <- paste0(output_folder, "/dea_tops.txt")
    write(paste0(contrast, "_UP_", n), file=filename, append=TRUE)
    write(paste(sort(gene_up_list$Symbol[1:n]), collapse=","), file=filename, append=TRUE)
    
    write(paste0(contrast, "_DN_", n), file=filename, append=TRUE)
    write(paste(sort(gene_dn_list$Symbol[1:n]), collapse=","), file=filename, append=TRUE)
  }
}
