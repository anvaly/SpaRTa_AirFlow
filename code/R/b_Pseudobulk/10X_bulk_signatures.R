# Author: Anna Lyubetskaya. Date: 22-02-16

# This script:
# - Calculates pseudo-bulk expression profiles for a set of ST samples
# - Calculates select gene signatures for each sample
# - Compares samples as though they are a bulk cohort


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

library(Biobase)
library(Seurat)

source("code/utils/utils_ggplot.R")
source("code/utils/utils_heatmap.R")
source("code/utils/utils_signatures.R")
source("code/utils/utils_rna_diff_expr.R")

source("code/R/Utils/utils_10X_matrix.R")


## PARAMETERS ----


# Run name to name the output folder
run_name <- "FERMT"

# Seurat RDS files are tagged as follows
cohort_name <- "PDAC"
# Samples to include
cohort_regex <- "PDAC"

# Break data by pathology
break_path_classes <- TRUE


# Exclude certain samples from analysis
samples_exclude <- NULL

# CPM thresholds: standard is 1 and 3. When the raw data has too many zeros, these thresholds need fine-tuning due to CPM normalization
cpm_threshold <- 0.5
rep_threshold <- 5

# Signatures to select
sig_select <- c("PDAC.collisson.classical", "PDAC.moffitt.classical", #"PDAC.P19.Ductal_2",
                "PDAC.moffitt.basal", "PDAC.collisson.quasimesenchymal", 
                "PDAC.moffitt.activatedstroma", #"PDAC.P19.Fibroblast", 
                "PDAC.moffitt.normalstroma", #"PDAC.P19.Stellate", 
                "BMS.PDAC.Hypoxia.CL", "Syng.U.State.Hypoxia.metabolism", #"PDAC.P19.Endothelial", 
                #"PDAC.P19.Acinar", "PDAC.P19.Endocrine", "PDAC.P19.Ductal_1", 
                #"PDAC.P19.Bcell", "PDAC.P19.Macrophage", "PDAC.P19.Tcell",
                "BMS.Pathway.TNFa", "BMS.Pathway.TGFB", "BMS.Pathway.IFNg", "BMS.Pathway.IFNa", 
                "Oxford.Tcell_exhausted",
                #"PDAC.CosMx.Bcell", "PDAC.CosMx.TumorBasal",
                #"PDAC.CosMx.TumorClassical", "PDAC.CosMx.TumorClassicalCycling",
                #"PDAC.CosMx.Ductal", "PDAC.CosMx.Endocrine",
                #"PDAC.CosMx.Endothelial", "PDAC.CosMx.Exocrine",
                #"PDAC.CosMx.Fibroblast", "PDAC.CosMx.Lowquality",
                #"PDAC.CosMx.Macrophage", "PDAC.CosMx.Mast",
                #"PDAC.CosMx.MDSC", "PDAC.CosMx.TumorMixed",
                #"PDAC.CosMx.TumorOther", "PDAC.CosMx.Plasma", "PDAC.CosMx.Tcell",
                "PDAC.Chiji.Acinar", "PDAC.Chiji.Bcell", "PDAC.Chiji.Ductal_1",
                "PDAC.Chiji.Ductal_2", "PDAC.Chiji.Endocrine", "PDAC.Chiji.Endothelial",
                "PDAC.Chiji.Fibroblast", "PDAC.Chiji.Macrophage", "PDAC.Chiji.Stellate",
                "PDAC.Chiji.Tcell", "PDAC.Chiji.TumorClassical", "PDAC.Chiji.TumorBasal"
)


# Rename signatures for visualizations
sig_rename <- c("collisson.classical", "moffitt.classical", #"P19.Ductal_2",
                "moffitt.basal", "collisson.quasimesenchymal", 
                "moffitt.activatedstroma", #"P19.Fibroblast", 
                "moffitt.normalstroma", #"P19.Stellate", 
                "Hypoxia.CL", "Hypoxia.metabolism", #"P19.Endothelial", 
                #"P19.Acinar", "P19.Endocrine", "P19.Ductal_1", 
                #"P19.Bcell", "P19.Macrophage", "P19.Tcell",
                "Pathway.TNFa", "Pathway.TGFB", "Pathway.IFNg", "Pathway.IFNa", 
                "Tcell_exhausted",
                #"CosMx.Bcell", "CosMx.TumorBasal",
                #"CosMx.TumorClassical", "CosMx.TumorClassicalCycling",
                #"CosMx.Ductal", "CosMx.Endocrine",
                #"CosMx.Endothelial", "CosMx.Exocrine",
                #"CosMx.Fibroblast", "CosMx.Lowquality",
                #"CosMx.Macrophage", "CosMx.Mast",
                #"CosMx.MDSC", "CosMx.TumorMixed",
                #"CosMx.TumorOther", "CosMx.Plasma", "CosMx.Tcell",
                "Chiji.Acinar", "Chiji.Bcell", "Chiji.Ductal_1",
                "Chiji.Ductal_2", "Chiji.Endocrine", "Chiji.Endothelial",
                "Chiji.Fibroblast", "Chiji.Macrophage", "Chiji.Stellate",
                "Chiji.Tcell", "Chiji.TumorClassical", "Chiji.TumorBasal"
)

# Meta data to sum across columns and add to sample-level meta data
meta_select <- paste0("Pathology.", c("IntestineAdj", "MuscleAdj", "NormalAdj",  "LymphNode", "TLSAggregate", "TLSImmature", "TLSMature", 
                                      "Adipose", "Vessel", "Muscle", "Nerve", "NonEpi", "BenignEpi", "Blood", "ExoEndo", "LuminalNec", "Tumor"), ".percent")

# Rename meta data for vis
meta_rename <- c("IntestineAdj", "MuscleAdj", "NormalAdj",  "LymphNode", "TLSAggregate", "TLSImmature", "TLSMature", 
                 "Adipose", "Vessel", "Muscle", "Nerve", "NonEpi", "BenignEpi", "Blood", "ExoEndo", "LuminalNec", "Tumor")

# Enforce upper threshold across meta columns in the heatmap
meta_ceiling <- 50  # NULL

# Add one-off genes to the plots
gene_select <- c("FERMT1", "FERMT2", "FHL2")


## ADJUST PARAMETERS ----


# Set absent signature names
if(is.null(sig_rename)){
  sig_rename <- sig_select
}

length(sig_select) == length(sig_rename)

# Create a dictionary to rename signatures
sig_dict <- c(sig_rename, gene_select)
names(sig_dict) <- c(sig_select, gene_select)


# Set absent meta fields names
if(is.null(meta_rename)){
  meta_rename <- meta_select
}

length(meta_select) == length(meta_rename)

# Create a dictionary to rename meta data
meta_dict <- meta_rename
names(meta_dict) <- meta_select


## PATHS ----


# Path to a signature file
# Signature_name	Gene_list
sig_path <- "data/import/Signatures/signatures_pdac_20230720.txt"

# Input folder
input_paths <- c("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_name, "_", run_name, "/")

dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Find all files following the regex
file_list <- unname(unlist(sapply(input_paths, function(p) dir(p, pattern=paste0(cohort_regex, ".*.rds"), full.names=TRUE))))

# Identify samples user wants to exclude and remove them from the file list
if(!is.null(samples_exclude)){
  samples_exclude_indx <- c(unname(sapply(samples_exclude, function(x) grep(x, file_list))))
  file_list <- setdiff(file_list, file_list[samples_exclude_indx])
}

# Ingest RDS objects with Seurat data and tranform them into a pseudo-count vector
pseudo_count_list <- list()
meta_data_list <- list()
for(rds_file in file_list){
  
  # Load the dataset
  data_seurat <- readRDS(rds_file)
  
  # Sample name
  sample_name <- data_seurat@misc[["user.Sample_Name"]]
  

  # Find major pathology classes in the sample
  if(break_path_classes == TRUE){
    path_group_counts <- table(data_seurat@meta.data$Pathology.Group)
    
    # Note: Only groups with enough spots will be included
    path_groups <- names(path_group_counts[which(path_group_counts >= 100)])
  } else{
    path_groups <- c("skip")
  }
  
  
  # Calculate pseudo bulk counts for a 10X sample
  for(p in path_groups){
    
    if(p == "skip"){
      
      # Sum spot level meta data to sample level meta data
      meta_spot_list <- Matrix::colSums(data_seurat@meta.data[intersect(colnames(data_seurat@meta.data), meta_select)]) / nrow(data_seurat@meta.data)
      names(meta_spot_list) <- meta_dict[names(meta_spot_list)]
      
      # Extract meta data information
      meta_data_list[[sample_name]] <- c(data_seurat@misc, as.list(meta_spot_list))
      
      pseudo_count_list[[sample_name]] <- pseudo_bulk_counts_my(data_seurat, log_transform=FALSE, assay="Spatial", slot="data")
      
    } else{
      
      # Select spots corresponding to a specific pathology compartment
      data_seurat_subset <- subset(data_seurat, cells=rownames(data_seurat@meta.data)[which(data_seurat@meta.data$Pathology.Group == p)])
      
      # Sum spot level meta data to sample level meta data
      meta_spot_list <- Matrix::colSums(data_seurat_subset@meta.data[intersect(colnames(data_seurat_subset@meta.data), meta_select)]) / nrow(data_seurat_subset@meta.data)
      names(meta_spot_list) <- meta_dict[names(meta_spot_list)]
      
      # Extract meta data information
      meta_data_list[[paste(sample_name, p)]] <- c(data_seurat_subset@misc, as.list(meta_spot_list))
      
      pseudo_count_list[[paste(sample_name, p)]] <- pseudo_bulk_counts_my(data_seurat_subset, log_transform=FALSE, assay="Spatial", slot="data")
      
    }
    
  }
}


## WRANGLE DATA ----


# Create a meta data tibble
meta_df <- list2tibble_my(meta_data_list)
meta_df$user.Sample_Name <- names(meta_data_list)

# Adjust meta data if 
# Combine a named list of lists into an expression wide tibble
data_wide_df <- pseudo_bulk_vectors_combine_my(pseudo_count_list, names(pseudo_count_list)) %>%
  dplyr::mutate_at(vars(dplyr::all_of(names(pseudo_count_list))), as.numeric)


## CREATE AND WRANGLE ESET ----


# Create an annotated data frame object of gene annotations
featureData <- new("AnnotatedDataFrame", data = data.frame(Symbol = data_wide_df$Symbol,
                                                           row.names = data_wide_df$Symbol))

# Create an annotated data frame object of sample phenotypes
phenoData <- new("AnnotatedDataFrame", data = meta_df %>% 
                   tibble::column_to_rownames("user.Sample_Name") %>%
                   as.data.frame())

# Create an eSet object
data_eset <- Biobase::ExpressionSet(assayData = as.matrix(data_wide_df %>% 
                                                            tibble::column_to_rownames("Symbol")),
                                    phenoData = phenoData,
                                    featureData = featureData)


# Calculate CPM values
data_norm_eset <- edgeR_normalize_my(data_eset)

# Extract the CPM matrix from the eSet
data_norm_matrix <- Biobase::exprs(data_norm_eset)

# Identify genes in the eset that pass the CPM and replicate thresholds
genes_present <- Reduce(intersect, sapply(pseudo_count_list, function(x) names(x)))
genes_represented <- names(which(rowSums(data_norm_matrix >= cpm_threshold) >= rep_threshold))
gene_list <- intersect(genes_present, genes_represented)

# Create a long tibble of expression CPMs
expr_df <- tibble::as_tibble(data_norm_matrix, rownames="Symbol") %>%
  df_wide2long_my(key="Sample", val="CPM")

# Calculate z-scores for each gene
zscore_df <- df_zscore_my(expr_df, col_by="Symbol", value="CPM")


## INGEST AND WRANGLE SIGNATURES ----


# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, genes_represented, sig_names=sig_select)

# Check missing signatures
setdiff(sig_select, names(signature_list))


# Calculate gene signature scores
sig_val_list <- list()
for(s in names(signature_list)){
  sig_val_list[[s]] <- zscore_df %>%
    dplyr::filter(Symbol %in% signature_list[[s]]) %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarise(Score = round(median(CPM_zscore), 3)) %>%
    dplyr::mutate(Signature = s)
}

# Add individual genes if user-defined
if(!is.null(gene_select)){
  for(s in gene_select){
    sig_val_list[[s]] <- zscore_df %>%
      dplyr::filter(Symbol == s) %>%
      dplyr::group_by(Sample) %>%
      dplyr::summarise(Score = round(median(CPM_zscore), 3)) %>%
      dplyr::mutate(Signature = s)
  }
}


# Create a single table of signature scores
sig_val_df <- dplyr::bind_rows(sig_val_list) %>%
  dplyr::inner_join(meta_df, by=c("Sample" = "user.Sample_Name"))

# Print data to file
filename <- paste0(output_path, "table_", cohort_name, "_", run_name, ".txt")
readr::write_delim(sig_val_df, filename, delim="\t")


## VISUALIZE SIGNATURE SCORES ----


# Rename signatures for visualizations
sig_val_df[["Signature"]] <- unname(sig_dict[sig_val_df$Signature])


# Create a bar and a box plot of gene signature levels by experiment group
for(s in unique(sig_val_df$Signature)){
  
  # Wrangle data for vis
  sig_loc_df <- sig_val_df %>%
    dplyr::filter(Signature == s) %>%
    dplyr::mutate(Pathology = gsub(".+ ", "", Sample))
  
  p <- create_bar_plot_my(sig_loc_df, 
                          x_label="Sample", y_label="Score", fill_label="Pathology", 
                          filename=NULL, labels=c("Sample", paste(s, "signature/gene score"), paste(s, "pseudobulk gene expression")), reorder_x=TRUE)
  
  filename <- paste0(output_path, "bar_", s, "_", cohort_name, "_", run_name)
  write_plot2file_my(p, filename, num_row=1, num_col=ceiling(length(unique(sig_val_df$Sample))/100))
  
  if(break_path_classes == TRUE){
    filename <- paste0(output_path, "box_", s, "_", cohort_name, "_", run_name)
    p <- create_box_plot_my(sig_loc_df, 
                            x_label="Pathology", y_label="Score", fill_label="Pathology", 
                            filename=filename, labels=c("Pathology niche", paste(s, "signature/gene score"), paste(s, "pseudobulk gene expression")), reorder_x=FALSE)
  }
}


# Adjust sample name column
meta_df <- meta_df %>%
  dplyr::rename(Sample = user.Sample_Name)

# Set an upper ceiling on epi fields
if(!is.null(meta_ceiling)){
  for(m in names(meta_dict)){
    indx <- which(meta_df[[meta_dict[m]]] > meta_ceiling)
    if(length(indx) > 0){
      meta_df[indx, meta_dict[m]] <- meta_ceiling
    }
  }
}

# Round pathology annotation columns
# Remove NAs
for(m in names(meta_dict)){
  if(grepl("Pathology.*percent", m)){
    meta_df[meta_dict[m]] <- round(meta_df[[meta_dict[m]]])
    
    indx <- is.na(meta_df[[meta_dict[m]]])
    meta_df[indx, meta_dict[m]] <- 0
  }
}

# Create a heatmap of signature scores by sample
params <- list(cell_value = "Score",
               row_label = "Signature", 
               col_label = "Sample", 
               distance = "pearson",
               row_annotation = NULL,
               col_annotation = intersect(c(meta_select, meta_rename), colnames(meta_df)),
               range = c(-2, 0, 2),
               colors = c("royalblue4", "white", "red3"),
               show_column_dend = TRUE,
               show_row_dend = TRUE)

# Order the heatmap by a pathology factor
# column_order = order(meta_df[[meta_rename]], meta_df$user.Sample_ID)

filename <- paste0(output_path, "hm_", cohort_name, "_", run_name, ".png")
create_heatmap_my(sig_val_df, params, row_list=NULL, col_list=NULL, col_meta_df=meta_df, row_meta_df=NULL, filename=filename,
                  width=20, height=10)
