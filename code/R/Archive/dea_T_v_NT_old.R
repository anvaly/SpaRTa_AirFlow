# Author: Anna Lyubetskaya. Date: 21-07-20
# Perform DEA of spots with a set of properties to the rest
# Spots are selected using signature scores, pathology, and clusters
# In this setting, the more data the better


## ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_signatures.R")

source("code/R/b_Cluster_Biomarkers/Seurat_cluster_utils.R")


## PARAMETERS ----


# Path to processed Seurat data
#sample_name <- "HumanPanc_merge"
#sample_exclude <- c("HumanPanc_ROI2_FFPE_D_Dec20_T")
sample_name <- "pdac_cohort_spotclean_harmony_integrated"
sample_exclude <- NULL

# Name for the signature group being plotted
sig_select <- "collisson.pdac.classical"
sig_pval_threshold <- NULL  # 0.0005

# Name of the pathology field to filter by groups
#pathology_select <- "Pathology.Epithelium.percent"
pathology_select <- "nCount_Spatial"
pathology_threshold <- 50

# # List of clusters to select - assumes there is a Clustering_Preferred column in the object to do this
# clusters_select <- tibble::tibble(user.Sample_Name = c("HumanPanc_FF_A_Jun20_F", "HumanPanc_FF_B_Jun20_T", "HumanPanc_FF_C_Jun20_F", "HumanPanc_FF_D_Jun20_F",
#                                                        "HumanPanc_FF_A_Jun20_F", "HumanPanc_FF_B_Jun20_T", "HumanPanc_FF_C_Jun20_F", "HumanPanc_FF_D_Jun20_F",
#                                                        "HumanPanc_ROI1_FFPE_A_Apr21", "HumanPanc_ROI1_FFPE_A_Jun20_F", "HumanPanc_ROI1_FFPE_B_Dec20_T", 
#                                                        "HumanPanc_ROI1_FFPE_B_Jun20_F", "HumanPanc_ROI1_FFPE_C_Jun20_F", "HumanPanc_ROI1_FFPE_D_Jun20_F", 
#                                                        "HumanPanc_S001890_FFPE_C_Apr21", "HumanPanc_S001891_FFPE_D_Apr21"),
#                                   Clustering_Preferred = c("0", "0", "0", "0", "1", "1", "1", "1", "4", "4", "4", "4", "4", "4", "2", "1"),
#                                   Cluster_Select = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
clusters_select <- NULL

# SCT and spot number thresholds for gene filtering
sct_threshold <- 0.5
spot_threshold <- 100

# Assay and slot to use to calculate signature score
assay <- "SCT"
slot <- "data"

# Reference signatures to use for plotting and annotation
sig_reference <- c("bailey.pdac.progenitor", "bailey.pdac.squamous",  # "bailey.pdac.adex", "bailey.pdac.immunogenic", 
                   "collisson.pdac.classical", "collisson.pdac.quasimesenchymal",  # "collisson.pdac.exocrine",
                   "moffitt.pdac.classical", "moffitt.pdac.basal")

# Feature # in a spot threshold
feature_threshold <- 3000


## PATHS ----


# Path to a signature file
sig_path <- "data/import/Signatures/signatures_pdac_t1000_May21.txt"

# Location of pre-processed data
# input_path <- paste0("XXXX")
input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, sample_name, "_", sig_select, "_3K/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Seurat data
data_seurat <- readRDS(input_path)

# Remove samples if necessary
if(!is.null(sample_exclude)){
  barcode_list <- rownames(data_seurat@meta.data[which(!data_seurat@meta.data[["user.Sample_Name"]] %in% sample_exclude),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# Remove spots if necessary
if(!is.null(feature_threshold)){
  barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data$nFeature_SCT >= feature_threshold),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# Cleanup the old signature data if necessary
data_seurat@meta.data <- data_seurat@meta.data[, !grepl("sig.", colnames(data_seurat@meta.data))]


## CALCULATE TARGET SIGNATURE SCORES ----


## Calculate signature scores and add them to Seurat meta data

# Find genes abundant in this sample
# gene_list <- unique(seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold, 
                                                    # assay=assay, slot=slot, split_by="user.Sample_Name"))
gene_list <- data_seurat@assays$SCT@var.features
intersect(c("CEACAM5", "MIF"), gene_list)
# Load signatures and filter them down to only well represented genes
signature_init_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_reference,
                                                 sig_length_min=1, sig_length_max=1000, ratio_threshold=0)

# Find all genes in signatures and their signature assignations
sig_invert_df <- invert_list_my(signature_init_list)

# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=c(sig_select), sig_length_min=1, sig_length_max=1000, ratio_threshold=0)
cat("Target signature length =", length(signature_list[[sig_select]]))

# Add signature scores to a seurat object
data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.", assay=assay)

# Find column names in the Seurat object - Seurat can rename original signatures
sig_names <- colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))])

# Find the properties of the genes in the target signature
sig_expr_df <- seurat_gene_properties_my(data_seurat, signature_list[[1]], assay="SCT", slot="data", sct_threshold=sct_threshold)

# Write stats to file
filename <- paste0(output_path, "/table_sig_by_gene_", sample_name, "_", sig_names, "_", pathology_select, ".txt")
readr::write_delim(sig_expr_df, filename, delim="\t")


## CREATE A RANDOM SIGNATURE SCORE BACKGROUND ----


if(!is.null(sig_pval_threshold)){
  # Write/read random gene expression scores to file
  filename_random <- paste0(output_path, "/table_sig_random_", sample_name, "_", sig_names, "_", pathology_select, ".txt")
  
  # Create a random signature background distribution and score the target signature against it
  data_seurat <- signature_empirical_pvalue_my(data_seurat, signature_list[[1]], gene_list, output_path, sample_name, sig_names[1], sig_random_filename=filename_random,
                                               n_simulations=10000, sct_threshold=sct_threshold, spot_threshold=spot_threshold, assay=assay, slot=data)
} else{
  # Dummy empirical p-values with just the signature scores
  data_seurat@meta.data[["Sig_Emp_Pvalue"]] <- data_seurat@meta.data[[sig_names[1]]]
}


## WRANGLE META DATA ----


# Extract coordinate meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

# Add user.Sample_Name to meta data if absent: happens in stand-alone samples
if(!"user.Sample_Name" %in% colnames(meta_df)){
  meta_df["user.Sample_Name"] = data_seurat@misc$user.Sample_Name
}


## IDENTIFY SIGNATURE THRESHOLDS ----


# Global signature stats
sig_median_global <- median(meta_df[[sig_names]]) 
sig_sd_global <- sd(meta_df[[sig_names]]) 
sig_threshold <- sig_median_global + sig_sd_global

# Calculate signature mean and SD for each sample in the cohort
meta_sig_df <- meta_df %>%
  dplyr::select(dplyr::all_of(c("user.Sample_Name", sig_names, pathology_select))) %>%
  dplyr::group_by(user.Sample_Name) %>%
  dplyr::summarise(Sig_Mean = mean(!!rlang::sym(sig_names)),
                   Sig_SD = sd(!!rlang::sym(sig_names)),
                   Path_Mean = median(!!rlang::sym(pathology_select)),
                   Path_SD = sd(!!rlang::sym(pathology_select)),
                   Spot_Count = sum(!!rlang::sym(sig_names) >= -1000),
                   Sig_Threshold = Sig_Mean + Sig_SD,
                   Spot_Count_IsHigh = sum(!!rlang::sym(sig_names) >= Sig_Threshold),
                   Spot_Perc_IsHigh = round(Spot_Count_IsHigh / Spot_Count * 100))

# Write stats to file
filename <- paste0(output_path, "/table_sig_by_sample_", sample_name, "_", sig_names, "_", pathology_select, ".txt")
readr::write_delim(meta_sig_df, filename, delim="\t")

# Add preferred clusters if user-defined
if(!is.null(clusters_select)){
  # Defactor column
  meta_df["Clustering_Preferred"] <- as.character(meta_df[["Clustering_Preferred"]])
  
  meta_df <- meta_df %>%
    dplyr::left_join(clusters_select, by=c("user.Sample_Name", "Clustering_Preferred")) %>%
    dplyr::mutate(Sig_IsHigh = Sig_Emp_Pvalue <= sig_pval_threshold,
                  Path_IsHigh = !!rlang::sym(pathology_select) >= pathology_threshold,
                  Cluster_IsHigh = !is.na(Cluster_Select),
                  All_IsHigh =  Sig_IsHigh & Path_IsHigh & Cluster_IsHigh)
} else{
  meta_df["Cluster_Select"] <- 1
}

# Bin spots into two groups by their cumulative status
if(is.null(sig_pval_threshold)){
  # Dummy signature empirical p-value threshold with a mean + 1SD of signature scores
  sig_pval_threshold <- sig_median_global + sig_sd_global * 2
  
  meta_df <- meta_df %>%
    dplyr::mutate(Sig_IsHigh = Sig_Emp_Pvalue >= sig_pval_threshold,
                  Path_IsHigh = !!rlang::sym(pathology_select) >= pathology_threshold,
                  Cluster_IsHigh = !is.na(Cluster_Select),
                  All_IsHigh =  Sig_IsHigh & Path_IsHigh & Cluster_IsHigh)
} else{
  meta_df <- meta_df %>%
    dplyr::mutate(Sig_IsHigh = Sig_Emp_Pvalue <= sig_pval_threshold,
                  Path_IsHigh = !!rlang::sym(pathology_select) >= pathology_threshold,
                  Cluster_IsHigh = !is.na(Cluster_Select),
                  All_IsHigh =  Sig_IsHigh & Path_IsHigh & Cluster_IsHigh)
}

print("Number of spots that have high (TRUE) signature vs low (FALSE):")
table(meta_df$Sig_IsHigh)
print("Number of spots that have high (TRUE) pathology vs low (FALSE):")
table(meta_df$Path_IsHigh)
print("Number of spots that have high (TRUE) cluster vs low (FALSE):")
table(meta_df$Cluster_IsHigh)

print("Number of spots that have high (TRUE) signature vs low (FALSE) and pass the pathology threshold:")
table(meta_df$All_IsHigh)

# Add updated meta data back to the Seurat object
data_seurat@meta.data <- meta_df %>%
  tibble::column_to_rownames("Coordinate")


## PLOT SIGNATURE DISTRIBUTIONS ----


# # Plot signature score split by sample with global thresholds
# filename <- paste0(output_path, "/hist_sig_scores_by_sample_", sample_name, "_", sig_names)
# create_hist_plot_my(meta_df, x_label=sig_names, fill_label="All_IsHigh",
#                     facet_var=c("user.Sample_Name", "fixed"),
#                     intercept=c(sig_median_global, sig_threshold), binwidth=0.05, add_density=FALSE,
#                     filename=filename, labels=c(sig_names, "Spot Number", "Signature score"))
# 
# # Plot signature empirical p-values split by sample with global thresholds
# filename <- paste0(output_path, "/hist_sig_pval_by_sample_", sample_name, "_", sig_names)
# create_hist_plot_my(meta_df, x_label="Sig_Emp_Pvalue", fill_label="All_IsHigh",
#                     facet_var=c("user.Sample_Name", "fixed"),
#                     intercept=c(0.01), binwidth=0.01, add_density=FALSE,
#                     filename=filename, labels=c(sig_names, "Spot Number", "Empirical P-value"), log_scale=TRUE)
# 
# # Plot pathology feature percent split by sample with global thresholds
# filename <- paste0(output_path, "/hist_path_by_sample_", sample_name, "_", pathology_select)
# create_hist_plot_my(meta_df, x_label=pathology_select, fill_label="All_IsHigh",
#                     facet_var=c("user.Sample_Name", "fixed"),
#                     intercept=c(75), binwidth=5, add_density=FALSE,
#                     filename=filename, labels=c(sig_names, "Spot Number", "Pathology class, %"))


# Plot spots selected for DEA using all critera
p <- spatial_dim_plot_my(data_seurat, group.by="All_IsHigh", title=sample_name)

# Write combo plot to file
filename <- paste0(output_path, "/spatial_status_final_", sample_name, "_", sig_names)
write_plot2file_my(p, filename, num_row=1, num_col=length(names(data_seurat@images))*1.5)

# Plot spot signature pvalues
p <- spatial_feature_plot_my(data_seurat, "Sig_Emp_Pvalue", min.cutoff="q0", max.cutoff="q100", name=sample_name)

# Write plot to file
filename <- paste0(output_path, "/spatial_value_sig_", sample_name, "_", sig_names)
write_plot2file_my(p, filename, num_row=1, num_col=length(names(data_seurat@images)))

# Plot pathology percent
p <- spatial_feature_plot_my(data_seurat, pathology_select, min.cutoff="q0", max.cutoff="q100", name=sample_name)

# Write plot to file
filename <- paste0(output_path, "/spatial_value_path_", sample_name, "_", pathology_select)
write_plot2file_my(p, filename, num_row=1, num_col=length(names(data_seurat@images)))


## FIND MARKERS ----


# Define DEA parameters
params <- cluster_params_my()

# Add latent parameters
params[["latent_vars"]] <- c("user.Sample_Name")

data_seurat@meta.data$All_IsHigh <- factor(data_seurat@meta.data$All_IsHigh, 
                                               levels=unique(data_seurat@meta.data$All_IsHigh))

# Put the clustering of interest into seurat clusters
data_seurat <- Seurat::SetIdent(data_seurat, value="All_IsHigh")

# Find markers of each cluster against the rest
# Only two groups defined so find_all=TRUE still works
markers_df <- seurat_find_markers_my(data_seurat, assay=params[["assay"]], 
                                     find_all=TRUE, group.by="All_IsHigh", 
                                     min.pct=params[["pct_min"]], 
                                     logfc.threshold=params[["logfc_threshold"]],
                                     test.use=params[["test_use"]], 
                                     latent.vars=params[["latent_vars"]])

# Find significant markers
filename_prefix <- paste0("dea_", sample_name, "_", sig_names)
markers_filt_df <- marker_analysis_my(markers_df, params, "All_IsHigh", 
                                      paste0(sample_name, "_", sig_names), output_path, filename_prefix) %>%
  dplyr::left_join(sig_invert_df, by="Symbol") %>%
  dplyr::filter(!grepl("^MT-|^RP[SL]", Symbol)) %>%
  tidyr::replace_na(list(InSignature = FALSE))

readr::write_delim(markers_filt_df, paste0(output_path, filename_prefix, "_markers_significant_ann.txt"), delim="\t")

# Select top biomarkers for the group
markers_top_df <- markers_filt_df %>%
  dplyr::filter(cluster == TRUE & 
                  pct_1 >= 75 & 
                  direction == "UP")

print(markers_top_df$Symbol)

readr::write_delim(markers_top_df, paste0(output_path, filename_prefix, "_mostrelevant.txt"), delim="\t")


## VISUALIZATIONS ----


# Dot plot of the most relevant genes
data_seurat@meta.data[["Group"]] <- sapply(1:nrow(meta_df), function(x) paste0(data_seurat@meta.data[x, c("All_IsHigh", "user.Sample_Name")], collapse=" "))

# Identify the width of the plot
num_col <- round(log2(length(markers_top_df$Symbol)) / 2) + 2
if(num_col < 1){
  num_col <- 1
}

# Go through every sample and create a dotplot of most significant genes
for(s in unique(data_seurat@meta.data$user.Group_Name)){
  # Identify barcodes of interest
  cells_select <- rownames(data_seurat@meta.data[which(data_seurat@meta.data$user.Group_Name == s),])
  data_seurat_loc <- subset(data_seurat, cells=cells_select)
  
  # Create a dot plot for a sample
  p <- Seurat::DotPlot(data_seurat_loc, features = markers_top_df$Symbol, group.by = "Group",
                       assay = assay, col.min = 0, col.max = 0.5, scale=TRUE) +
    Seurat::RotatedAxis()
  
  filename <- paste0(output_path, "/dea_dot_", s, "_", sig_names)
  write_plot2file_my(p, filename, num_row=1, num_col=num_col)
}

# Create a dot plot for the integrated cohort
p1 <- Seurat::DotPlot(data_seurat, features = markers_top_df$Symbol, group.by = "All_IsHigh",
                     assay = assay, col.min = 0, col.max = 0.5, scale=TRUE) +
  Seurat::RotatedAxis()

filename <- paste0(output_path, "/dea_dot_", sig_names)
write_plot2file_my(p1, filename, num_row=1, num_col=num_col)

# Format the filtered data for plotting
markers_filt_plot_df <- markers_filt_df %>%
  dplyr::mutate(avg_logFC = log10(2^avg_logFC),
                random = round(runif(nrow(markers_filt_df), 1, 60)),
                p_val_adj_neg_log10 = ifelse(p_val_adj_neg_log10 == Inf, 340 + random, p_val_adj_neg_log10),
                labels = ifelse(Symbol %in% markers_top_df$Symbol, Symbol, NA),
                size = round(pct_1/20, 1),
                TopHit = Symbol %in% markers_top_df$Symbol) %>%
  dplyr::filter(cluster == TRUE & direction == "UP" & pct_1 >= 60 & avg_logFC >= 1)

# Create a scatter plot illustrating significant gene distribution
p2 <- create_scatter_plot_my(markers_filt_plot_df, 
                            x_label="avg_logFC", y_label="p_val_adj_neg_log10", fill_label="TopHit",
                            shape=19, size="size", filename=filename, dot_labels="labels", do_fit=NULL,
                            labels=c("Fold change, log10", "Negative adjusted p-value, log10", "Tumor v Non-tumor in PDAC")) + 
  ggplot2::geom_hline(yintercept=320, color="black", linetype="dashed", size=0.2)

filename <- paste0(output_path, "/dea_scatter_pval_", sig_names)
write_plot2file_my(p2, filename, num_row=2, num_col=2)


# Format the filtered data for plotting
markers_filt_plot_df <- markers_filt_df %>%
  dplyr::mutate(labels = ifelse(pct_1 >= 99 | avg_logFC >= 180, Symbol, NA),
                TopHit = Symbol %in% markers_top_df$Symbol) %>%
  dplyr::filter(cluster == TRUE & direction == "UP" & p_val_adj_neg_log10 >= 300 & avg_logFC >= 2 & pct_1 >= 75)

# Create a scatter plot illustrating significant gene distribution
p3 <- create_scatter_plot_my(markers_filt_plot_df, 
                            x_label="avg_logFC", y_label="pct_1", fill_label="InSignature",
                            shape=21, size=6, filename=filename, dot_labels="labels", do_fit=NULL,
                            labels=c("Fold change, log2", "Percent spots expressing the gene in tumor", "Tumor v Non-tumor in PDAC")) + 
  ggplot2::geom_hline(yintercept=90, color="black", linetype="dashed", size=0.2)

filename <- paste0(output_path, "/dea_scatter_pct1_", sig_names)
write_plot2file_my(p3, filename, num_row=2, num_col=2)
