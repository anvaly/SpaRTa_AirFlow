# Author: Anna Lyubetskaya. Date: 20-08-06


##_ SETUP ENVIRONMENT _##


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_in_out.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_ggplot.R")
#source("code/utils/utils_signatures.R")
#source("code/utils/utils_specialized_plots.R")

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_vis.R")

limit_val_mb <- 4000
options(future.globals.maxSize = limit_val_mb*1024^2)


## PARAMETERS ----


# The cohort of interest regex ID
cohort_name <- "Syngeneic_FF"
# Number of features to use in integration
feature_num <- 3000


## PATHS ----


# Input folder
input_folder <- "C:/Users/lyubetsa/Documents/Projects_NewTech/10X_ST_Reports/"

# Seurat RDS files are tagged as follows
data_regex_files <- "_all.Sobj.rds"

# Output folder
output_folder <- "C:/Users/lyubetsa/Documents/Projects_NewTech/10X_ST_Integration/"

# DEA thresholds
fc_threshold <- 1
fdr_threshold <- 2  # in -log10(FDR) space


## INGEST DATA ----


# Find all files following the regex
file_list <- dir(input_folder, pattern=paste0(cohort_name, ".*", data_regex_files), full.names=TRUE)

# Ingest clean Seurat data
seurat_list <- list()  
for(rds_file in file_list){
  seurat_list[gsub(paste0("^.+/|", data_regex_files), "", rds_file)] <- read_seurat_rds_my(rds_file, do_subset=FALSE, do_cluster=FALSE, marker_file=NULL)
}


## UMAP UNNORMALIZED DATA ----


combined_seurat <- seurat_list[[1]]
for(s in seurat_list[2:length(seurat_list)]){
  combined_seurat <- merge(combined_seurat, s, project = cohort_name)
}

# Perform SCTransform normalization
combined_seurat <- Seurat::SCTransform(combined_seurat, return.only.var.genes = FALSE)

combined_seurat <- cluster_dea_my(combined_seurat, features=NULL, num_dimensions=30, resolution=0.5)

# Plot a classic two-dimensional UMAP projection of integrated data
for(var in c("meta.SampleName", "seurat_clusters")){
  p1 <- dim_plot_my(combined_seurat, var)
  filename <- paste0(output_folder, "/dim_plot_init_", cohort_name, "_", feature_num, "_", var)
  write_plot2file_my(p1, filename, num_row=2, num_col=2)
}


## JOINTLY NORMALIZE DATA ----


# Select features for integration
feature_list <- Seurat::SelectIntegrationFeatures(object.list = seurat_list, nfeatures = feature_num)
# Calculate Pearson residuals
seurat_list <- Seurat::PrepSCTIntegration(object.list = seurat_list, anchor.features = feature_list, verbose = FALSE)


# Identify anchors
anchors_list <- Seurat::FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = feature_list, verbose = FALSE)
integrated_seurat <- Seurat::IntegrateData(anchorset = anchors_list, normalization.method = "SCT", verbose = FALSE)


## CLUSTER DATA ----


integrated_seurat <- cluster_dea_my(integrated_seurat, features=NULL, num_dimensions=30, resolution=0.5)

integrated_labels <- integrated_seurat@meta.data %>% 
  dplyr::select(seurat_clusters, meta.SampleID) %>% 
  tibble::rownames_to_column()
  
count_clusters <- integrated_labels %>%
  dplyr::group_by(meta.SampleID, seurat_clusters) %>% 
  dplyr::tally() %>%
  df_long2wide_my(rows="meta.SampleID", cols="seurat_clusters", value="n")


# Integrated cohort data to RDS
filename <- paste0(output_folder, "/integrated_", cohort_name, "_", feature_num, ".rds")
saveRDS(integrated_seurat, file = filename)


## CLUSTER DATA ----


# Plot a classic two-dimensional UMAP projection of integrated data
for(var in c("seurat_clusters", "meta.Protocol", "meta.SampleName")){
  p1 <- dim_plot_my(integrated_seurat, var)
  filename <- paste0(output_folder, "/dim_plot_", cohort_name, "_", feature_num, "_", var)
  write_plot2file_my(p1, filename, num_row=2, num_col=2)
}


# Plot clusters on the original slice
for(sample in names(seurat_list)){
  integrated_clusters <- integrated_labels %>%
    dplyr::filter(meta.SampleID == sample) %>%
    dplyr::mutate(rowname = gsub("_\\d+$", "", rowname))
  
  index <- match(integrated_clusters$rowname, rownames(seurat_list[[sample]]@meta.data))
  seurat_list[[sample]]@meta.data[index, "integrated_clusters"] <- integrated_clusters$seurat_clusters
  
  # Plot clusters on the original slice
  p2 <- Seurat::SpatialDimPlot(seurat_list[[sample]], group.by="integrated_clusters")
  
  filename <- paste0(output_folder, "/sp_dim_plot_", cohort_name, "_", feature_num, "_", sample)
  write_plot2file_my(p2, filename)
}


## FIND MARKERS ----


markers_df <- compare_clusters_my(integrated_seurat, TRUE)

# DEA tibble
filename <- paste0(output_folder, "/dea_", cohort_name, "_", feature_num, ".csv")
readr::write_delim(markers_df, path=filename, delim="\t", col_names=TRUE)

# Select rows for significant genes
dea_table <- markers_df %>%
  dplyr::filter(abs(avg_logFC) >= 2 & 
                  p_val_adj_neg_log10 >= fdr_threshold &
                  (pct_1 >= 75 | pct_2 >= 75)) %>%
  dplyr::mutate(p_val_adj_neg_log10 = ifelse(p_val_adj_neg_log10 == Inf, 1000, p_val_adj_neg_log10)) %>%
  dplyr::arrange(desc(fdr_threshold), desc(abs(avg_logFC)))


## VISUALIZE MARKERS ----


markers_list <- unique(dea_table$Symbol)

filename <- paste0(output_folder, "/markers_hm_", cohort_name, "_", feature_num)
seurat_heatmap_my(integrated_seurat, markers_list, plot_title=paste("Integrated", cohort_name), filename=filename)
