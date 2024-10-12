# Author: Anna Lyubetskaya. Date: 20-06-03
# Using signatures, partition data into groups and perform DEA analysis on the groups


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_signatures.R")
source("code/utils/utils_rna_diff_expr.R")
source("code/utils/utils_ggplot.R")
source("code/utils/utils_tibble.R")

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_signatures.R")
source("code/R/Utils/utils_10X_dea.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_10X_matrix.R")
#source("code/R/Utils/utils_10X_vis.R")

source("code/R/e_DEA/dea_signature_groups_utils.R")


## PARAMETERS ----


# Cohort name
cohort_name <- "Syng_late"
# Subset cohort to this specific set
cohort_subset <- "Syng"

# Seurat RDS files are tagged as follows
# Syng_late_mt-ribo-mouse_merge_annotated.rds
data_regex_files <- "_mt-ribo-mouse_merge_annotated.rds"

# Two signature of interest, one is the focus, the other is a condition
sig_names_list <- list("CellLine_B16F10" = c(0, 1),
                       "CellLine_MC38" = c(0, 0.5),
                       "Immune_Macrophage" = c(0, 0.5))
# Category names
category_names_list <- paste0("Category_", names(sig_names_list))

# Contrasts
contrast_list <- list()

# List of misc fields to keep in meta data
misc_field_list <- c("Coordinate", "Batch", "user.Tissue")
group_by_var <- "user.Tissue"

## Parameters for DEA ----

# Minimum % in spots
pct_min <- 0.1
# Minimum FC difference to test
logfc_threshold <- 0.15

# FC signficance threshold
sign_fc_threshold <- 0.25
# Signficance p-value threshold
sign_pval_adj_neglog10_threshold <- 2


## PATHS ----


# Location of annotated Seurat objects
input_path <- paste0("XXXX", cohort_name, data_regex_files)

# Path to a signature file
sig_path <- "data/import/Signatures/signatures_syngeneics_mixed_dea_Jan21.txt"

# Output information
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_name, "/")

# Create output folder
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, sig_names=names(sig_names_list), sig_length_min=1, sig_length_max=1000, ratio_threshold=0)

# Find a signatures list for every gene
sig_inverted_list <- invert_list_my(signature_list)

# Ingest an RDS object
data_seurat <- readRDS(input_path)


## WRANGLE DATA ----


# Add signature scores to the integrated dataset
data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="")

# Meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate") %>%
  dplyr::select(dplyr::all_of(misc_field_list))

# Wide tibble of signature scores
sig_wide_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate") %>%
  dplyr::select(dplyr::all_of(c("Coordinate", names(sig_names_list))))

# Make a long tibble of signature scores and add meta data
sig_df <- sig_wide_df %>%
  df_wide2long_my("Signature", val="Score") %>%
  dplyr::inner_join(meta_df, by="Coordinate")

# Plot signature distribution for reference
filename <- paste0(output_path, "hist_sig_scores")
create_hist_plot_my(sig_df, x_label="Score", fill_label="Signature", facet_var=c("Batch", "fixed"), 
                    intercept=c(0, 0.5, 1), binwidth=0.1, filename=filename, labels=NULL)

# Add meta data to the wide tibble
sig_wide_df <- sig_wide_df %>%
  dplyr::inner_join(meta_df, by="Coordinate") %>% 
  dplyr::filter(grepl(cohort_subset, Batch))

# Partition spots into hi-lo quadrants based on scores of each signature
for(sig in names(sig_names_list)){
  sig_wide_df <- sig_wide_df %>%
    vector_categories_hi_lo_my(sig, dn_threshold=sig_names_list[[sig]][1], up_threshold=sig_names_list[[sig]][2])
  
  colnames(sig_wide_df) <- gsub("^Category$", paste0("Category_", sig), colnames(sig_wide_df))
}

# Add a joint category name
sig_wide_df[["Category"]] <- sapply(1:nrow(sig_wide_df), function(x) paste0(sig_wide_df[x, category_names_list], collapse="-"))
sig_wide_df[["Category"]] <- paste0(sig_wide_df[[group_by_var]], "_", sig_wide_df[["Category"]])

# Count number of categorical labels
bin_stat_wide_df <- sig_wide_df %>% 
  dplyr::group_by(Category) %>% 
  dplyr::summarise(n = dplyr::n_distinct(Coordinate))

filename <- paste0(output_path, "group_sizes.txt")
readr::write_delim(bin_stat_wide_df, filename, delim="\t")


## DEA: FindMarkers ----


# Check that data can be merged
table(sig_wide_df$Coordinate == rownames(data_seurat@meta.data))
# Add categories to the Seurat object
data_seurat@meta.data[["Category"]] <- sig_wide_df$Category

# Put the clustering of interest into seurat clusters
data_seurat <- Seurat::SetIdent(data_seurat, value="Category")

cluster_pairs <- list(c("Syng_B16_Hi-Lo-Med", "Syng_B16_Hi-Lo-Lo"),  # Immune infiltration in pure B16
                      c("Syng_B16_Hi-Med-Lo", "Syng_B16_Hi-Lo-Lo"),  # "MC38" infiltration in pure B16
                      
                      c("Syng_MC38_Lo-Hi-Med", "Syng_MC38_Lo-Hi-Lo"),  # Immune infiltration in pure B16: Med / Lo
                      c("Syng_MC38_Lo-Hi-Hi", "Syng_MC38_Lo-Hi-Med"),  # Immune infiltration in pure B16: Hi / Med
                      c("Syng_MC38_Lo-Hi-Hi", "Syng_MC38_Lo-Hi-Lo"),  # Immune infiltration in pure B16: Hi / Lo
                      
                      c("Syng_B16-MC38_Hi-Lo-Lo", "Syng_B16_Hi-Lo-Lo"),  # Mixed B16 v pure B16
                      c("Syng_B16-MC38_Lo-Hi-Med", "Syng_MC38_Lo-Hi-Med"),  # Mixed MC38 v pure MC38
                      
                      c("Syng_B16-MC38_Hi-Lo-Med", "Syng_B16-MC38_Hi-Lo-Lo"),  # Immune infiltration in mixed B16
                      c("Syng_B16-MC38_Lo-Hi-Hi", "Syng_B16-MC38_Lo-Hi-Med"),  # Immune infiltration in mixed MC38
                      
                      c("Syng_B16-MC38_Hi-Lo-Med", "Syng_B16-MC38_Hi-Hi-Med"),  # MC38 infiltration in mixed B16
                      c("Syng_B16-MC38_Lo-Hi-Med", "Syng_B16-MC38_Hi-Hi-Med"))  # B16 infiltration in mixed MC38

# Perform basic FindMarkers on select categories
markers_df <- seurat_find_markers_my(data_seurat, find_all=FALSE, group.by="Category", cluster_pairs=do.call("cbind", cluster_pairs),
                                     test.use="wilcox", min.pct=pct_min, logfc.threshold=logfc_threshold)

# Write all marker data to a file
filename <- paste0(output_path, "/findmarkers_simple")
readr::write_delim(markers_df, delim="\t", paste0(filename, ".txt"))


# Filter marker tibble to significant markers only
filename <- paste0(output_path, "/findmarkers_simple_filt")
markers_filt_df <- filter_and_vis_markers_my(markers_df, fc_threshold=sign_fc_threshold, 
                                             pval_adj_neglog10_threshold=sign_pval_adj_neglog10_threshold, 
                                             output_file=filename)

# Gather stats for DE markers for each clustering
stat_df <- markers_filt_df %>%
  dplyr::group_by(cluster) %>% 
  dplyr::summarise(UP = sum(avg_logFC > 0),
                   DN = sum(avg_logFC < 0),
                   COUNT = dplyr::n_distinct(Symbol))

# Write marker stats to a file
filename <- paste0(output_path, "/findmarkers_simple_stats.txt")
readr::write_delim(stat_df, delim="\t", filename)
