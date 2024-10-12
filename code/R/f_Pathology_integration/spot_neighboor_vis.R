# Author: Anna Lyubetskaya. Date: 22-09-02
# Create various visualizations of a set of signatures as well as specific genes within them


## ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_tibble.R")
source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_10X_vis.R")


## PARAMETERS ----


# Name of the processed Seurat object
sample_name <- "PDAC108_path14_harmonystr_rpca_sct"

# Name for the signature group being plotted
run_name <- "str_nbs_path"

# User-defined clustering resolution
resolution <- "integrated_snn_res.0.2"

# Select specific clusters
cluster_select <- NULL  # c(3, 10, 0, 4, 6, 9)


## PATHS ----


# Location of pre-processed data
input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, run_name, "_", sample_name, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)

# Import additional meta data for the object
input_meta_path <- c("XXXX",
                     "XXXX")


## INGEST DATA ----


# Seurat data
data_seurat <- readRDS(input_path)

# Import additional meta data from a table
meta_list <- list()
for(f in input_meta_path){
  meta_list[[f]] <- readr::read_delim(f, delim="\t") %>%
    dplyr::select(NbDistance, Feature, Coordinate, ScoreMean)
}

# Combined additional meta data
meta_extra_df <- dplyr::bind_rows(meta_list) %>%
  dplyr::mutate(Feature = gsub("sig.|Pathology.|.percent|P19.|pdac.", "", Feature)) %>%
  tidyr::drop_na()

# Meta parameters select
meta_select <- unique(meta_extra_df$Feature)


## WRANGLE DATA ----


# Extract coordinate meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

if(!is.null(cluster_select)){
  meta_df <- meta_df %>%
    dplyr::filter(!!rlang::sym(resolution) %in% cluster_select)
}

# Crete a long tibble for meta data
meta_long_df <- meta_df %>%
  dplyr::select(dplyr::all_of(c("Coordinate", resolution))) %>%
  dplyr::inner_join(meta_extra_df, by="Coordinate") %>%
  dplyr::mutate(NbDistance = as.character(NbDistance)) %>%
  tidyr::drop_na()

# Identify common lists of coordinates and features
barcode_count <- table(meta_long_df$Coordinate)
barcode_list <- names(barcode_count[which(barcode_count == median(barcode_count))])
feature_count <- table(meta_long_df$Feature)
feature_list <- names(feature_count[which(feature_count == median(feature_count))])

# Find the smallest score
score_min <- min(meta_long_df$ScoreMean)

# Filter down long meta data tibble
meta_long_df <- meta_long_df %>%
  dplyr::filter(Coordinate %in% barcode_list & Feature %in% feature_list)

# Average features by cluster and feature
meta_sum_df <- meta_long_df %>%
  dplyr::rename(Score = ScoreMean) %>%
  dplyr::group_by(!!rlang::sym(resolution), NbDistance, Feature) %>%
  dplyr::summarise(ScoreMean = mean(Score),
                   ScoreSD = sd(Score))

# Create a set of wide spot tibbles by neighborhood
meta_wide_list <- list()
for(nb in sort(unique(meta_long_df$NbDistance))){
  meta_wide_list[[nb]] <- meta_long_df %>%
    dplyr::filter(NbDistance == nb) %>%
    dplyr::select(-NbDistance) %>%
    df_long2wide_my(rows="Coordinate", cols="Feature", value="ScoreMean") %>%
    tibble::column_to_rownames("Coordinate")
}


## VISUALIZE DATA ----


# Define colors
cols <- define_cols_my(n=length(sort(unique(meta_long_df$NbDistance))), col_type="viridis")
names(cols) <- sort(unique(meta_long_df$NbDistance))


# Cycle through features
for(feature in meta_select){
  
  p <- create_box_plot_my(meta_long_df %>%
                            dplyr::filter(Feature == feature), 
                          x_label="NbDistance", y_label="ScoreMean", fill_label="NbDistance",
                          facet_var=c(resolution, "fixed"), outlier_shape=19,
                          filename=NULL, labels=c(resolution, "Score", feature), cols=cols)
  
  filename <- paste0(output_path, "/box_", feature, "_", sample_name)
  write_plot2file_my(p, filename, num_row=2, num_col=0.5)
  
  
  p <- create_bar_plot_my(meta_sum_df %>%
                            dplyr::filter(Feature == feature), 
                          x_label="NbDistance", y_label="ScoreMean", fill_label="NbDistance",
                          facet_var=c(resolution, "fixed"), # error_label="ScoreSD",
                          filename=NULL, labels=c(resolution, "Score", feature), cols=cols)
  
  filename <- paste0(output_path, "/bar_", feature, "_", sample_name)
  write_plot2file_my(p, filename, num_row=2, num_col=0.5)
  
}


## CALCULATE DISTANCE FROM CENTER SPOT ----


# Difference in signal between neighborhoods
diff_wide_df <- (meta_wide_list[[2]] - meta_wide_list[[1]]) %>%
  tibble::rownames_to_column("Coordinate") %>%
  dplyr::inner_join(meta_df, by="Coordinate")


for(col in unique(meta_long_df$Feature)){
  
  mean_diff <- mean(diff_wide_df[[col]])
  sd_diff <- sd(diff_wide_df[[col]])
  threshold1 <- mean_diff + sd_diff * 2
  threshold2 <- mean_diff - sd_diff * 2
  
  diff_wide_df <- diff_wide_df %>%
    dplyr::mutate(Outlier = ifelse(!!rlang::sym(col) >= threshold1, "Up", NA),
                  Outlier = ifelse(!!rlang::sym(col) <= threshold2, "Dn", Outlier))
  
  filename <- paste0(output_path, "/hist_dist_", col, "_", sample_name)
  p <- create_hist_plot_my(diff_wide_df, x_label=col, fill_label="Outlier", binwidth=0.01, intercept=c(0, mean_diff, -0.5, 0.5),
                           filename=filename, labels=c(paste(col, "signature score difference"), "Spot number", col))
  
  table(diff_wide_df$Outlier)
  
  diff_loc1_df <- diff_wide_df %>%
    dplyr::filter(Outlier %in% c("Up"))
  
  diff_loc2_df <- diff_wide_df %>%
    dplyr::filter(Outlier %in% c("Dn"))
  
  a = table(diff_wide_df[[resolution]])
  b = table(diff_loc1_df[[resolution]])
  c = table(diff_loc2_df[[resolution]])
  
  round(b/a*100)
  round(c/a*100)
  
}
