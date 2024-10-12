# Author: Anna Lyubetskaya. Date: 20-04-22

# Measure distance between two signature-defined features


## SETUP ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_signatures.R")
source("code/utils/utils_ggplot.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_heatmap.R")

source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_image.R")
source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_signatures.R")


## PARAMETERS ----


# Run name
run_name <- "FAP_Mac"

# The cohort of interest regex ID
# Seurat RDS files are tagged as follows
cohort_name <- "PDAC84_path12_merge"

# Gene abundance filters
sct_threshold <- 0.5
spot_threshold <- 5

# Signatures to select
sig1 <- "PDAC.P19.Macrophage"
sig2 <- "BMS.Program.FAP"

# Signature list
sig_select <- c(sig1, sig2)


## PATHS ----


# Input folder
input_path <- "XXXX"

# Output information
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_name, "_", run_name, "/")

# Create output folder
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)

# Path to a signature file
sig_path <- "data/import/Signatures/signatures_forDatta.txt"


## INGEST DATA ----


# Ingest a merged RDS Seurat objects
data_seurat <- readRDS(input_path)


## WRANGLE SIGNATURES ----


## Calculate signature scores and add them to Seurat meta data ----


# Find genes abundant in this sample
gene_list <- seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold,
                                             assay="SCT", slot="data", split_by="user.Sample_Name")

# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select)

# Signature names in the object
sig_names <- names(signature_list)

# Add signature scores to a seurat object
data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.", assay="SCT")

# Seurat renames column names, this step finds new signature names
sig_names_upd <- colnames(data_seurat@meta.data)[grep("sig.", colnames(data_seurat@meta.data))]


# Threshold signature distributions
for(s in sig_names_upd){
  # Set any negative signature score value to zero
  data_seurat@meta.data[data_seurat@meta.data[[s]] < 0, s] <- 0
  
  # Get rid of the strongest positive outlier signature scores
  data_seurat@meta.data[data_seurat@meta.data[[s]] > quantile(data_seurat@meta.data[[s]], 0.95), s] <- quantile(data_seurat@meta.data[[s]], 0.95)
}


## Calculate signature stats ----


# Wide tibble of signature scores
sig_wide_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate") %>%
  dplyr::select(dplyr::all_of(c("Coordinate", sig_names_upd)))

# Create a long tibble of signature scores
sig_df <- sig_wide_df %>%
  df_wide2long_my(key="Signature_name", val="Score") %>%
  dplyr::mutate(Sample_Name = gsub(":.+$", "", Coordinate))

# Calculate global sample stats
sample_stat_df <- sig_df %>%
  dplyr::group_by(Sample_Name) %>%
  dplyr::summarise(spot_num = dplyr::n_distinct(Coordinate)) %>%
  dplyr::ungroup()

# Calculate global signature stats
sig_stat_df <- sig_df %>%
  dplyr::group_by(Signature_name) %>%
  dplyr::summarise(score_median = median(Score),
                   score_sd = sd(Score)) %>%
  dplyr::mutate_if(is.numeric, round, 3) %>%
  dplyr::ungroup()

# Classify each spot
sig_status_df <- sig_df %>%
  dplyr::inner_join(sig_stat_df, by="Signature_name") %>%
  dplyr::mutate(Status = ifelse(Score >= score_median + score_sd, "T", "F"))

# Create a wide tibble with each signature score in each spot
sig_score_wide_df <- sig_status_df %>%
  df_long2wide_my(rows="Coordinate", cols="Signature_name", value="Score") %>%
  dplyr::mutate(Sample_Name = gsub(":.+$", "", Coordinate))

# Create a wide tibble with each signature status in each spot
sig_status_wide_df <- sig_status_df %>%
  df_long2wide_my(rows="Coordinate", cols="Signature_name", value="Status") %>%
  dplyr::mutate(Sample_Name = gsub(":.+$", "", Coordinate))

# Count the number of spots selected by both signatures and their intersection
spot_list1 <- sig_status_df %>%
  dplyr::filter(Status == "T" & Signature_name == sig_names_upd[1]) %>%
  dplyr::pull(Coordinate)

spot_list2 <- sig_status_df %>%
  dplyr::filter(Status == "T" & Signature_name == sig_names_upd[2]) %>%
  dplyr::pull(Coordinate)

spot_list3 <- intersect(spot_list1, spot_list2)

print(c(length(spot_list1), length(spot_list2), length(spot_list3)))


## CALCULATE DISTANCE BETWEEN SPOTS ----


spot_count_list <- list()
spot_dist_list <- list()

for(sample in unique(data_seurat@meta.data$user.Sample_Name)){
  
  # Select spots that are positive one or the other signature in this sample
  spots_select1 <- sig_status_df %>%
    dplyr::filter(Status == "T" & Sample_Name == sample & Signature_name == sig_names_upd[1]) %>%
    dplyr::pull(Coordinate) %>%
    unique
  
  # Select spots that are positive one or the other signature in this sample
  spots_select2 <- sig_status_df %>%
    dplyr::filter(Status == "T" & Sample_Name == sample & Signature_name == sig_names_upd[2]) %>%
    dplyr::pull(Coordinate) %>%
    unique
  
  # Number of spots corresponding to each of the signatures in the given sample
  spot_count_list[[sample]] <- c(sample, length(spots_select1), length(spots_select2))
  
  # Find spot coordinates for the 1st signature
  image1_df <- tibble::as_tibble(data_seurat@images[[sample]]@coordinates[c("imagerow", "imagecol")], rownames="Coordinate") %>%
    dplyr::filter(Coordinate %in% spots_select1) %>%
    dplyr::rename(X = imagecol, Y = imagerow)
  
  # Find spot coordinates for the 2nd signature
  image2_df <- tibble::as_tibble(data_seurat@images[[sample]]@coordinates[c("imagerow", "imagecol")], rownames="ClosestCoordinate") %>%
    dplyr::filter(ClosestCoordinate %in% spots_select2) %>%
    dplyr::rename(X = imagecol, Y = imagerow)
  
  # Find the closest barcode from one tibble in another
  spot_dist_list[[sample]] <- nearest_neighbor_df_my(image1_df, image2_df, value="ClosestCoordinate")
}

# Bind together spot number data
spot_count_df <- t(dplyr::bind_rows(spot_count_list)) %>%
  tibble::as_tibble() %>%
  dplyr::rename(Sample_Name = V1, Sig1_Spots = V2, Sig2_Spots = V3) %>%
  dplyr::inner_join(sample_stat_df, by="Sample_Name") %>%
  dplyr::mutate(Sig1_Perc = round(as.numeric(Sig1_Spots) / spot_num * 100, 1),
                Sig2_Perc = round(as.numeric(Sig2_Spots) / spot_num * 100, 1))

# Bind together spot distance data
spot_dist_df <- dplyr::bind_rows(spot_dist_list) %>%
  dplyr::mutate(Sample_Name = gsub(":.+$", "", Coordinate),
                Distance_log2 = round(log2(Distance + 1), 2))

# Summarize distance by sample
spot_dist_sum_df <- spot_dist_df %>%
  dplyr::group_by(Sample_Name) %>%
  dplyr::summarise(DistMean = mean(Distance_log2),
                   DistSD = sd(Distance_log2))


## VISUALIZE DATA ----


# Signature score distribution
filename <- paste0(output_path, "hist_sig_scores")
create_hist_plot_my(sig_status_df, x_label="Score", fill_label="Status", facet_var=c("Signature_name", "fixed"),
                    intercept=c(0, 1), binwidth=0.05, filename=filename, add_density=FALSE, log_scale=FALSE,
                    labels=c("Signature Score", "# Spots", sig_names_upd[1]))

# Percent spots corresponding to the pair of signatures in each sample
filename <- paste0(output_path, "bar_sig_spot_percent")
p <- create_bar_plot_my(spot_count_df %>%
                          dplyr::select(Sample_Name, Sig1_Perc, Sig2_Perc) %>%
                          df_wide2long_my(key="Signature_Name", val="SpotPercent") %>%
                          dplyr::mutate(Signature_Name = gsub("Sig1_Perc", sig1, Signature_Name),
                                        Signature_Name = gsub("Sig2_Perc", sig2, Signature_Name)),
                        x_label="Sample_Name", y_label="SpotPercent", fill_label="Signature_Name", 
                        position="stack", facet_var=c("Signature_Name", "fixed"), 
                        filename=NULL, labels=c("Sample_Name", "% spots", ""), reorder_x=TRUE)
write_plot2file_my(p, filename, num_row=1, num_col=4, width=NULL, height=NULL)

# Distances between two sets of signature positive spots in each sample
filename <- paste0(output_path, "bar_sig_spot_distance")
p <- create_bar_plot_my(spot_dist_sum_df, x_label="Sample_Name", y_label="DistMean", fill_label="Sample_Name", 
                        filename=NULL, labels=c("Sample Name", "Distance, log2, mean", paste(sig1, sig2)), reorder_x=TRUE,
                        error_label="DistSD")
write_plot2file_my(p, filename, num_row=1, num_col=4, width=NULL, height=NULL)


# Spot status
data_seurat@meta.data[which(rownames(data_seurat@meta.data) %in% spot_dist_df$Coordinate), "Status"] <- sig1
data_seurat@meta.data[which(rownames(data_seurat@meta.data) %in% spot_dist_df$ClosestCoordinate), "Status"] <- "Proximal"
data_seurat@meta.data[which(rownames(data_seurat@meta.data) %in% setdiff(spot_list2, spot_dist_df$ClosestCoordinate)), "Status"] <- sig2
data_seurat@meta.data[which(is.na(data_seurat@meta.data$Status)), "Status"] <- "Exclude"


# Spatial plot of selected spots
cols <- c("darkblue", "red", "darkgoldenrod1", "lightgrey")
names(cols) <- c(sig1, "Proximal", sig2, "Exclude")
filename <- paste0(output_path, "spatial_spots")
p <- batch_spatial_feature_plot_my(list(Data = data_seurat), c("Status"), output_file=filename, title=run_name,
                                   plot_type=rep("d", 1), col_list=list(Status=cols))


## SAVE DATA ----


# Signature scores, calls, and distances joined
data_summary_df <- sig_status_wide_df %>%
  dplyr::inner_join(sig_score_wide_df, by=c("Coordinate", "Sample_Name")) %>%
  dplyr::left_join(spot_dist_df, by=c("Coordinate", "Sample_Name"))

# Write signature and distance data to a file
filename <- paste0(output_path, "table_", cohort_name, ".txt")
readr::write_delim(data_summary_df, filename, delim="\t")


# Add distance information to the original object
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate") %>%
  dplyr::inner_join(data_summary_df %>%
                      dplyr::select("Coordinate", "Distance"), by="Coordinate")

# Update meta.data information
data_seurat@meta.data <- meta_df %>%
  tibble::column_to_rownames("Coordinate")

# # Subset and save object
# data_seurat <- subset(data_seurat, cells=spot_dist_df$Coordinate)
# filename <- gsub(".rds", "_dist.rds", input_path)
# saveRDS(data_seurat, filename)
