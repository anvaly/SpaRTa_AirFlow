# Author: Anna Lyubetskaya. Date: 22-08-23
# Investigate various annotation factors relative to each other


## ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_signatures.R")


## PARAMETERS ----


# Name of the analysis to use in folder/file names
run_name <- "DC_FAP_stroma"

# Path to processed Seurat data
sample_name <- "PDAC84_path12_merge"
sample_exclude <- NULL


# ## STROMA

# Name for the signature group being plotted
sig_select <- c("BMS.Program.FAP", "PDAC.P19.Macrophage")
sig_rename <- c("FAP", "macrophage")
names(sig_rename) <- sig_select

# Name of the pathology field to filter by groups
pathology_select <- c("Pathology.Stroma.percent")
pathology_name <- "StromaTME"


# SCT and spot number thresholds for gene filtering
sct_threshold <- 0.5
spot_threshold <- 100

# Assay and slot to use to calculate signature score
assay <- "SCT"
slot <- "data"

# Feature # in a spot threshold
feature_threshold <- 2000

# Subset and write out a Seurat object
subset_def <- NULL


## PATHS ----


# Path to a signature file
sig_path <- "data/import/Signatures/signatures_forDatta.txt"

# Location of pre-processed data
input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, run_name, "_", sample_name, "/")

# Path to which to write a subset Seurat object
output_path_subset <- "XXXX"

# Log information
log_file <- paste0(output_path, "log.txt")

# Create output folders
dir.create(output_path_subset, showWarnings = FALSE)
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
  barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data$nFeature_Spatial >= feature_threshold),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}


## CALCULATE TARGET SIGNATURE SCORES ----


## Calculate signature scores and add them to Seurat meta data

# Find genes abundant in this sample
gene_list <- unique(seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold, 
                                                    assay=assay, slot=slot))

# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select, 
                                            sig_length_min=1, sig_length_max=1000, ratio_threshold=0)

# Rename signatures to shorter names
if(!is.null(sig_rename)){
  names(signature_list) <- unname(sig_rename[names(signature_list)])
}

# Find all genes in signatures and their signature assignations
sig_invert_df <- invert_list_my(signature_list)

# Add signature scores to a Seurat object
data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.", assay=assay)

# Find column names in the Seurat object - Seurat can rename original signatures
sig_names <- colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))])

# Find the properties of the genes in the target signature
sig_gene_stat_list <- list()
for(s in names(signature_list)){
  sig_gene_stat_list[[s]] <- seurat_gene_properties_my(data_seurat, signature_list[[s]], assay=assay, sct_threshold=sct_threshold) %>%
    dplyr::mutate(Signature_Name = s)
  
  # Write stats to file
  filename <- paste0(output_path, "/table_sig_by_gene_", sample_name, "_", s, ".txt")
  readr::write_delim(sig_gene_stat_list[[s]], filename, delim="\t")
}


## WRANGLE META DATA ----


# Extract coordinate meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

# Create a joint pathology label
meta_df[[pathology_name]] <- matrixStats::rowMaxs(as.matrix(meta_df[pathology_select]))

# Add user.Sample_Name to meta data if absent: happens in stand-alone samples
if(!"user.Sample_Name" %in% colnames(meta_df)){
  meta_df["user.Sample_Name"] = data_seurat@misc$user.Sample_Name
}

# Generalize the pathology percent
interval_factors <- sort(unique(round(meta_df[[pathology_name]]/10)*10))
meta_df <- meta_df %>%
  dplyr::mutate(PathologyIntervals = factor(round(!!rlang::sym(pathology_name)/10)*10, levels=interval_factors),
                Path1 = !!rlang::sym(pathology_name) >= 50,
                Path2 = !!rlang::sym(pathology_name) <= 25)

# Rename column with pathology calls
colnames(meta_df) <- gsub("^Path1$", paste0(pathology_name, "_True"), colnames(meta_df))
colnames(meta_df) <- gsub("^Path2$", paste0(pathology_name, "_False"), colnames(meta_df))


# Threshold signatures
for(s in sig_names){
  # Set any negative signature score value to zero
  meta_df[meta_df[[s]] < 0, s] <- 0
  # Get rid of the strongest positive outlier signature scores
  meta_df[meta_df[[s]] > quantile(meta_df[[s]], 0.95), s] <- quantile(meta_df[[s]], 0.95)
  
  # Global signature stats
  sig_mean_global <- round(mean(meta_df[[s]]), 3)
  sig_sd_global <- round(sd(meta_df[[s]]) , 3)
  sig_threshold_h <- sig_mean_global + sig_sd_global * 2
  sig_threshold_l <- sig_mean_global
  
  cat(s, sig_mean_global, sig_sd_global, sig_threshold_h, sig_threshold_l, "\n")
  
  # New column names
  name1 <- paste0(s, "_True")
  name2 <- paste0(s, "_False")
  
  # ID tumor categories
  meta_df <- meta_df %>%
    dplyr::mutate(Sig1 = !!rlang::sym(s) >= sig_threshold_h,
                  Sig2 = !!rlang::sym(s) <= sig_threshold_l,
                  PathSig1 = !!rlang::sym(pathology_name) >= 50 & !!rlang::sym(s) >= sig_threshold_h,
                  PathSig2 = !!rlang::sym(pathology_name) <= 25 & !!rlang::sym(s) <= sig_threshold_l)
  
  colnames(meta_df) <- gsub("^Sig2$", paste0(s, "_False"), colnames(meta_df))
  colnames(meta_df) <- gsub("^Sig1$", paste0(s, "_True"), colnames(meta_df))
  colnames(meta_df) <- gsub("^PathSig2$", paste0(s, "_", pathology_name, "_False"), colnames(meta_df))
  colnames(meta_df) <- gsub("^PathSig1$", paste0(s, "_", pathology_name, "_True"), colnames(meta_df))
}

# Find samples where pathology seems to disagree with signatures a lot
# Find samples where signatures seem to disagree with pathology a lot
meta_df <- meta_df %>%
  dplyr::mutate(PathologyFail = !!rlang::sym(paste0(pathology_name, "_False")) == TRUE  & (!!rlang::sym(paste0(sig_names[[1]], "_True")) == TRUE | !!rlang::sym(paste0(sig_names[[2]], "_True")) == TRUE),
                SignatureFail = !!rlang::sym(paste0(pathology_name, "_True")) == TRUE & (!!rlang::sym(paste0(sig_names[[1]], "_False")) == TRUE & !!rlang::sym(paste0(sig_names[[2]], "_False")) == TRUE))

write("Pathology FAIL", log_file, append="FALSE")
write(paste(names(table(meta_df$PathologyFail)), collapse="\t"), log_file, append="TRUE")
write(table(meta_df$PathologyFail), log_file, append="TRUE")

write("\nSignature FAIL", log_file, append="TRUE")
write(paste(names(table(meta_df$SignatureFail)), collapse="\t"), log_file, append="TRUE")
write(table(meta_df$SignatureFail), log_file, append="TRUE")


## ID OUTLIER SPOTS ----


# Count spots by sample
count_spots_df <- meta_df %>%
  dplyr::group_by(user.Sample_Name) %>%
  dplyr::summarise(CountAll = dplyr::n_distinct(Coordinate))


# Count group sizes for each category
stat_list <- list()

category_col_list <- c("PathologyFail", "SignatureFail", 
                       paste0(pathology_name, "_True"), paste0(pathology_name, "_False"),
                       paste0(sig_names, "_True"), paste0(sig_names, "_False"),
                       paste0(sig_names, "_", pathology_name, "_True"), paste0(sig_names, "_", pathology_name, "_False"))

for(g in category_col_list){
  stat_list[[g]] <- meta_df %>%
    dplyr::filter(!!rlang::sym(g) == TRUE) %>%
    dplyr::group_by(user.Sample_Name) %>%
    dplyr::summarise(Count = dplyr::n_distinct(Coordinate)) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(count_spots_df, by="user.Sample_Name") %>%
    dplyr::mutate(Percent = round(Count / CountAll * 100),
                  Group = g) %>%
    dplyr::filter(Count >= 20)
}

intersect(stat_list[[1]]$user.Sample_Name, stat_list[[2]]$user.Sample_Name)

write(paste0("\n", pathology_name, "_True", collapse=""), log_file, append="TRUE")
write(paste(names(table(meta_df[paste0(pathology_name, "_True")])), collapse="\t"), log_file, append="TRUE")
write(table(meta_df[paste0(pathology_name, "_True")]), log_file, append="TRUE")

write(paste0("\n", paste(sig_names, collapse=" "), "\t", pathology_name, "_True", collapse=""), log_file, append="TRUE")
write.table(table(meta_df[paste0(sig_names, "_", pathology_name, "_True")]), log_file, append="TRUE")


# Write meta data to file
filename <- paste0(output_path, "/table_meta_data_", sample_name, ".txt")
readr::write_delim(meta_df %>%
                     dplyr::select(dplyr::all_of(c("Coordinate", "user.Sample_Name", sig_names, pathology_select, pathology_name, category_col_list))), filename, delim="\t")

# Add meta data to the Seurat object
data_seurat@meta.data <- meta_df %>%
  tibble::column_to_rownames("Coordinate")


# Subset and write out a Seurat object
if(!is.null(subset_def)){
  table(meta_df[subset_def])
  
  # Find samples and corresponding barcodes to subset
  subset_df <- meta_df %>%
    dplyr::filter(!!rlang::sym(subset_def) == TRUE) %>%
    dplyr::select(Coordinate, user.Sample_Name)
  
  for(s in unique(subset_df$user.Sample_Name)){
    # Barcodes for the select samples
    cell_list <- subset_df %>%
      dplyr::filter(user.Sample_Name == s) %>%
      dplyr::pull(Coordinate)
    
    filename_rds <- paste0(output_path_subset, sample_name, "_", s, "_", subset_def, ".rds")
    
    if(length(cell_list) >= 100 && !file.exists(filename_rds)){
      data_subset_seurat <- subset(data_seurat, cells=cell_list)
      
      # Only leave the relevant image
      data_subset_seurat@images <- data_subset_seurat@images[unique(data_subset_seurat@meta.data$user.Sample_Name)]
      
      saveRDS(data_subset_seurat, filename_rds)
    }
  }
  
}


## PLOT CATEGORY SPOT COUNTS ----


# Spot counts at various thresholds
stat_df <- dplyr::bind_rows(stat_list)

# Wide tibble of spot counts
stat_wide_df <- stat_df %>%
  df_long2wide_my(rows="user.Sample_Name", cols="Group", value="Percent")

# Write meta data to file
filename <- paste0(output_path, "/table_spotcounts_", sample_name, ".txt")
readr::write_delim(stat_wide_df, filename, delim="\t")


# Barplot of percent spots in each selected category
for(c in category_col_list){
  stat_loc_df <- stat_df %>%
    dplyr::filter(Group == c)
  
  p <- create_bar_plot_my(stat_loc_df, x_label="user.Sample_Name", y_label="Percent", fill_label="Group",
                          facet_var=c("Group", "fixed"), filename=NULL, reorder_x=TRUE,
                          labels=c("Sample Names", "Percent of spots in the category", c))
  
  filename <- paste0(output_path, "/bar_spotnum_", c, "_", sample_name)
  write_plot2file_my(p, filename, num_row=2, num_col=4)
  
}


## PLOT FACTOR DISTRIBUTIONS ----


# Define categorical color palette for pathology labels
cols <- define_cols_my(n=length(interval_factors))

# Plot signature scores by pathology class category
for(s in sig_names){
  filename <- paste0(output_path, "/box_sig_v_path_", sample_name, "_", s, "_", pathology_name)
  create_box_plot_my(meta_df, x_label="PathologyIntervals", y_label=s, fill_label="PathologyIntervals",
                     filename=filename, labels=c(pathology_name, s, ""), cols=cols)
}

# Compare two signatures by pathology category
filename <- paste0(output_path, "/scatter_sig_v_path_", sample_name, "_", run_name, "_", pathology_name)
create_scatter_plot_my(meta_df, x_label=sig_names[[1]], y_label=sig_names[[2]], fill_label="PathologyIntervals",
                       facet_var=c("PathologyIntervals", "fixed"),
                       filename=filename, labels=c(sig_names[[1]], sig_names[[2]], ""), 
                       cols=define_cols_my(n=length(interval_factors)), stroke=0)


## PLOT OUTLIER SPOTS ----


# Plot spots that have strongly discordant pathology and signature
for(g in c("PathologyFail", "SignatureFail")){
  if(length(stat_list[[g]]$user.Sample_Name) > 0){
    
    # List of relevant images
    im_list <- stat_list[[g]]$user.Sample_Name
    
    # Highlight problematic spots
    filename <- paste0(output_path, "/spatial_", g, "_", sample_name)
    p <- batch_spatial_feature_plot_my(list(Cohort = data_seurat), c(g), output_file=filename, title=g, plot_type=c("d"), image_list=im_list)
    
    # List pathology annotations
    filename <- paste0(output_path, "/spatial_", g, "_pathgroups_", sample_name)
    p <- batch_spatial_feature_plot_my(list(Cohort = data_seurat), c("Pathology.Group"), output_file=filename, title=g, plot_type=c("d"), image_list=im_list)
    
    # Show relevant annotations
    for(a in c(pathology_name, sig_names)){
      
      # List pathology annotations
      filename <- paste0(output_path, "/spatial_", g, "_", a, "_", sample_name)
      p <- batch_spatial_feature_plot_my(list(Cohort = data_seurat), c(a), output_file=filename, title=g, plot_type=c("f"), image_list=im_list)
      
    }
  }
}


## INVESTIGATE SIGNATURE FAILURES ----


# Pull outlier coordinates
barcode_list <- meta_df %>%
  dplyr::filter(SignatureFail == TRUE) %>%
  dplyr::pull(Coordinate)

# Subset Seurat object
data_subset_seurat <- subset(data_seurat, cells=barcode_list)

# Check if SignatureFail spots have unusually low library size
filename <- paste0(output_path, "/box_signaturefail_v_libsize_", sample_name)
create_box_plot_my(meta_df, x_label="SignatureFail", y_label="nFeature_Spatial", fill_label="SignatureFail",
                   filename=filename)

# Find the properties of the genes in the target signature only within the SignatureFail spots
for(s in names(signature_list)){
  sig_gene_stat_list[[paste0(s, "_SigFail")]] <- seurat_gene_properties_my(data_subset_seurat, signature_list[[s]], assay=assay, sct_threshold=sct_threshold) %>%
    dplyr::mutate(Signature_Name = paste0(s, "_SigFail"))
  
  # Write stats to file
  filename <- paste0(output_path, "/table_signaturefail_by_gene_", sample_name, "_", s, ".txt")
  readr::write_delim(sig_gene_stat_list[[paste0(s, "_SigFail")]], filename, delim="\t")
}

# Stack all gene stats into a tibble
sig_gene_stat_df <- dplyr::bind_rows(sig_gene_stat_list)

# Check individual gene levels within SigFail spots compared to the whole dataset
filename <- paste0(output_path, "/bar_signaturefail_genes_", sample_name)
create_bar_plot_my(sig_gene_stat_df, x_label="Symbol", y_label="ExpressionMean", fill_label="Signature_Name",
                   facet_var=c("Signature_Name", "fixed"), filename=filename, error_label="ExpressionSD", reorder_x=TRUE)
