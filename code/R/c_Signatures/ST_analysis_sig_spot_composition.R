# Author: Anna Lyubetskaya. Date: 23-07-27
# Call cell types in each Visium spot using the simplistic signature scoring method
# Use pathology niches to score signature better


## ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_heatmap.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_signatures.R")


## PARAMETERS ----


# Name of the processed Seurat object
sample_name <- "PDAC108_path14_merge"
remove_string <- ""

# Name for the signature group being plotted
sig_group_name <- "Peng_UCell"

# SCT and spot number thresholds for gene filtering
sct_threshold <- 0.5
spot_threshold <- 5

# Assay and slot to use to calculate signature score
assay <- "SCT"
slot <- "data"

# Name for the signatures to be plotted
sig_select <- c("PDAC.collisson.classical","PDAC.collisson.exocrine","PDAC.collisson.quasimesenchymal",
                "PDAC.moffitt.activatedstroma","PDAC.moffitt.basal","PDAC.moffitt.classical","PDAC.moffitt.normalstroma",
                "PDAC.P19.Acinar","PDAC.P19.Bcell","PDAC.P19.Ductal_1","PDAC.P19.Ductal_2",
                "PDAC.P19.Endocrine","PDAC.P19.Endothelial","PDAC.P19.Fibroblast",
                "PDAC.P19.Macrophage","PDAC.P19.Stellate","PDAC.P19.Tcell",
                "PDAC.CosMx.Mast","PDAC.CosMx.Plasma","PDAC.U.Nervous",
                "BMS.Pathway.IFNa","BMS.Pathway.IFNg",
                "BMS.Pathway.TGFB","BMS.Pathway.TNFa",
                "BMS.CL.Hypoxia", "Syng.U.State.Hypoxia_Hallmark",
                "PDAC.Elyada19.panCAF","PDAC.Elyada19.iCAF","PDAC.Elyada19.myCAF")

# User-defined clustering resolution
resolution <- "Pathology.Group"

# Meta parameters select
# Names of pathology fields to be plotted
meta_select <- NULL


# Which test to use to identify anchor classes
# Diff (diff), KS-test (ks), T-tes (t), or manual
anchor_method <- "manual"

# Number of SDs to threshold signatures
sd_threshold <- 2


## PATHS ----


# File of manually set anchor pathology niches for each signature
anchor_file <- "XXXX"

# Path to a signature file
sig_path <- "data/import/Signatures/signatures_pdac_20230720.txt"

# Location of pre-processed data
input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, sample_name, "_", sig_group_name, "_v2_", anchor_method, "_SD", sd_threshold, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Seurat data
data_seurat <- readRDS(input_path)

# Read anchors if set
if(!is.null(anchor_file)){
  anchors_df <- readr::read_delim(anchor_file) %>%
    tibble::column_to_rownames("Sig_Name")
}


## CALCULATE SIGNATURE SCORES ----


## Calculate signature scores and add them to Seurat meta data

# Find genes abundant in this sample
gene_list <- seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold,
                                             assay=assay, slot=slot)

# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select,
                                            sig_length_min=1, sig_length_max=1000, ratio_threshold=0)

# Add signature scores to a seurat object
data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.", assay=assay)

# Find column names in the Seurat object - Seurat can rename original signatures
sig_names <- colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))])

# Find all genes in signatures and their signature assignations
sig_invert_df <- invert_list_my(signature_list)


## WRANGLE DATA ----


# Extract coordinate meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

# Add user.Sample_Name to meta data if absent: happens in stand-alone samples
if(!"user.Sample_Name" %in% colnames(meta_df)){
  meta_df["user.Sample_Name"] = data_seurat@misc$user.Sample_Name
}

# Crete a long tibble for meta data
meta_long_df <- meta_df %>%
  dplyr::select(dplyr::all_of(c("Coordinate", "user.Sample_Name", resolution, sig_names, meta_select))) %>%
  df_wide2long_my(key="Sig_Name", val="Score", start_col=4)


## THRESHOLD SIGNATURE ----


# Cycle through signatures
# Calculate the threshold using pathology information
# Collect spot status
status_list <- list()
for(sig in c(sig_names, meta_select)){
  
  data_seurat@meta.data[which(is.na(data_seurat@meta.data[sig])), sig] <- 0
  
  # Select data for a particular feature
  meta_long_select_df <- meta_long_df %>%
    dplyr::filter(Sig_Name == sig)
  
  # Mean and SD of the feature distribution by selected group (pathology or clusters)
  meta_sum_df <- meta_long_select_df %>%
    dplyr::group_by(!!rlang::sym(resolution)) %>%
    dplyr::summarize(Mean = round(mean(Score), 4),
                     SD = round(sd(Score), 4)) %>%
    dplyr::arrange(Mean)
  
  # Identify the difference between consequtive elements in a vector of means
  if(anchor_method == "diff"){
    
    meta_sum_diff <- diff(meta_sum_df$Mean)
    max_indx <- which.max(meta_sum_diff)
    # max_indx <- which(meta_sum_diff == Rfast::nth(meta_sum_diff, 2, descending = T))
    meta_sum_df[["Difference"]] <- c(meta_sum_diff, 0)
    anchor_list <- as.character(meta_sum_df$Pathology.Group[1:max_indx])
    
  } else if(anchor_method == "ks"){
    
    ks_test <- c()
    for(i in 1:(nrow(meta_sum_df)-1)){
      ks_test_loc <- ks.test(meta_long_select_df %>%
                               dplyr::filter(!!rlang::sym(resolution) == meta_sum_df[[i+1, resolution]]) %>%
                               # dplyr::filter(!!rlang::sym(resolution) %in% c(unlist(unname(meta_sum_df[i+1:nrow(meta_sum_df), resolution])))) %>%
                               dplyr::pull(Score), 
                             meta_long_select_df %>%
                               dplyr::filter(!!rlang::sym(resolution) == meta_sum_df[[i, resolution]]) %>%
                               # dplyr::filter(!!rlang::sym(resolution) %in% c(unlist(unname(meta_sum_df[1:i, resolution])))) %>%
                               dplyr::pull(Score)
      )
      
      if(ks_test_loc$p.value <= 0.01){
        ks_test[i] <- unlist(unname(ks_test_loc$statistic))
      } else{
        ks_test[i] <- 0
      }
    }
    
    max_indx <- which.max(ks_test)
    # max_indx <- Rfast::nth(ks_test, 2, descending = T)
    meta_sum_df[["Difference"]] <- c(ks_test, 0)
    anchor_list <- as.character(meta_sum_df$Pathology.Group[1:max_indx])
    
  } else if(anchor_method == "t"){
    
    t_test <- c()
    for(i in 1:(nrow(meta_sum_df)-1)){
      t_test_loc <- t.test(meta_long_select_df %>%
                             dplyr::filter(!!rlang::sym(resolution) == meta_sum_df[[i+1, resolution]]) %>%
                             # dplyr::filter(!!rlang::sym(resolution) %in% c(unlist(unname(meta_sum_df[i+1:nrow(meta_sum_df), resolution])))) %>%
                             dplyr::pull(Score),
                           meta_long_select_df %>%
                             dplyr::filter(!!rlang::sym(resolution) == meta_sum_df[[i, resolution]]) %>%
                             # dplyr::filter(!!rlang::sym(resolution) %in% c(unlist(unname(meta_sum_df[1:i, resolution])))) %>%
                             dplyr::pull(Score)
      )
      
      if(t_test_loc$p.value <= 0.01){
        t_test[i] <- unlist(unname(t_test_loc$statistic))
      } else{
        t_test[i] <- 0
      }
    }
    
    max_indx <- which.max(t_test)
    # max_indx <- Rfast::nth(t_test, 2, descending = T)
    meta_sum_df[["Difference"]] <- c(t_test, 0)
    anchor_list <- as.character(meta_sum_df$Pathology.Group[1:max_indx])
  } else if(anchor_method == "manual"){
    anchor_list <- strsplit(anchors_df[gsub("sig.", "", sig),1], ",")[[1]]
  }
  
  # Identify signature mean and SD just in the low signature populations
  meta_long_select_anchor_df <- meta_long_select_df %>% 
    dplyr::filter(Pathology.Group %in% anchor_list)
  
  sig_mean <- mean(meta_long_select_anchor_df$Score)
  sig_sd <- sd(meta_long_select_anchor_df$Score)
  
  
  # Print the anchors
  print(paste0(sig, "    ", paste0(anchor_list, collapse="; ")))
  
  # Signature threshold
  sig_threshold <- sig_mean + sig_sd * sd_threshold
  
  # Some signatures are so over abundant this logic doesn't work at all
  # This is a dumb patch to lower the threshold in these cases
  if(sig_threshold > 1){
    sig_threshold <- sig_mean + sig_sd * (sd_threshold - 1)
  }
  
  # Threshold the signature
  status_list[[sig]] <- meta_long_select_df %>%
    dplyr::mutate(SigStatus = Score >= sig_threshold,
                  NicheAnchorStatus = Pathology.Group %in% anchor_list)
  
  # Calculate pass/fail ratio of the signature
  sig_stat <- table(status_list[[sig]]$SigStatus)
  sig_percent <- round(sig_stat[[2]] / (sig_stat[[2]] + sig_stat[[1]]) * 100, 2)
  
  # Signature histogram
  filename <- paste0(output_path, "/hist_", sig, "_", sample_name)
  create_hist_plot_my(status_list[[sig]], x_label="Score", fill_label="SigStatus", 
                      filename=filename, labels=c(sig, "Score", paste0(sig, ".  Success rate = ", sig_percent, "%")), 
                      log_scale=FALSE, add_density=TRUE,
                      binwidth=0.01, intercept=c(0,1))
  
  # Average signature score vs cluster
  filename <- paste0(output_path, "/box_", sig, "_", sample_name)
  create_box_plot_my(status_list[[sig]], x_label=resolution, y_label="Score", fill_label="NicheAnchorStatus",
                     facet_var=c("Sig_Name", "free_x"), filename=filename, labels=c(resolution, "Score", sig),
                     outlier_shape=NA, reorder_x=TRUE)
  
}

# Join signature call results
status_df <- dplyr::bind_rows(status_list) %>%
  dplyr::mutate(SigStatus = ifelse(SigStatus == TRUE, 1, 0))

# Write signature status information to file
filename <- paste0(output_path, "/table_sig_status_", sig_group_name, ".txt")
readr::write_delim(status_df %>%
                     df_long2wide_my(rows="Coordinate", cols="Sig_Name", value="SigStatus"), 
                   filename, delim="\t")


# Summarize sig status by sample and pathology
status_sample_path_df <- status_df %>%
  dplyr::group_by(user.Sample_Name, Sig_Name, Pathology.Group) %>%
  dplyr::summarise(SpotCount = sum(SigStatus)) %>%
  dplyr::ungroup()

# Write signature status information to file
filename <- paste0(output_path, "/table_sig_status_sample_path_", sig_group_name, ".txt")
readr::write_delim(status_sample_path_df, filename, delim="\t")


# Summarize sig status by sample
status_sample_df <- status_df %>%
  dplyr::group_by(user.Sample_Name, Sig_Name) %>%
  dplyr::summarise(SampleCount = sum(SigStatus)) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(meta_df %>% 
                      dplyr::group_by(user.Sample_Name) %>% 
                      dplyr::summarise(Count = dplyr::n_distinct(Coordinate)), by="user.Sample_Name") %>%
  dplyr::mutate(SamplePercent = round(SampleCount / Count * 100, 1))

# Write signature status information to file
filename <- paste0(output_path, "/table_sig_sample_", sig_group_name, ".txt")
readr::write_delim(status_sample_df, filename, delim="\t")


# Summarize sig status by pathology
status_path_df <- status_df %>%
  dplyr::group_by(Pathology.Group, Sig_Name) %>%
  dplyr::summarise(PathCount = sum(SigStatus)) %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(meta_df %>% 
                      dplyr::group_by(Pathology.Group) %>% 
                      dplyr::summarise(Count = dplyr::n_distinct(Coordinate)), by="Pathology.Group") %>%
  dplyr::mutate(PathPercent = round(PathCount / Count * 100, 1))

# Write signature status information to file
filename <- paste0(output_path, "/table_sig_path_", sig_group_name, ".txt")
readr::write_delim(status_path_df, filename, delim="\t")


## VISUALIZE PER-SPOT SIGNATURE CALLS ----


# Plot pathology classes by signature definition of spot
p <- create_bar_plot_my(status_path_df, x_label="Pathology.Group", y_label="PathPercent", fill_label="Sig_Name", 
                        position="stack", filename=NULL, 
                        labels=c("Pathology Group", "% spots assigned to signature", "Pathology niches by signature")) +
  scale_fill_manual(values=define_cols_my(n=length(sig_names)))

filename <- paste0(output_path, "bar_path_", sig_group_name)
write_plot2file_my(p, filename, 1, 2)


# Plot samples by signature definition of spot
p <- create_bar_plot_my(status_sample_df, x_label="user.Sample_Name", y_label="SamplePercent", fill_label="Sig_Name", 
                        position="stack", filename=NULL, 
                        labels=c("Sample Name", "% spots assigned to signature", "Samples by signature")) +
  scale_fill_manual(values=define_cols_my(n=length(sig_names)))

filename <- paste0(output_path, "bar_sample_", sig_group_name)
write_plot2file_my(p, filename, 1, 5)


# Create a heatmap of signature scores by sample
params <- list(cell_value = "PathPercent",
               row_label = "Sig_Name", 
               col_label = "Pathology.Group", 
               distance = "pearson",
               row_annotation = NULL,
               col_annotation = NULL,
               range = seq(0,50,10),
               colors = RColorBrewer::brewer.pal(n=6, name="PuBuGn"),
               show_column_dend = TRUE,
               show_row_dend = TRUE)

filename <- paste0(output_path, "hm_path_", sig_group_name, ".png")
create_heatmap_my(status_path_df, params, row_list=NULL, col_list=NULL, 
                  col_meta_df=NULL, row_meta_df=NULL, filename=filename,
                  width=9, height=6)


## CO-OCCURRENCE OF SIGNATURES IN SPOTS ----


sig_list <- unique(status_df$Sig_Name)
sig_num <- length(sig_list)
tibble_size <- sig_num*sig_num*length(unique(status_df$Pathology.Group))

count <- 0
cooccurr_list <- list()
for(path in unique(status_df$Pathology.Group)){
  path_num <- status_df %>%
    dplyr::filter(Pathology.Group == path) %>%
    dplyr::pull(Coordinate) %>%
    unique() %>%
    length()
  
  for(i in 1:sig_num){
    sig1_num <- status_df %>%
      dplyr::filter(Pathology.Group == path & Sig_Name == sig_list[[i]] & SigStatus == 1) %>%
      dplyr::pull(Coordinate) %>%
      unique() %>%
      length()
    
    for(j in 1:sig_num){
      count <- count+1
      
      # Calculate number and percent of spots where the two signatures co-occurr for every pathology niceh
      intersect_val <- sum(status_df %>%
                             dplyr::filter(Pathology.Group == path & Sig_Name == sig_list[[i]]) %>%
                             dplyr::pull(SigStatus) +
                             status_df %>%
                             dplyr::filter(Pathology.Group == path & Sig_Name == sig_list[[j]]) %>%
                             dplyr::pull(SigStatus) == 2)
      intersect_perc_path <- round(intersect_val / path_num * 100, 2)
      intersect_perc_sig1 <- round(intersect_val / sig1_num * 100, 2)
      
      cooccurr_list[[count]] <- c(a=path, b=sig_list[[i]], c=sig_list[[j]], d=intersect_val, e=intersect_perc_path, f=intersect_perc_sig1)
    }
  }
}


# Co-occurrence tibble
cooccurr_df <- dplyr::bind_rows(cooccurr_list)
colnames(cooccurr_df) <- c("Pathology.Group", "Sig_Name1", "Sig_Name2", "OccurrenceNum", 
                           "OccurrencePathPerc", "OccurrenceSig1Perc")

# Write co-occurrence information to file
filename <- paste0(output_path, "/table_cooccurr_", sig_group_name, ".txt")
readr::write_delim(cooccurr_df, filename, delim="\t")


# Create a heatmap of co-occurrence
for(path in unique(status_df$Pathology.Group)){
  
  cooccurr_loc_df <- cooccurr_df %>%
    dplyr::filter(Pathology.Group == path) %>%
    dplyr::mutate(OccurrenceSig1Perc = as.numeric(OccurrenceSig1Perc),
                  OccurrenceNum = as.numeric(OccurrenceNum))
  
  cooccurr_loc_meta_df <- cooccurr_loc_df %>%
    dplyr::filter(Sig_Name1 == Sig_Name2)
  
  params <- list(cell_value = "OccurrenceSig1Perc",
                 row_label = "Sig_Name2", 
                 col_label = "Sig_Name1", 
                 distance = "pearson",
                 row_annotation = NULL,
                 col_annotation = c("OccurrenceNum"),
                 range = seq(0,50,10),
                 colors = RColorBrewer::brewer.pal(n=6, name="YlGnBu"),
                 row_order = paste0("sig.", sig_select),
                 column_order = paste0("sig.", sig_select))
  
  filename <- paste0(output_path, "hm_cooccurr_", path, "_", sig_group_name, ".png")
  create_heatmap_my(cooccurr_loc_df, params, row_list=NULL, col_list=NULL, 
                    col_meta_df=cooccurr_loc_meta_df, row_meta_df=NULL, filename=filename,
                    width=11, height=7)
}
