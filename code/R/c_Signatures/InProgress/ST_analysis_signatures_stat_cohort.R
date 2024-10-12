# Author: Anna Lyubetskaya. Date: 20-04-22

# Calculate various statistics for a set of signatures in a merged cohot-level Seurat objects
# Stats include: min, mean, sd, max


## SETUP ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`


source("code/utils/utils_signatures.R")
source("code/utils/utils_ggplot.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_heatmap.R")

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_signatures.R")


## PARAMETERS ----


# The cohort of interest regex ID
# Seurat RDS files are tagged as follows
cohort_name <- "PDAC"

# Gene abundance filters
sct_threshold <- 1
spot_threshold <- 5

# Signatures to select
sig_select <- c("PDAC.P19.Fibroblast", "PDAC.moffitt.activatedstroma", "PDAC.moffitt.normalstroma",
                "PDAC.U.Immune.Macrophage", "BMS.WZ.CD8.effector.CuiCell2021", "BMS.WZ.Bcells.GC.CuiCell2021", "BMS.WZ.TFH.CuiCell2021", "BMS.WZ.Plasma", "PDAC.P19.Tcell",
                "PDAC.U.Exocrine.Acinar", "PDAC.P19.Endocrine",
                "Syng.U.State.Hypoxia.metabolism", "BMS.Pathway.IFNg", "BMS.Pathway.IFNa", "BMS.Pathway.TGFB",
                "BMS.MoCR.CEACAM5", "BMS.MoCR.NIK", "BMS.MoCR.DGKA", "BMS.MoCR.DGKZ", "BMS.MoCR.KRAS", "BMS.MoCR.ADAR", "BMS.MoCR.PTPN2",
                "BMS.MoCR.ATF4", "BMS.NFKB2.NIK.pathway2",
                "PDAC.collisson.classical", "PDAC.collisson.quasimesenchymal", "PDAC.collisson.exocrine",
                "PDAC.moffitt.classical", "PDAC.moffitt.basal",
                "PDAC.bailey.progenitor", "PDAC.bailey.squamous", "PDAC.bailey.adex", "PDAC.bailey.immunogenic"
)

# Rename signatures for visualizations
sig_rename <- c("Fibroblast", "Acitvated stroma", "Normal stroma",
                "Macrophage U", "CD8", "B-cell Cui", "CD4 TFH", "Plasma", "T-cell",
                "Acinar", "Endocrine",
                "Hypoxia metabolism", "Pathway IFNg", "Pathway IFNa", "Pathway TGFb",
                "CEACAM5", "NIK", "DGKA", "DGKZ", "KRAS", "ADAR", "PTPN2",
                "Pathway ATF4", "Pathway NFKB2",
                "PDAC classical C", "PDAC basal C", "PDAC exocrine C", 
                "PDAC classical M", "PDAC basal M", 
                "PDAC classical B", "PDAC basal B", "PDAC exocrine B", "PDAC immunogenic B"
)

# Select signature pairs to correlate
sig_pairs <- list(c("PDAC.collisson.exocrine", "PDAC.U.Exocrine.Acinar"),
                  c("PDAC.bailey.adex", "PDAC.U.Exocrine.Acinar"),
                  c("BMS.MoCR.CEACAM5", "PDAC.collisson.classical"),
                  c("BMS.MoCR.CEACAM5", "PDAC.moffitt.classical"),
                  c("BMS.MoCR.CEACAM5", "PDAC.moffitt.basal"))

# A word of focus in signature names
keyword <- "CEACAM5"


## PATHS ----


# Input folder
input_path <- "XXXX"

# Output information
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_name, "/")

# Create output folder
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)

# Path to a signature file
sig_path <- "data/import/Signatures/signatures_220210.txt"


## INGEST DATA ----


# Ingest a merged RDS Seurat objects
data_seurat <- readRDS(input_path)


## SINGLE-SIGNATURE ANALYSIS ----


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
sig_names_upd <- sort(colnames(data_seurat@meta.data)[grep("sig.", colnames(data_seurat@meta.data))])


## Calculate signature stats ----


# Wide tibble of signature scores
sig_wide_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate") %>%
  dplyr::select(dplyr::all_of(c("Coordinate", sig_names_upd)))

# Create a long tibble of signature scores
sig_df <- sig_wide_df %>%
  df_wide2long_my(key="Signature_name", val="Score") %>%
  dplyr::mutate(Sample_Name = gsub("_[^_]+$", "", Coordinate))

# Calculate global signature stats
sig_stat_df <- sig_df %>%
  dplyr::group_by(Signature_name) %>%
  dplyr::summarise(score_min = min(Score),
                   score_mean = mean(Score),
                   score_median = median(Score),
                   score_sd = sd(Score),
                   score_max = max(Score)) %>%
  dplyr::mutate_if(is.numeric, round, 3) %>%
  dplyr::arrange(desc(score_mean), desc(score_max)) %>%
  dplyr::ungroup()

# Calculate global sample stats
sample_stat_df <- sig_df %>%
  dplyr::group_by(Sample_Name) %>%
  dplyr::summarise(spot_num = dplyr::n_distinct(Coordinate)) %>%
  dplyr::ungroup()

# Classify each spot
sig_status_df <- sig_df %>%
  dplyr::inner_join(sig_stat_df, by="Signature_name") %>%
  dplyr::inner_join(sample_stat_df, by="Sample_Name") %>%
  dplyr::mutate(Status = ifelse(Score >= score_median + score_sd, "T", "F"))

# Calculate signature stats for each sample
sig_sample_stat_df <- sig_status_df %>%
  dplyr::group_by(Signature_name, Sample_Name, score_median, score_sd, spot_num) %>%
  dplyr::summarise(spot_above1sd_num = round(sum(Status == "T") / spot_num * 100),
                   score_min_sample = min(Score),
                   score_mean_sample = mean(Score),
                   score_median_sample = median(Score),
                   score_sd_sample = sd(Score),
                   score_max_sample = max(Score)) %>%
  dplyr::ungroup() %>%
  unique() %>%
  dplyr::mutate(Keyword = grepl(keyword, Signature_name))

# Rename signatures for visualizations
if(!is.null(sig_rename)){
  names(sig_rename) <- paste0("sig.", sig_select)
  print(sig_rename)
  
  sig_sample_stat_df[["Signature_name_user"]] <- sig_rename[sig_sample_stat_df$Signature_name]
}

filename <- paste0(output_path, "cohort_joint_table_sig_stats.txt")
readr::write_delim(sig_sample_stat_df, filename, delim="\t")


## VISUALIZE SINGLE-SIGNATURE DATA ----


# Percent of spots with a signature score higher than median + 1SD
filename <- paste0(output_path, "cohort_joint_signature_boxplot_spot_percent")
create_box_plot_my(sig_sample_stat_df, x_label="Signature_name_user", y_label="spot_above1sd_num", fill_label="Keyword", 
                   labels=c("Signature", "% spot in section above threshold", "% spots with signature score higher than median + 1SD"), 
                   filename=filename, reorder_x=TRUE)

# Signature medians
filename <- paste0(output_path, "cohort_joint_signature_boxplot_sig_median")
create_box_plot_my(sig_sample_stat_df, x_label="Signature_name_user", y_label="score_median_sample", fill_label="Keyword", 
                   labels=c("Signature", "Median signature score", "Median score across spots in a sample"), 
                   filename=filename, with_dots=FALSE, outlier_shape=NULL, reorder_x=TRUE)


# Create a heatmap of signature scores by sample
params <- list(cell_value = "spot_above1sd_num",
               row_label = "Signature_name_user", 
               col_label = "Sample_Name", 
               distance = "pearson",
               row_annotation = NULL,
               col_annotation = NULL,
               range = c(0, 5, 10, 15, 20, 30, 40),
               colors = c("white", "#ffe119", "#bfef45", "#3cb44b", "forestgreen", "#4363d8", "#002F6C"))

filename <- paste0(output_path, "cohort_joint_signature_heatmap.png")
create_heatmap_my(sig_sample_stat_df, params, row_list=NULL, col_list=NULL, col_meta_df=NULL, row_meta_df=sig_sample_stat_df, filename=filename)


## SIGNATURE-PAIRS ANALYSIS ----


# Create a wide tibble with each signature score in each spot
sig_score_wide_df <- sig_status_df %>%
  df_long2wide_my(rows="Coordinate", cols="Signature_name", value="Score") %>%
  dplyr::mutate(Sample_Name = gsub("_[^_]+$", "", Coordinate))

# Create a wide tibble with each signature status in each spot
sig_status_wide_df <- sig_status_df %>%
  df_long2wide_my(rows="Coordinate", cols="Signature_name", value="Status") %>%
  dplyr::mutate(Sample_Name = gsub("_[^_]+$", "", Coordinate))


# Calculate pair-wise stats ----

# Calculate sample-wise signature-pair correlations
# Calculate spot occupancy for every pair of signatures
sig_pairs_list <- list()

for(i in 1:length(sig_names_upd)){
  s1 <- sig_names_upd[i]
  for(j in 1:length(sig_names_upd)){
    s2 <- sig_names_upd[j]
    
    if(i < j){
      
      # Count number of spots in each sample where two sigantures co-occurr or occurr without the other signature
      sig_loc1_df <- sig_status_wide_df %>% 
        dplyr::select(dplyr::all_of(c("Coordinate", "Sample_Name", s1, s2))) %>%
        dplyr::rename(S1 = !!rlang::sym(s1), S2 = !!rlang::sym(s2)) %>%
        dplyr::group_by(Sample_Name, S1, S2) %>%
        dplyr::summarise(Count = dplyr::n_distinct(Coordinate)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(Sig_Pair = paste(sig_rename[s1], "-", sig_rename[s2]),
                      Sig1 = sig_rename[s1],
                      Sig2 = sig_rename[s2],
                      Status = paste0(S1, S2),
                      Pvalue = 1,
                      R = 0) %>%
        dplyr::select(-S1, -S2) %>%
        dplyr::arrange(Status)
      
      # Extract 4 counts for the Fisher's exact test for each sample and signature and calculate Fisher's exact p-value
      # Calculate signature-pair correlation coefficient for each sample
      for(sample in unique(sig_loc1_df$Sample_Name)){
        # Fisher's exact
        index <- which(sig_loc1_df$Sample_Name == sample)
        
        sig_loc1_df[index, "Pvalue"] <- fisher.test(matrix(unname(unlist(sig_loc1_df[index, "Count"])), nrow=2, ncol=2))$p
        
        # Correlation coefficient
        sig_loc2_df <- sig_score_wide_df %>% 
          dplyr::filter(Sample_Name == sample)
        
        sig_loc1_df[index, "R"] <- cor(sig_loc2_df[[s1]], sig_loc2_df[[s2]])
      }
      
      sig_pairs_list[[paste0(s1, s2)]] <- sig_loc1_df
    }
  }
}


# Tibble of pairwise spot occupancies
sig_pairs_df <- dplyr::bind_rows(sig_pairs_list)

# Number of tests
ntests <- nrow(unique(sig_pairs_df[c("Sample_Name", "Sig_Pair")]))

# Count the number of spots in each sample for each signature that have at least one of two signatures "present
sample_stat2_df <- sig_pairs_df %>%
  dplyr::filter(Status != "FF") %>%
  dplyr::group_by(Sample_Name, Sig_Pair) %>%
  dplyr::summarise(spot_num_pos = sum(Count)) %>%
  dplyr::ungroup()

# Bonferroni correction on the Fisher's exact p-values
sig_pairs_df <- sig_pairs_df %>%
  dplyr::mutate(Pvalue = Pvalue * ntests) %>%
  dplyr::inner_join(sample_stat_df, by="Sample_Name") %>%
  dplyr::inner_join(sample_stat2_df, by=c("Sample_Name", "Sig_Pair")) %>%
  dplyr::mutate(PercentSpot = round(Count / spot_num * 100, 1),
                PercentPos = round(Count / spot_num_pos * 100, 1))

# Write signature-pair data to file
filename <- paste0(output_path, "cohort_joint_table_sig_pairs_stats.txt")
readr::write_delim(sig_pairs_df, filename, delim="\t")


# Correlated pairs
sig_pair_cor_select_df <- sig_pairs_df %>%
  dplyr::filter(abs(R) >= 0.25) %>%
  dplyr::group_by(Sig_Pair) %>%
  dplyr::summarise(Count = dplyr::n_distinct(Sample_Name))

# Co-occurring pairs
sig_pair_ocurr_select_df <- sig_pairs_df %>%
  dplyr::filter(Pvalue <= 0.05 & Status == "TT" & Count >= 10) %>%
  dplyr::group_by(Sig_Pair) %>%
  dplyr::summarise(Count = dplyr::n_distinct(Sample_Name),
                   CountVal = dplyr::n_distinct(PercentPos)) %>%
  dplyr::filter(CountVal > 1)


## VISUALIZE SIGNATURE-PAIR STATS ----


# Find percent occupancy for reference signatures
meta_df <- sig_sample_stat_df %>%
  dplyr::filter(Signature_name_user %in% c("PDAC classical C", "PDAC basal", "Acitvated stroma")) %>%
  df_long2wide_my(rows="Sample_Name", cols="Signature_name_user", value="spot_above1sd_num")

# Heatmap of correlation coefficients
for(s in unique(sig_pairs_df$Sig1)){

  # Find percent occupancy for the signature in question for every sample as well as reference signatures
  meta_sample_df <- sig_sample_stat_df %>%
    dplyr::filter(Signature_name_user == s) %>%
    dplyr::select(Sample_Name, spot_above1sd_num) %>%
    dplyr::inner_join(meta_df, by="Sample_Name") %>%
    dplyr::mutate(Sample_Name = gsub("HumanPanc_", "", Sample_Name))
    
  # Create a heatmap of sample-wise signature-pair correlation coefficicents
  params <- list(cell_value = "R",
                 row_label = "Sample_Name", 
                 col_label = "Sig_Pair", 
                 distance = "pearson",
                 row_annotation = c("spot_above1sd_num", "PDAC classical C", "PDAC basal", "Acitvated stroma"),
                 col_annotation = NULL,
                 range = c(-0.5, 0, 0.5),
                 colors = c("#002F6C", "white", "darkred"))
  
  # Select samples and signatures for the correlation heatmap
  sig_loc_df <- sig_pairs_df %>%
    dplyr::filter(Sig_Pair %in% sig_pair_cor_select_df$Sig_Pair &
                    (Sig1 == s | Sig2 == s)) %>%
    dplyr::mutate(Sample_Name = gsub("HumanPanc_", "", Sample_Name))
  
  if(length(unique(sig_loc_df$Sig_Pair)) > 1){
    filename <- paste0(output_path, "cohort_heatmap_cor_", s, ".png")
    create_heatmap_my(sig_loc_df, params, row_list=NULL, col_list=NULL, col_meta_df=NULL, row_meta_df=meta_sample_df, filename=filename)
  }
  
}

# Meta data annotation
meta_df <- sig_pairs_df %>%
  dplyr::filter()

# Heatmap of co-ocurrence rate (%)
for(s in unique(sig_pairs_df$Sig1)){
  
  # Select data for the heatmap
  sig_loc_df <- sig_pairs_df %>%
    dplyr::filter(Sig_Pair %in% sig_pair_ocurr_select_df$Sig_Pair &
                    (Sig1 == s | Sig2 == s) & 
                    Status == "TT") %>%
    dplyr::filter(PercentPos > 0) %>%
    dplyr::mutate(CountLog2 = round(log2(Count + 1), 1))
  
  # Create a heatmap of signature pair co-occurrences by sample
  params <- list(cell_value = "PercentPos",
                 row_label = "Sample_Name",
                 col_label = "Sig_Pair",
                 distance = "pearson",
                 row_annotation = c("CountLog2"),
                 col_annotation = NULL,
                 range = c(0, 5, 10, 15, 20, 30, 40),
                 colors = c("white", "#ffe119", "#bfef45", "#3cb44b", "forestgreen", "#4363d8", "#002F6C"))
  
  if(length(unique(sig_loc_df$Sig_Pair)) > 2){
    filename <- paste0(output_path, "cohort_heatmap_coocurr_", s, ".png")
    create_heatmap_my(sig_loc_df, params, row_list=NULL, col_list=NULL, col_meta_df=NULL, row_meta_df=sig_loc_df, filename=filename)
  }
}


# Create heatmaps of spot intersections between two signatures
if(!is.null(sig_pairs)){
  for(pair in sig_pairs){

    # Name for gene signature pair
    pair_name <- paste(c(sig_rename[[paste0("sig.", pair[1])]], sig_rename[[paste0("sig.", pair[2])]]), collapse=" - ")

    # Select data for the heatmap
    sig_loc_df <- sig_pairs_df %>%
      dplyr::filter(Sig_Pair == pair_name & Status != "FF") %>%
      dplyr::filter(PercentPos > 0) %>%
      dplyr::mutate(Status = gsub("TT", pair_name, Status),
                    Status = gsub("FT", sig_rename[[paste0("sig.", pair[2])]], Status),
                    Status = gsub("TF", sig_rename[[paste0("sig.", pair[1])]], Status))

    # Create a heatmap of signature pair co-occurrences by sample
    params <- list(cell_value = "PercentPos",
                   row_label = "Sample_Name",
                   col_label = "Status",
                   distance = "pearson",
                   row_annotation = c("R"),
                   col_annotation = NULL,
                   range = c(0, round(max(sig_loc_df$PercentPos - 10), -1)),
                   colors = c("white", "#002F6C"))

    filename <- paste0(c(output_path, "cohort_heatmap_pair_coocurr_", paste(pair, collapse="-"), ".png"), collapse="")
    create_heatmap_my(sig_loc_df, params, row_list=NULL, col_list=NULL, col_meta_df=NULL, row_meta_df=sig_loc_df, filename=filename)

  }
}


# Scatter plots of pairs of signatures
if(!is.null(sig_pairs)){
  for(pair in sig_pairs){

    # Prepare a tibble of two signature scores to plot
    sig_loc_df <- sig_status_df %>%
      dplyr::filter(Signature_name %in% paste0("sig.", pair)) %>%
      df_long2wide_my(rows="Coordinate", cols="Signature_name", value="Score") %>%
      dplyr::mutate(Sample_Name = gsub("_[^_]+$", "", Coordinate),
                    Flag = TRUE)

    # Signature lengths and intersection
    l1 <- length(signature_list[[pair[1]]])
    l2 <- length(signature_list[[pair[2]]])
    l3 <- length(intersect(signature_list[[pair[1]]], signature_list[[pair[2]]]))
    # Correlation coefficient
    r <- round(cor(sig_loc_df[[paste0("sig.", pair[1])]], sig_loc_df[[paste0("sig.", pair[2])]]), 3)

    filename <- paste0(c(output_path, "scatter_pair_", paste(pair, collapse="-")), collapse="")
    create_scatter_plot_my(sig_loc_df, x_label=paste0("sig.", pair[1]), y_label=paste0("sig.", pair[2]),
                           fill_label="Flag",
                           # fill_label="Sample_Name", facet_var=c("Sample_Name", "fixed"),
                           filename=filename, do_fit=FALSE, stroke=0, shape=21, size=1,
                           labels=c(pair[1], pair[2],
                                    paste0("SigX L=", l1, ". SigY L=", l2, ". Intersection=", l3, "\n", "R = ", r)))

    filename <- paste0(c(output_path, "scatter_pair_", paste(pair, collapse="-"), "_bysample"), collapse="")
    create_scatter_plot_my(sig_loc_df, x_label=paste0("sig.", pair[1]), y_label=paste0("sig.", pair[2]),
                           fill_label="Flag", facet_var=c("Sample_Name", "fixed"),
                           # fill_label="Sample_Name"
                           filename=filename, do_fit=FALSE, stroke=0, shape=21, size=1,
                           labels=c(pair[1], pair[2],
                                    paste0("SigX L=", l1, ". SigY L=", l2, ". Intersection=", l3, "\n", "R = ", r)))

  }
}
