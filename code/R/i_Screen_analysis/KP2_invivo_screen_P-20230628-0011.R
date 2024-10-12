# Author: Anna Lyubetskaya. Date: 23-07-25
# Analyze KP2 in vivo screen
# https://biodoc.pri.bms.com:8443/display/MGS/P-20230628-0011-2023-Q1-KP2-Invivo-screen

# Data location: s3://XXXX
# Relevant files:
# - X_merged_mage.count.txt
# - X_merged_mage.count_normalized.txt
# - X_merged_mage.countsummary.txt
# - X_merged_read_stat.txt


## ENVIRONMENT ----


# Allow piping throughout the package
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")
source("code/utils/utils_heatmap.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_pca.R")
source("code/utils/utils_stats.R")
source("code/utils/utils_specialized_plots.R")


## PARAMETERS ----


# pDNA
group_pDNA <- c("KP2syng_pDNA")

# 12 doublings of the cell line
group_double12 <- c("KP2syng_12double_r1",
                    "KP2syng_12double_r2",
                    "KP2syng_12double_r3")

# Day zero pre-implantation
group_d0 <- c("KP2syng_preimpl_r1",
              "KP2syng_preimpl_r2",
              "KP2syng_preimpl_r3")

# Syngeneic days 21 and 28
group_d21 <- paste0("KP2syng_d21_m", 1:25)
group_d28 <- paste0("KP2syng_d28_m", setdiff(26:50, c(26, 28, 47)))

# Comparisons
comparison_list <- list(
  ## Controls
  D0_pDNA = list(group_d0, group_pDNA),
  Double12_pDNA = list(group_double12, group_pDNA),
  D0_Double12 = list(group_d0, group_double12),
  ## Comparisons
  D21_pDNA = list(group_d21, group_pDNA),
  D28_pDNA = list(group_d28, group_pDNA),
  D21_D0 = list(group_d21, group_d0),
  D28_D0 = list(group_d28, group_d0),
  D21_Double12 = list(group_d21, group_double12),
  D28_Double12 = list(group_d28, group_double12))

# List of pairs of comparisons to plot
pair_list <- list(c("D0_pDNA", "Double12_pDNA"),
                  c("D21_pDNA", "D28_pDNA"),
                  
                  c("D21_D0", "D28_D0"),
                  c("D21_Double12", "D28_Double12"),
                  
                  c("D21_pDNA", "D21_Double12"),
                  c("D28_pDNA", "D28_Double12"),
                  
                  c("Double12_pDNA", "D28_pDNA"),
                  c("Double12_pDNA", "D21_pDNA"),
                  
                  c("Double12_pDNA", "D28_Double12"),
                  c("Double12_pDNA", "D21_Double12"))

# List of pairs of comparisons to plot
pair_mle_list <- list(c("D21", "D28"),
                      c("D0", "Double12"),
                      c("Double12", "D21"),
                      c("Double12", "D28"))

# Organism
organism <- "mmu" # hsa

# Threshold guides: Number of standard deviations below mean to use to exclude guides
# NULL or 1 or 2
threshold_guides <- 1


## PATHS ----


# MAGeCK MLE design file
mle_design_file <- "XXXX"

# Sample meta data
meta_file <- "XXXX"

# Raw screen data
input_files <- c("XXXX",
                 "XXXX")

# Screen reference data - Brunello library
reference_file <- "XXXX"

# Mouse to human homolog file
homologs_file <- "data/import/References/ensembl_human_mouse_homologs_230111.txt"

# Guide plate layout
plate_file <- "XXXX"

# Common essentials from DepMap
essential_file <- "XXXX"
# Nonessentials from DepMap
nonessential_file <- "XXXX"


# The name of the final RDS object intended for GECO - the name has to be descriptive
output_path <- paste0("XXXX")
output_path_qc <- paste0(output_path, "QC/")
output_path_rra <- paste0(output_path, "RRA/")
output_path_mle <- paste0(output_path, "MLE/")
output_path_result <- paste0(output_path, "Result/")

dir.create(output_path, showWarnings = FALSE)
dir.create(output_path_qc, showWarnings = FALSE)
dir.create(output_path_rra, showWarnings = FALSE)
dir.create(output_path_mle, showWarnings = FALSE)
dir.create(output_path_result, showWarnings = FALSE)

file_log <- paste0(output_path_qc, "log.txt")


## INGEST DATA ----


## Read lib size information


lib_size_pipe_df <- readr::read_delim(libsize_file, delim="\t")


## Read gene informaiton


# Read reference data
control_list <- readr::read_delim(reference_file, delim="\t") %>% 
  dplyr::filter(Group == "Negative Control") %>% 
  dplyr::pull(Gene_Name) %>% 
  unique()

# Read in homolog information
homolog_df <- readr::read_delim(homologs_file, delim="\t")

# Read DepMap essential and nonessential gene lists
# Translate from human to mouse
essential_list <- readr::read_delim(essential_file, delim=",") %>%
  dplyr::mutate(Gene = gsub(" .+", "", Gene)) %>%
  dplyr::inner_join(homolog_df, by=c("Gene" = "Symbol_human")) %>%
  dplyr::rename(Symbol = Symbol_mouse) %>%
  dplyr::pull(Symbol)

nonessential_list <- readr::read_delim(nonessential_file, delim=",") %>%
  dplyr::mutate(Gene = gsub(" .+", "", Gene)) %>%
  dplyr::inner_join(homolog_df, by=c("Gene" = "Symbol_human")) %>%
  dplyr::rename(Symbol = Symbol_mouse) %>%
  dplyr::pull(Symbol)

# Read guide plate information
plate_df <- readxl::read_xlsx(plate_file) %>%
  dplyr::rename(Symbol = Target, sgRNA = Guide) %>%
  dplyr::select(Symbol, sgRNA, Plate, Row, Column) %>%
  dplyr::mutate(Essential = ifelse(Symbol %in% essential_list, "Essential", ""),
                NonEssential = ifelse(Symbol %in% nonessential_list, "", ""),
                TargetingControl = ifelse(Symbol %in% control_list, "Cntrl", ""),
                GeneCategory = paste0(Essential, NonEssential, TargetingControl),
                GeneCategory = ifelse(GeneCategory == "", "Target", GeneCategory),
                Plate = gsub("BR_2023Q1_in-vivo_", "", Plate)) %>%
  dplyr::select(-Essential, -NonEssential, -TargetingControl)


## Read sample informaiton


# Read meta data
meta_df <- readr::read_delim(meta_file, delim="\t") %>%
  dplyr::select(Sample_ID, Sample_Name) %>%
  dplyr::mutate(Double12 = ifelse(Sample_Name %in% group_double12, "Double12", ""),
                D0 = ifelse(Sample_Name %in% group_d0, "D0", ""),
                pDNA = ifelse(Sample_Name %in% group_pDNA, "pDNA", ""),
                D21 = ifelse(Sample_Name %in% group_d21, "Syng_d21", ""),
                D28 = ifelse(Sample_Name %in% group_d28, "Syng_d28", ""),
                SampleCategory = paste0(pDNA, Double12, D0, D21, D28),
                Replicate = gsub(".+_(r\\d+)$|.+_(m\\d+)$", "\\1\\2", Sample_Name)) %>%
  dplyr::select(-pDNA, -Double12, -D0, -D21, -D28)

# Factorize groups to enforce order in vis
meta_df[["SampleCategory"]] <- factor(meta_df[["SampleCategory"]], levels=c("pDNA", "D0", "Double12", "Syng_d21", "Syng_d28"))


## Read data


# Read screen data
screen_init_df <- readr::read_delim(input_files[1], delim="\t") %>%
  # Note: next line is completely specific to this dataset and the current bugs in the production CRISPR pipeline
  dplyr::select(-`KP2_2023Q1_plasmidrepeat-QC`) %>%
  dplyr::rename(Symbol = Gene) %>%
  dplyr::mutate(Symbol = gsub("_sg\\d+", "", sgRNA)) %>%
  dplyr::inner_join(readr::read_delim(input_files[2], delim="\t") %>%
                      dplyr::select(-gRNA_Seq), by=c("sgRNA" = "gRNA_ID", "Symbol" = "Gene_Name"))


## WRANGLE DATA ----


# Sample level stats
lib_size_df <- screen_init_df %>%
  df_wide2long_my(key="Sample_Name", val="Guide", start_col=3) %>%
  dplyr::group_by(Sample_Name) %>%
  dplyr::summarise(LibSize = sum(Guide)) %>%
  dplyr::inner_join(meta_df, by="Sample_Name")

# Library sizes for each sample
lib_sizes <- c(0, 0, colSums(screen_init_df[3:ncol(screen_init_df)]))

# Calculate normalized read counts by adjusting levels to the total library size
screen_init_norm_df <- screen_init_df
for(col in 3:ncol(screen_init_df)){
  screen_init_norm_df[,col] <- log2((as.numeric(c(unlist(screen_init_df[,col]))) + 1) / lib_sizes[[col]] * 1000000)
}

# Generate normalized screen data
screen_init_norm_df <-  screen_init_norm_df %>%
  dplyr::inner_join(plate_df, by=c("sgRNA", "Symbol"))

table(screen_init_norm_df$GeneCategory)
table(screen_init_norm_df$Plate)

# Rename sample columns
for(i in 1:nrow(meta_df)){
  colnames(screen_init_df) <- gsub(meta_df[i,1], meta_df[i,2], colnames(screen_init_df))
  colnames(screen_init_norm_df) <- gsub(meta_df[i,1], meta_df[i,2], colnames(screen_init_norm_df))
}


# Identify guides with unusually low levels in pDNA sample
if(!is.null(threshold_guides)){
  guide_threshold <- mean(screen_init_norm_df$KP2syng_pDNA) - sd(screen_init_norm_df$KP2syng_pDNA) * threshold_guides
  guide_list <- screen_init_norm_df$sgRNA[which(screen_init_norm_df$KP2syng_pDNA < guide_threshold)]
  
  write("GUIDES REMOVED:", file_log, append=TRUE)
  write(paste(guide_list, collapse="; "), file_log, append=TRUE)
  readr::write_delim(as.data.frame(sort(table(gsub("_sg\\d+", "", guide_list)))), file_log, append=TRUE)
  
  # Apply filter
  screen_df <- screen_init_df %>%
    dplyr::filter(!sgRNA %in% guide_list)
  
  screen_norm_df <- screen_init_norm_df %>%
    dplyr::filter(!sgRNA %in% guide_list)
  
} else{
  screen_df <- screen_init_df
  screen_norm_df <- screen_init_norm_df
}


# Write the table fit for MAGECK analysis to file
mageck_input_file <- paste0(gsub("QC/", "", output_path_qc), "mageck_input.txt")
readr::write_delim(screen_df, mageck_input_file, delim="\t")

# Write the table with normalized counts to file
filename <- paste0(output_path_qc, "count_norm_wide.txt")
readr::write_delim(screen_norm_df, filename, delim="\t")


# Make a tibble for vis
screen_norm_long_df <- screen_norm_df %>%
  df_wide2long_my(key="Sample_Name", val="GuideNormLog", start_col=3) %>%
  dplyr::mutate(GuideNormLog = as.numeric(GuideNormLog)) %>%
  dplyr::inner_join(meta_df, by="Sample_Name") %>%
  dplyr::inner_join(plate_df, by=c("sgRNA", "Symbol")) %>%
  dplyr::relocate(GuideNormLog, .after = last_col()) %>%
  tidyr::drop_na()

# Write the table of normalized annotated values to file
filename <- paste0(output_path_qc, "count_norm_long_ann.txt")
readr::write_delim(screen_norm_long_df, filename, delim="\t")


# Calculate guide mean and SD by sample group
guide_stat_df <- screen_norm_long_df %>%
  dplyr::group_by(sgRNA, Symbol, Plate, GeneCategory, SampleCategory) %>%
  dplyr::summarise(GuideNormLogMean = round(mean(GuideNormLog), 2),
                   GuideNormLogSD = round(sd(GuideNormLog), 2),
                   GuideNormLogMedian = round(median(GuideNormLog), 2)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(SampleCategory) %>%
  dplyr::mutate(GuideRank = dplyr::min_rank(dplyr::desc(GuideNormLogMean)))

filename <- paste0(output_path_qc, "guide_stat_long.txt")
readr::write_delim(guide_stat_df, filename, delim="\t")


# Create a wide tibble of ranks
guide_stat_wide_df <- guide_stat_df %>%
  df_long2wide_my(rows="sgRNA", cols="SampleCategory", value="GuideRank") %>%
  dplyr::group_by(sgRNA) %>%
  dplyr::mutate(RankDelta = max(D0, Double12, pDNA, Syng_d21, Syng_d28) - min(D0, Double12, pDNA, Syng_d21, Syng_d28))

filename <- paste0(output_path_qc, "guide_rank_wide.txt")
readr::write_delim(guide_stat_wide_df, filename, delim="\t")


## VISUALIZE DATA FOR QC ----


## LIBRARY SIZES


filename <- paste0(output_path_qc, "library_size")
create_violin_plot_my(lib_size_pipe_df %>%
                        dplyr::inner_join(meta_df, by="Sample_ID"), 
                      x_label="SampleCategory", y_label="P-20230628-0011", fill_label="SampleCategory", 
                      filename=filename, labels=c("Sample Category", "Guide Number", "P-20230628-0011"))


## GUIDE LEVELS


filename <- paste0(output_path_qc, "bar_guide_levels_by_genecategory")
p <- create_bar_plot_my(guide_stat_df, 
                        x_label="sgRNA", y_label="GuideNormLogMean", fill_label="GeneCategory", 
                        facet_var=c("SampleCategory", "fixed"),
                        filename=NULL, labels=c("Guide group", "Guide levels, log norm", "Guide levels by group"),
                        error_label="GuideNormLogSD", reorder_x=TRUE, cols=c("blue", "red", "grey"))

if(!is.null(threshold_guides)){
  p <- p  + 
    geom_hline(yintercept=c(guide_threshold), color="black", linetype="dashed", size=0.2)
}

write_plot2file_my(p, filename, num_row=3, num_col=10)


filename <- paste0(output_path_qc, "bar_guide_levels_by_plate")
p <- create_bar_plot_my(guide_stat_df, 
                        x_label="sgRNA", y_label="GuideNormLogMean", fill_label="GeneCategory", 
                        facet_var=c("Plate", "fixed"),
                        filename=NULL, labels=c("Guide group", "Guide levels, log norm", "Guide levels by group"),
                        error_label="GuideNormLogSD", reorder_x=TRUE, cols=c("blue", "red", "grey"))

write_plot2file_my(p, filename, num_row=7, num_col=10)


filename <- paste0(output_path_qc, "box_genecategory_by_samplecategory_levels")
p <- create_box_plot_my(screen_norm_long_df, 
                        x_label="GeneCategory", y_label="GuideNormLog", fill_label="Sample_Name", 
                        facet_var=c("SampleCategory", "fixed"),
                        filename=NULL, labels=c("Guide group", "Guide levels, log norm", "Guide levels by group"),
                        outlier_shape=NA)

write_plot2file_my(p, filename, num_row=5, num_col=3)


filename <- paste0(output_path_qc, "box_plate_by_genecategory_levels")
p <- create_box_plot_my(screen_norm_long_df, 
                        x_label="Plate", y_label="GuideNormLog", fill_label="GeneCategory", 
                        facet_var=c("SampleCategory", "fixed"),
                        filename=NULL, labels=c("Guide group", "Guide levels, log norm", "Guide levels by group"),
                        outlier_shape=NA)

write_plot2file_my(p, filename, num_row=5, num_col=4)


## PLATE MAPS


# Select only pDNA data
screen_norm_long_pdna_df <- screen_norm_long_df %>%
  dplyr::filter(Sample_Name == "KP2syng_pDNA")

# Establish range of the normalized guide levels (log2)
range_min <- round(min(screen_norm_long_pdna_df$GuideNormLog), 1)
range_mean <- round(mean(screen_norm_long_pdna_df$GuideNormLog), 1)
range_max <- round(max(screen_norm_long_pdna_df$GuideNormLog), 1)
range_sd <- round(sd(screen_norm_long_pdna_df$GuideNormLog), 1)

for(plate in unique(screen_norm_long_pdna_df$Plate)){
  
  # Create a heatmap of signature scores by sample
  params <- list(cell_value = "GuideNormLog",
                 row_label = "Row", 
                 col_label = "Column", 
                 distance = "pearson",
                 range = c(range_min, range_mean, range_max),
                 colors = c("royalblue4", "white", "red3"),
                 row_order = sort(unique(screen_norm_long_pdna_df$Row)),
                 column_order = as.character(sort(as.numeric(unique(screen_norm_long_pdna_df$Column)))),
                 cell_text = screen_norm_long_pdna_df %>%
                   dplyr::filter(Plate == plate) %>% 
                   df_long2wide_my(rows="Row", cols="Column", value="sgRNA") %>%
                   dplyr::select(c("Row", dplyr::all_of(as.character(sort(as.numeric(unique(screen_norm_long_pdna_df$Column)))))))
  )
  
  filename <- paste0(output_path_qc, "hm_", plate, ".png")
  create_heatmap_my(screen_norm_long_pdna_df %>%
                      dplyr::filter(Plate == plate), 
                    params, row_list=NULL, col_list=NULL, col_meta_df=meta_df, row_meta_df=NULL, filename=filename,
                    width=length(unique(screen_norm_long_pdna_df$Column))/2, height=length(unique(screen_norm_long_pdna_df$Row))/2)
  
  filename <- paste0(output_path_qc, "hm_", plate, ".txt")
  readr::write_delim(params$cell_text, filename, delim="\t")
}


## PCA


# PCA parameters
params <- list(sample_value = "Sample_Name", 
               sample_filter = unique(screen_norm_long_df$Sample_Name),
               feature_value = "sgRNA", 
               feature_filter = unique(screen_norm_long_df$sgRNA),
               cell_value = "GuideNormLog", 
               color_by = c("SampleCategory"))

# Perform PCA on the whole cohort
filename <- paste0(output_path_qc, "pca_of_samples")
pc_df <- pca_from_long_tibble_my(screen_norm_long_df, meta_df, 
                                 params, filename=filename, table2file=TRUE, extra_plot=TRUE, dot_labels=NULL)


# PCA parameters
params <- list(sample_value = "sgRNA", 
               sample_filter = unique(screen_norm_long_df$sgRNA),
               feature_value = "Sample_Name", 
               feature_filter = unique(screen_norm_long_df$Sample_Name),
               cell_value = "GuideNormLog", 
               color_by = c("GeneCategory"))

# Perform PCA on the whole cohort
filename <- paste0(output_path_qc, "pca_of_guides_by_genecategory")
pc_df <- pca_from_long_tibble_my(screen_norm_long_df, plate_df, 
                                 params, filename=filename, table2file=FALSE, extra_plot=FALSE, dot_labels=NULL)


# PCA parameters
params <- list(sample_value = "sgRNA", 
               sample_filter = unique(screen_norm_long_df$sgRNA),
               feature_value = "Sample_Name", 
               feature_filter = unique(screen_norm_long_df$Sample_Name),
               cell_value = "GuideNormLog", 
               color_by = c("Plate"))

# Perform PCA on the whole cohort
filename <- paste0(output_path_qc, "pca_of_guides_by_plate")
pc_df <- pca_from_long_tibble_my(screen_norm_long_df, plate_df, 
                                 params, filename=filename, table2file=FALSE, extra_plot=FALSE, dot_labels=NULL)


## CORRELATION


# Output file name for the pairwise correlation analysis
filename <- paste0(output_path_qc, "corr_guides")

# Calculate pairwise correlation between samples
corr_wide_df <- correlation_calculate_plot_my(screen_norm_long_df %>% 
                                                df_long2wide_my(rows="Sample_Name", cols="sgRNA", value="GuideNormLog"), 
                                              row_col="Sample_Name", output_file=filename)

correlation_plot_my(corr_wide_df, scale=c(0, 0.5, 0.25), cols=c("white", "darkblue", "orange"), filename=filename, rowname_col="term")

# Long correlation table
corr_long_df <- corr_wide_df %>%
  dplyr::rename(sgRNA1 = term) %>%
  df_wide2long_my(key="sgRNA2", val="R") %>%
  dplyr::mutate(Symbol1 = gsub("_sg\\d+", "", sgRNA1),
                Symbol2 = gsub("_sg\\d+", "", sgRNA2),
                Self = Symbol1 == Symbol2,
                R = round(R, 3)) %>%
  tidyr::drop_na()

filename <- paste0(output_path_qc, "box_guide_correlation")
p <- create_violin_plot_my(corr_long_df, x_label="Self", y_label="R", fill_label="Self", 
                           filename=filename, labels=c("Do guides belong to the same target?", "Correlation coefficient across samples between two guides", ""))

filename <- paste0(output_path_qc, "corr_guides_long.txt")
readr::write_delim(corr_long_df, filename, delim="\t")


## RUN MAGECK RRA ----


for(run_name in names(comparison_list)){
  
  # Define groups
  target_col <- comparison_list[[run_name]][[1]]
  control_col <- comparison_list[[run_name]][[2]]
  
  # https://sourceforge.net/p/mageck/wiki/usage/
  # mageck test -k sample.txt -t HL60.final,KBM7.final -c HL60.initial,KBM7.initial -n demo
  command_line <- paste("mageck test -k", mageck_input_file, 
                        "-t", paste0(target_col, collapse=","), 
                        "-c", paste0(control_col, collapse=","), 
                        "-n", paste0(output_path_rra, run_name))
  
  print(command_line)
  
  # Run MAGeCK
  system(command_line)
  
  # Run FluteRRA with both gene summary file and sgRNA summary file
  MAGeCKFlute::FluteRRA(paste0(paste0(output_path_rra, run_name, ".gene_summary.txt")), 
                        # paste0(paste0(output_path_rra, run_name, ".guide_summary.txt")), 
                        incorporateDepmap=FALSE,
                        organism=organism, outdir=output_path_rra,
                        proj=run_name, keytype = "Symbol", omitEssential=FALSE)
  
}


## RUN MAGECK MLE ----


run_name <- "MLE"

# https://sourceforge.net/p/mageck/wiki/usage/
# mageck mle -k leukemia.new.csv -d designmat.txt -n beta_leukemia --remove-outliers
command_line <- paste("mageck mle -k", mageck_input_file, 
                      "-d", mle_design_file, 
                      "-n", paste0(output_path_mle, run_name))

print(command_line)

# Run MAGeCK
system(command_line)

# # Run FluteMLE for each pair of comparisons
# for(pair in pair_mle_list){
#   MAGeCKFlute::FluteMLE(paste0(output_path_mle, run_name, ".gene_summary.txt"), keytype = "Symbol",
#                         treatname=pair[[2]], ctrlname=pair[[1]], proj=run_name, organism=organism,
#                         incorporateDepmap=FALSE, omitEssential=FALSE, outdir=output_path_mle)
# }


## LOAD EXTERNAL DEPENDENCY DATA ----


## Manguso KPC result

manguso_file <- "XXXX"

manguso_df <- readr::read_delim(manguso_file, "\t") %>%
  dplyr::filter(Manguso22_KPC_WT_v_Output_nlog10pv > 1.3) %>%
  dplyr::rename(Dependency = Manguso22_KPC_WT_v_Output_log2FC) %>%
  dplyr::select(-Manguso22_KPC_WT_v_Output_nlog10pv) %>%
  dplyr::rename(Symbol = Symbol_mouse) %>%
  dplyr::mutate(Dependency = ifelse(is.na(Dependency), 0, Dependency),
                Negative = ifelse(Dependency <= -0.5, "Neg", ""),
                Positive = ifelse(Dependency >= 0.5, "Pos", ""),
                KPC_Dep = paste0(Negative, Positive)  #, KPC_Dep = ifelse(KPC_Dep == "", "Other", KPC_Dep)
  ) %>%
  dplyr::select(-Negative, -Positive, -Dependency)


# DepMap PDAC dependency information

total_models <- 39 + 11

depmap2d_file <- "XXXX"
depmap3d_file <- "XXXX"

depmap_df <- readr::read_delim(depmap3d_file, "\t") %>%
  dplyr::inner_join(readr::read_delim(depmap2d_file, "\t"), by=c("Symbol" = "Hugo_Symbol")) %>%
  dplyr::select(Symbol, OrgLin_DepMap_SensitivityCount, ClLin_DepMap_SensitivityCount) %>%
  dplyr::mutate(DepMapStatus = paste0("DM:", ClLin_DepMap_SensitivityCount, "/", OrgLin_DepMap_SensitivityCount),
                DepMapStatus = round((ClLin_DepMap_SensitivityCount + OrgLin_DepMap_SensitivityCount) / total_models * 100, -1)) %>%
  dplyr::inner_join(homolog_df, by=c("Symbol" = "Symbol_human")) %>%
  dplyr::select(Symbol_mouse, DepMapStatus) %>%
  dplyr::rename(Symbol = Symbol_mouse) %>%
  dplyr::left_join(manguso_df, by="Symbol") %>%
  dplyr::mutate(KPC_Dep = ifelse(is.na(KPC_Dep), "", KPC_Dep),
                DepStatus = paste0(KPC_Dep, DepMapStatus, "%"))

depmap_df$DepStatus <- factor(depmap_df$DepStatus, levels=c(paste0("Pos", seq(0, 100, 10), "%"),
                                                            paste0(seq(0, 100, 10), "%"),
                                                            paste0("Neg", seq(0, 100, 10), "%")))
table(depmap_df$DepStatus)


## ANALYZE RESULT OF RRA ----


## Load MAGeCK result


# Read through MAGeCK results
result_list <- list()
for(r in names(comparison_list)){
  
  filename <- paste0(paste0(output_path_rra, "MAGeCKFlute_", r, "/RRA/", r, "_processed_data.txt"))
  
  result_df <- readr::read_delim(filename, delim="\t") %>%
    dplyr::select(Symbol, Score, FDR)
  
  colnames(result_df) <- c("Symbol", paste0(r, c("_Score", "_FDR")))
  
  result_list[[r]] <- result_df
  
}

# Read screen result data
result_df <- result_list %>%
  purrr::reduce(dplyr::full_join, by = "Symbol")

# Add DepMap data to the result
result_edit_df <- result_df %>%
  dplyr::left_join(depmap_df, by="Symbol")

# Identify significant hits
result_edit_df[["Significant"]] <- as.character(sapply(1:nrow(result_edit_df), function(x) min(result_edit_df[x, paste0(names(comparison_list), "_FDR")]) <= 0.1))

# Add labels to any standout genes
# result_edit_df[["Label"]] <- as.character(sapply(1:nrow(result_edit_df), function(x) max(abs(result_edit_df[x, paste0(names(comparison_list), "_Score")])) >= 1 &&
#                                             min(abs(result_edit_df[x, paste0(names(comparison_list), "_FDR")])) <= 0.1))
result_edit_df[["Label"]] <- result_edit_df[["Significant"]]

result_edit_df <- result_edit_df %>%
  dplyr::mutate(Label = ifelse(Label == TRUE, Symbol, "")) %>%
  dplyr::left_join(unique(plate_df[c("GeneCategory", "Symbol")]), by="Symbol")

table(result_edit_df$DepStatus)
table(result_edit_df$Label)
table(result_edit_df[["Significant"]])
table(result_edit_df$GeneCategory)


# Write data to file
filename <- paste0(output_path_result, "table_result.txt")
readr::write_delim(result_edit_df, filename, delim="\t")


# Visualize results for every pair
for(pair in pair_list){
  
  # Max value in the data
  value_max <- ceiling(max(abs(result_edit_df[paste0(pair, "_Score")])))
  
  for(f in c("DepStatus", "Significant", "GeneCategory")){
    
    filename <- paste0(output_path_result, "result_", pair[[1]], "_", pair[[2]], "_", f)
    p <- create_scatter_plot_my(result_edit_df, 
                                x_label=paste0(pair[[1]], "_Score"), y_label=paste0(pair[[2]], "_Score"), fill_label=f, 
                                stroke=0, shape=21, size=1.5,
                                filename=NULL, dot_labels="Label", do_fit=TRUE, cols=define_cols_my(n=length(unique(result_edit_df[[f]])))) + 
      geom_vline(xintercept=c(-value_max, 0, value_max), color="black", linetype="dashed", size=0.2) +
      geom_hline(yintercept=c(-value_max, 0, value_max), color="black", linetype="dashed", size=0.2)
    
    write_plot2file_my(p, filename, num_row=1, num_col=1)
    
  }
}


## ANALYZE RESULT OF MLE ----


# List of MLE columns
mle_cols <- paste0(unique(unlist(pair_mle_list)))

# Read MLE result
result_mle_df <- readr::read_delim(paste0(output_path_mle, run_name, ".gene_summary.txt"), delim="\t") %>%
  dplyr::select(c(Gene, dplyr::all_of(c(paste0(mle_cols, "|beta"), paste0(mle_cols, "|fdr"))))) %>%
  dplyr::rename(Symbol = Gene) %>%
  dplyr::left_join(depmap_df, by="Symbol")

# Edit column names
colnames(result_mle_df) <- gsub("\\|beta", "_Score", colnames(result_mle_df))

# Identify significant hits
result_mle_df[["Significant"]] <- as.character(sapply(1:nrow(result_mle_df), function(x) min(result_mle_df[x, paste0(mle_cols, "|fdr")]) <= 0.1))

# Add labels to any standout genes
# result_mle_df[["Label"]] <- as.character(sapply(1:nrow(result_mle_df), function(x) max(abs(result_mle_df[x, paste0(mle_cols, "_Score")])) >= 1 &&
#                                                    min(abs(result_mle_df[x, paste0(mle_cols, "|fdr")])) <= 0.05))
result_mle_df[["Label"]] <- result_mle_df[["Significant"]]

result_mle_df <- result_mle_df %>%
  dplyr::mutate(Label = ifelse(Label == TRUE, Symbol, "")) %>%
  dplyr::left_join(unique(plate_df[c("GeneCategory", "Symbol")]), by="Symbol")


table(result_mle_df$DepStatus)
table(result_mle_df$Label)
table(result_mle_df[["Significant"]])
table(result_mle_df$GeneCategory)

# Visualize results for every pair
for(pair in pair_mle_list){
  
  # Max value in the data
  value_max <- ceiling(max(abs(result_mle_df[paste0(pair, "_Score")])))
  
  for(f in c("DepStatus", "Significant", "GeneCategory")){
    
    filename <- paste0(output_path_result, "mle_result_", pair[[1]], "_", pair[[2]], "_", f)
    p <- create_scatter_plot_my(result_mle_df, 
                                x_label=paste0(pair[[1]], "_Score"), y_label=paste0(pair[[2]], "_Score"), fill_label=f, 
                                filename=NULL, dot_labels="Label", do_fit=TRUE, cols=define_cols_my(n=length(unique(result_edit_df[[f]])))) + 
      geom_vline(xintercept=c(-value_max, 0, value_max), color="black", linetype="dashed", size=0.2) +
      geom_hline(yintercept=c(-value_max, 0, value_max), color="black", linetype="dashed", size=0.2)
    
    write_plot2file_my(p, filename, num_row=1, num_col=1)
    
  }
}


## PLOT SELECT GENES ----



# Extract pDNA data only
pDNA_df <- screen_norm_long_df %>% 
  dplyr::filter(SampleCategory == "pDNA") %>%
  dplyr::select(sgRNA, GuideNormLog) %>%
  dplyr::rename(GuideNormLogBase = GuideNormLog)

# Normalize all guide data to pDNA
screen_norm_pdna_long_df <- screen_norm_long_df %>%
  dplyr::inner_join(pDNA_df, by="sgRNA") %>%
  dplyr::mutate(GuideNormLogpDNA = GuideNormLog - GuideNormLogBase)


write("TARGET HITS:", file_log, append=TRUE)
write(paste(sort(unique(result_edit_df$Label)), collapse="; "), file_log, append=TRUE)

write("TARGET HITS MLE:", file_log, append=TRUE)
write(paste(sort(unique(result_mle_df$Label)), collapse="; "), file_log, append=TRUE)


# Plot boxplot of normalized levels for select genes
for(s in unique(result_edit_df$Label)){
  
  if(s %in% screen_norm_pdna_long_df$Symbol){
    filename <- paste0(output_path_result, "box_normlevels_rra_", s)
    p <- create_box_plot_my(screen_norm_pdna_long_df %>%
                              dplyr::filter(Symbol == s), 
                            x_label="SampleCategory", y_label="GuideNormLogpDNA", fill_label="sgRNA",
                            filename=NULL, labels=c("Sample Group", "Guide levels, log norm", "Guide levels by group"),
                            outlier_shape=NA, with_dots=TRUE)
    
    write_plot2file_my(p, filename, num_row=1, num_col=1)
  }
  
}


# Plot boxplot of normalized levels for select genes
for(s in unique(result_mle_df$Label)){
  
  if(s %in% screen_norm_pdna_long_df$Symbol){
    filename <- paste0(output_path_result, "box_normlevels_mle_", s)
    p <- create_box_plot_my(screen_norm_pdna_long_df %>%
                              dplyr::filter(Symbol == s), 
                            x_label="SampleCategory", y_label="GuideNormLogpDNA", fill_label="sgRNA",
                            filename=NULL, labels=c("Sample Group", "Guide levels, log norm", "Guide levels by group"),
                            outlier_shape=NA, with_dots=TRUE)
    
    write_plot2file_my(p, filename, num_row=1, num_col=1)
  }
  
}
