# Author: Anna Lyubetskaya. Date: 20-11-09
# Parse and plot 10X CRISPR output described here:
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/crispr


## ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_in_out.R")
source("code/utils/utils_heatmap.R")
source("code/R/Perturb-Seq/qc_perturb_utils.R")


## PARAMETERS ----


# Project ID
project_id <- "P-20221024-0001"
# Cohort name
cohort_regex <- "MC38_Glides"


## PATHS ----


# Input folder
input_folder <- "XXXX"

# Metadata Sheet
metadata_file <- paste0(input_folder, "PerturbSeq_Metadata_Sheet.txt")


# Input file list
file_list <- c("crispr_analysis/perturbation_efficiencies_by_feature.csv",
               "crispr_analysis/perturbation_efficiencies_by_target.csv",
               "crispr_analysis/protospacer_umi_thresholds.csv",
               "crispr_analysis/protospacer_calls_per_cell.csv", 
               "crispr_analysis/protospacer_calls_summary.csv")

# Output folder
output_folder <- "XXXX"
dir.create(output_folder)


## INGEST META DATA ----


# Read sample file
meta_df <- read_file2df_my(metadata_file, delim="\t") %>%
  tidyr::drop_na() %>%
  dplyr::filter(grepl(cohort_regex, Cohort))


## ANALYZE METRICS ----
for(i in 1:dim(meta_df)[1]){
  sample_name <- meta_df[[i, "Sample_Name"]]
  cohort_name <- meta_df[[i, "Cohort"]]
  data_dir <- meta_df[[i, "FullPath"]]
  filename_in <- paste0(data_dir, "metrics_summary.csv")
  if(i == 1){
    metrics_df <- readr::read_csv(file = filename_in) %>%
      dplyr::mutate_all(dplyr::funs(stringr::str_replace(., "%", ""))) %>%
      tibble::add_column(Sample_Name = sample_name, Cohort = cohort_name,  .before = 1)
    
  } else{
    temp_df <- readr::read_csv(file = filename_in)%>%
      dplyr::mutate_all(dplyr::funs(stringr::str_replace(., "%", ""))) %>%
      tibble::add_column(Sample_Name = sample_name, Cohort = cohort_name,  .before = 1)
    metrics_df <- tibble::add_row(metrics_df, temp_df)
  }
}

# Define a metric by which to sort the QC bar plot
data_order <- as.numeric(metrics_df[["CRISPR: Median UMIs per Cell (summed over all recognized protospacers)"]])
names(data_order) <- metrics_df$Sample_Name

# Factorize the sample names to follow the order of the metric
c["Sample_Name"] <- factor(metrics_df[["Sample_Name"]], levels=names(sort(data_order)))

# QC fields common across all 10X technologies we've seen so far
qc_field_groups <- colnames(metrics_df)[3:ncol(metrics_df)]

# Go through QC field groups and bar plot them
for(metric in unname(unlist(qc_field_groups))){
  metrics_df[metric] <- as.numeric(metrics_df[[metric]])
  metrics_df[metric] <- as.numeric(metrics_df[[metric]])
  
  filename <- paste0(output_folder, "box_", metric)
  create_box_plot_my(metrics_df, x_label="Cohort", y_label=metric, fill_label="Cohort", with_dots=TRUE, 
                     filename=filename, labels=c("", metric, metric))
  
  filename <- paste0(output_folder, "bar_", metric)
  create_bar_plot_my(metrics_df, x_label="Sample_Name", y_label=metric, fill_label="Cohort", 
                     filename=filename, labels=c("", metric, metric))
}


## ANALYZE perturbation efficiency data ----


pvalue_threshold <- 0.05

for(f in file_list[1:2]){
  
  output_file <- gsub("crispr_analysis/|.csv", "", f)
  
  # Ingest a standard 10X output file across a cohort
  data_df <- ingest_dataset(f, meta_df)
  data_df <- dplyr::filter(data_df, !grepl("\\|", data_df$Perturbation))
  
  # Generate various barplots for perturbation_efficiencies_by_feature.csv or perturbation_efficiencies_by_target.csv
  crispr_analysis_barplots(data_df, output_folder, output_file, pvalue_threshold=0.05)
  
  # Heatmap meta data
  hm_meta_df <- data_df %>% 
    dplyr::group_by(Perturbation) %>% 
    dplyr::summarise(IsSignificant = min(`p Value`) <= pvalue_threshold)
  
  params <- list(cell_value = "Log2 Fold Change",
                 row_label = "Perturbation", 
                 col_label = "Sample_Name", 
                 distance = "euclidean",
                 row_annotation = c("IsSignificant"),
                 col_annotation = NULL,
                 range = c(-2, -1, 0),
                 colors = c("royalblue4", "lightseagreen", "white"))
  
  # Create an annotated clustergram
  filename <- paste0(output_folder, output_file, "_hm.png")
  create_heatmap_my(data_df, params, col_meta_df=NULL, row_meta_df=hm_meta_df, filename=filename)
}

data_df <- ingest_dataset(file_list[1], meta_df)

data_df <- dplyr::mutate(data_df, IsSignificant = `p Value` <= pvalue_threshold) %>% dplyr::filter(`Target Gene` != "Nlrc5") %>% dplyr::filter(`Target Gene` != "Cd274")
#  ymin = `Log2 Fold Change Lower Bound`, ymax = `Log2 Fold Change Upper Bound`,

ggplot(data_df, aes(x = reorder(`Target Gene`, `Log2 Fold Change`), y = `Log2 Fold Change`,fill = IsSignificant, ymin = `Log2 Fold Change Lower Bound`, ymax = `Log2 Fold Change Upper Bound`))+
  geom_col()+
  # scale_x_discrete("Target Gene", labels = paste0("Gene ", seq(1,10,1)))+
  geom_hline(yintercept = 0, size = 2)+
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## ANALYZE protospacer data ----


## crispr_analysis/protospacer_umi_thresholds.csv


data_df <- ingest_dataset(file_list[3], meta_df) %>%
  dplyr::mutate(`UMI threshold` = ifelse(`UMI threshold` >= 200, 200, `UMI threshold`))

filename <- paste0(output_folder, gsub("crispr_analysis/|.csv", "", file_list[3]), "_hist")
create_hist_plot_my(data_df, x_label="UMI threshold", fill_label="Sample_Name", 
                    facet_var=c("Sample_Name", "free_y"), intercept=10, binwidth=5, 
                    filename=filename, labels=c("UMI threshold", "Number of gRNAs", "Protospacer UMI thresholds"), 
                    add_density=FALSE, log_scale=FALSE)


## crispr_analysis/protospacer_calls_per_cell.csv


data_df <- ingest_dataset(file_list[4], meta_df, col_types=list(readr::col_character(), readr::col_double(), readr::col_character(), readr::col_character())) %>%
  dplyr::mutate(num_features = ifelse(num_features >= 10, 10, num_features),
                symbol = gsub("-\\d+$", "", feature_call))

filename <- paste0(output_folder, gsub("crispr_analysis/|.csv", "", file_list[4]), "_hist1")
create_hist_plot_my(data_df, x_label="num_features", fill_label="Sample_Name", 
                    facet_var=c("Sample_Name", "free_y"), intercept=1, binwidth=1, 
                    filename=filename, labels=c("Number of features", "Number of cells", "Protospacer calls per cell, features"), 
                    add_density=FALSE, log_scale=FALSE)

filename <- paste0(output_folder, gsub("crispr_analysis/|.csv", "", file_list[4]), "_box")
create_box_plot_my(data_df %>%
                     dplyr::filter(num_features == 1) %>%
                     dplyr::mutate(num_umis = as.numeric(num_umis)), 
                   x_label="feature_call", y_label="num_umis", fill_label="symbol", 
                   facet_var=c("Sample_Name", "free_y"), filename=filename)

filename <- paste0(output_folder, gsub("crispr_analysis/|.csv", "", file_list[4]), "_violin")
create_violin_plot_my(data_df %>%
                        dplyr::filter(num_features == 1) %>%
                        dplyr::mutate(num_umis = as.numeric(num_umis)), 
                      x_label="feature_call", y_label="num_umis", fill_label="symbol", 
                      facet_var=c("Sample_Name", "free_y"), filename=filename)



## crispr_analysis/protospacer_calls_summary.csv


data_df <- ingest_dataset(file_list[5], meta_df)

# Factorize the sample names to follow the order of the metric
data_df["Sample_Name"] <- factor(data_df[["Sample_Name"]], levels=names(sort(data_order)))

# Summary of cell numbers by guide assignment category
filename <- paste0(output_folder, gsub("crispr_analysis/|.csv", "", file_list[5]), "_summary_bar")
create_bar_plot_my(data_df %>% 
                     dplyr::filter(median_umis == "None"), position="dodge",
                   x_label="feature_call", y_label="num_cells", fill_label="Sample_Name", filename=filename, 
                   labels=c("Feature group", "Number of cells", "Cell status summary"), reorder_x=TRUE)

# Number of cells per feature called in each sample
filename <- paste0(output_folder, gsub("crispr_analysis/|.csv", "", file_list[5]), "_summary_cell_bar")
create_bar_plot_my(data_df %>%
                     dplyr::filter(median_umis != "None" & !grepl("\\|", feature_call)), position="dodge",
                   x_label="feature_call", y_label="num_cells", fill_label="Sample_Name", filename=filename, 
                   labels=c("Feature group", "Number of cells", "Number of cells per feature"), reorder_x=TRUE)

# Median number of UMIs assigned to each feature called in each sample
filename <- paste0(output_folder, gsub("crispr_analysis/|.csv", "", file_list[5]), "_summary_umi_bar")
create_bar_plot_my(data_df %>%
                     dplyr::filter(median_umis != "None" & !grepl("\\|", feature_call)) %>%
                     dplyr::mutate(median_umis = as.numeric(median_umis)), position="dodge",
                   x_label="feature_call", y_label="median_umis", fill_label="Sample_Name", filename=filename, 
                   labels=c("Feature group", "Median UMIs", "Median UMI per feature per cell"), reorder_x=TRUE)


single_features <- data_df[!(grepl("\\|",  data_df$feature_call)), ]

