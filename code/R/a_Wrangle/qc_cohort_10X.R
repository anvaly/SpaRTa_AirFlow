# Author: Anna Lyubetskaya. Date: 20-02-04

# Create QC visualizations using 10X reports

# Gather and visualize 10X QC data across mutliple samples
# Samples should be accompanied by a sample_names.txt meta data file containing the following columns:
# Sample_ID (matching folder names), Sample_Name (user-friendly names), FullPath (full path to data folders), meta data columns for groupings

# Note:
# 1. 10X introduced a new QC metric for Visium between spaceranger-1.1.0 and spaceranger-1.2.1: "Mean Reads Under Tissue per Spot" (in addition to "Mean Reads per Spot")
# New new metric is the only one that makes sense.
# Mean Reads per Spot = The number of reads, both under and outside of tissue, divided by the number of barcodes associated with a spot under tissue.
# Mean Reads Under Tissue per Spot = The number of reads under tissue divided by the number of barcodes associated with a spot under tissue.


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_in_out.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_ggplot.R")


## PARAMETERS ----


meta_col_x <- "Block_ID"  # "Cohort"
meta_col_f <- "Block_ID"  # "Set"

# Flag to parse technology-specific fields
technology <- "Visium_probes"  # Visium_polyA, Visium_probes, Perturb, other


## PATHS ----


# File with a list of samples
sample_file <- "XXXX"

# Read sample file
meta_df <- read_file2df_my(sample_file, delim="\t")

# QC file name
qc_file_name <- "metrics_summary.csv"

# A list of full paths to summary metrics file
qc_files <- paste0(meta_df$Pipeline.FullPath, "/", qc_file_name)

# Output folder
output_folder <- "XXXX"
dir.create(output_folder, showWarnings=FALSE)


## Setup plotting groups ----


# QC fields common across all 10X technologies we've seen so far
qc_field_groups <- list("Totals" = c("Number of Reads", "Genes Detected"), 
                        "Reads Mapped Confidently" = c("Reads Mapped Confidently to Intergenic Regions", "Reads Mapped Confidently to Intronic Regions", "Reads Mapped Confidently to Exonic Regions"), 
                        "Reads Mapped" = c("Reads Mapped to Genome", "Reads Mapped Confidently to Genome", "Reads Mapped Confidently to Transcriptome", "Reads Mapped Antisense to Gene"),
                        "Seq Checks" = c("Valid Barcodes", "Sequencing Saturation"),
                        "Q30 Bases" = c("Q30 Bases in Barcode", "Q30 Bases in UMI", "Q30 Bases in RNA Read"))

# QC fields specific to Visium technology; just look for the word "Spot" in the field name
if(grepl("Visium", technology)){
  qc_field_groups <- c(qc_field_groups, list("Spot Totals" = c("Number of Spots Under Tissue", "Fraction of Spots Under Tissue", "Fraction Reads in Spots Under Tissue"), 
                                             "Counts Per Spot" = c("Median Genes per Spot", "Median UMI Counts per Spot", "Mean Reads per Spot", "Mean Reads Under Tissue per Spot", "Valid UMIs")))
  
  if(technology == "Visium_probes"){
    qc_field_groups[["Reads Mapped Confidently"]] <- c("Reads Mapped Confidently to the Filtered Probe Set", "Reads Mapped Confidently to Probe Set")
    qc_field_groups[["Reads Mapped"]] <- c("Reads Mapped to Probe Set")
    qc_field_groups[["Q30 Bases"]] <- c("Q30 Bases in Barcode", "Q30 Bases in UMI", "Q30 Bases in Probe Read")
  }
} else{
  qc_field_groups <- c(qc_field_groups, list("Cell Totals" = c("Estimated Number of Cells", "Fraction Reads in Cells"), 
                                             "Counts Per Cell" = c("Median Genes per Cell", "Median UMI Counts per Cell", "Mean Reads per Cell")))
}

# QC fields specific to Perturb-Seq technology; just look for the word "CRISPR" in the field name
if(technology == "Perturb"){
  qc_field_groups <- c(qc_field_groups, list("CRISPR Totals" = c("CRISPR: Number of Reads", "CRISPR: Mean Reads per Cell", "CRISPR: Median UMIs per Cell (summed over all recognized protospacers)"), 
                                             "CRISPR Protospacers" = c("CRISPR: Cells with 1 or more protospacers detected", "CRISPR: Cells with 2 or more protospacers detected"), 
                                             "CRISPR Fractions" = c("CRISPR: Fraction Guide Reads", "CRISPR: Fraction Guide Reads Usable", 
                                                                    "CRISPR: Fraction Protospacer Not Recognized", "CRISPR: Fraction Reads with Putative Protospacer Sequence"),
                                             "CRISPR Guide" = c("CRISPR: Guide Reads in Cells", "CRISPR: Guide Reads Usable per Cell"),
                                             "CRISPR Seq Checks" = c("CRISPR: Valid Barcodes", "CRISPR: Sequencing Saturation"),
                                             "CRISPR Q30 Bases" = c("CRISPR: Q30 Bases in Barcode", "CRISPR: Q30 Bases in UMI", "CRISPR: Q30 Bases in CRISPR Read")))
}


## INGEST DATA ----


# Read in all QC files
qc_wide_df <- read_dir2file_my(qc_files, delim=",") %>%
  dplyr::mutate(File = gsub("/metrics_summary.csv", "", File)) %>%
  dplyr::mutate_all(~ gsub("%", "", .))

# Turn columns numeric
col_list <- setdiff(colnames(qc_wide_df), c("File", "Sample ID"))
qc_wide_df <- qc_wide_df %>%
  dplyr::mutate_at(vars(dplyr::all_of(col_list)), as.numeric)

if(technology == "Visium"){
  qc_wide_df <- qc_wide_df %>%
    dplyr::mutate(`Sample ID` = gsub("_spaceranger_110", "", `Sample ID`))
}


## WRANGLE DATA ----


# Merge meta data with QC data
qc_wide_df <- meta_df %>%
  dplyr::inner_join(qc_wide_df, by=c("Pipeline.FullPath" = "File")) %>%
  dplyr::select(-Pipeline.FullPath)

# Make QC easier read by switching fractions to percentages
for(c in colnames(qc_wide_df)){
  min_val <- min(qc_wide_df[[c]])
  max_val <- max(qc_wide_df[[c]])
  
  if(!is.na(min_val) & !is.na(max_val) & min_val <= 1 & max_val <= 1 & min_val >= 0 & max_val >= 0){
    qc_wide_df[[c]] <- qc_wide_df[[c]] * 100
  }
}

# Round all numbers
qc_wide_df <- qc_wide_df %>%
  dplyr::mutate_if(is.numeric, round, 1)

# Write joint QC data to file
readr::write_delim(qc_wide_df, paste0(output_folder, "stats.txt"), delim = "\t")


## VISUALIZE DATA ----


# Go through QC field groups and bar plot them
for(metric in unname(unlist(qc_field_groups))){

  filename <- paste0(output_folder, "box_single_", gsub(" ", "", metric))
  labels <- c(meta_col_x, metric, metric)
  
  create_box_plot_my(qc_wide_df, x_label=meta_col_x, y_label=metric, fill_label=meta_col_f,
                     labels=labels, with_dots=TRUE, filename=filename, reorder_x=FALSE)
}


# Go through QC field groups and bar plot them
for(group in names(qc_field_groups)){
  
  # Select columns and transform QC from wide to long tibble
  qc_filt_df <- qc_wide_df %>%
    dplyr::select(dplyr::all_of(c("Sample_Name", qc_field_groups[[group]]))) %>%
    df_wide2long_my(key="Metric", val="Value") %>%
    dplyr::arrange(Sample_Name) %>%
    dplyr::inner_join(meta_df, by="Sample_Name") %>%
    dplyr::mutate(Metric = gsub(" \\(.+\\)", "", Metric))
  
  filename <- paste0(output_folder, "box_combo_", gsub(" ", "", group))
  labels <- c("Sample", "QC Metric", group)

  create_box_plot_my(qc_filt_df, x_label=meta_col_x, y_label="Value", fill_label=meta_col_f, 
                     facet_var=c("Metric", "free_y"), labels=labels, with_dots=TRUE,
                     filename=filename, reorder_x=FALSE)
}
