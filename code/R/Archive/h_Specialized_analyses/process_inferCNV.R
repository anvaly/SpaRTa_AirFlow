# Author: Anna Lyubetskaya. Date: 21-01-13


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_tibble.R")


## PARAMETERS ----


## PATHS ----


# Input folder
input_path <- "XXXX"

# Syngeneics CNV reference
ref_cnv_file <- "XXXX"

# Output folder
output_path <- "XXXX"


## INGEST DATA ----


infercnv_obj <- readRDS(input_path)
ref_df <- readr::read_delim(ref_cnv_file, delim="\t")


## WRANGLE DATA ----


# Extract inferCNV data
data_df <- tibble::as_tibble(infercnv_obj@expr.data, rownames="Symbol") %>%
  df_wide2long_my(key="Sample", val="Score") %>%
  dplyr::mutate(Group = gsub("\\..+$", "", Sample))

# Find inferCNV mean and sd
data_score_mean <- mean(data_df$Score)
data_score_sd <- sd(data_df$Score)

# Summarize inferCNV score by group and gene
summary_df <- data_df %>%
  dplyr::group_by(Group, Symbol) %>%
  dplyr::summarise(Score_mean = mean(Score)) %>%
  dplyr::filter(Score_mean >= data_score_mean + data_score_sd*2 |
                  Score_mean <= data_score_mean - data_score_sd*2) %>%
  dplyr::mutate(Change = ifelse(Score_mean >= data_score_mean + data_score_sd*2, "Gain", "Loss"))

intersect(ref_df$Gene, summary_df$Symbol)


## VISUALIZE DATA ----

