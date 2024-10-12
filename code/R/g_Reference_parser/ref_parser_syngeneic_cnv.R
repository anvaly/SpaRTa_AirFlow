# Author: Anna Lyubetskaya. Date: 21-01-13

# Parse syngeneics_SM_Mosely_CancerImmRes_2017.xlsx which contains CNV measurements for MC38 and B16
# 1a. Read sheet: Exome-seq_copy_number
# Columns: Sample (filter for MC38 and B16), Gene (edit to upper-case), Log2ratio, and Broad copy changes (values=Gain/Loss)
# 1b. Read sheet: arrayCGH_copy_number
# Columns: Gene_Symbol_Final (edit to upper-case), CT26_CNVfinal, B16F0*_CNVfinal


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

#source("code/utils/utils_tibble.R")
#source("code/utils/utils_gene_lists.R")


## PARAMETERS ----


## PATHS ----


# Input folder
input_path <- "XXXX"

# Output folder
output_file <- "XXXX"


## "syngeneics_SM_Mosely_CancerImmRes_2017.xlsx" ----


# Parse CNV data
data_df <- readxl::read_excel(input_path, sheet="Exome-seq_copy_number", skip=0) %>%
  dplyr::select(Sample, Gene, Chr, Log2ratio, `Broad copy changes`) %>%
  dplyr::mutate(Gene = toupper(Gene)) %>%
  dplyr::rename(Change = `Broad copy changes`) %>%
  dplyr::filter(grepl("B16|MC38", Sample)) %>%
  tidyr::drop_na() %>%
  dplyr::arrange(Sample, Chr, Change, Gene)

# Summarize CNV data
summary_df <- data_df %>%
  dplyr::group_by(Sample, Chr, Change) %>%
  dplyr::summarise(Count = dplyr::n_distinct(Gene)) %>%
  dplyr::arrange(desc(Count))


readr::write_delim(data_df, output_file, delim="\t")
