# Author: Anna Lyubetskaya. Date: 20-02-21


## ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_s3.R")


## PARAMETERS ----


# The work folder

# S3 path information
s3_bucket <- "s3://XXXX"

## The list of all projects run are in ~/data/import/project_list.txt

project_id <- c("P-20211020-0002", "P-20210915-0008", "P-20211108-0001",
                "P-20211213-0006", "P-20220510-0005",
                "P-20211201-0002", "P-20220427-0002", "P-20220228-0002",
                "P-20220322-0008", "P-20220310-0001")


sample_select <- c("HumanPanc_1284133B_ROI1_FFPE_s4_A",
                   "Human_FFPE_PDAC_1255908B_s4_D",
                   "Human_PDAC_FFPE_1255880B_ROI1_s1_C",
                   "Human_PDAC_FFPE_1275233B_ROI1_s1_C",
                   "Human_PDAC_FFPE_E27038_ROI1_s1_A",
                   "Human_FFPE_PDAC_1275301B_s6_D",
                   "Human_FFPE_PDAC_Mega_Pt1_ROI1_s1_A",
                   "Human_FFPE_PDAC_Mega_Pt2_ROI1_s1_B",
                   "Human_FFPE_PDAC_Mega_Pt11_ROI1_s1_A",
                   "Human_FFPE_PDAC_E5058_ROI1_s1_A",
                   "Human_FFPE_PDAC_E1265_ROI1_s1_A")

extensions <- c("fastq.gz")


## S3 TO STASH ----


# Investigate the bucket before copying it
file_list <- list()
for(i in project_id){
  file_list[[i]] <- paste0(s3_bucket, ping_bucket_my(s3_bucket, s3_folder=i, extensions=extensions))
  print(file_list[[i]])
}

file_df <- data.frame(Path = c(unname(unlist(file_list)))) %>%
  dplyr::mutate(Sample = gsub("^.*\\/|_S\\d+_L00\\d+_[RI]\\d+_001.fastq.gz", "", Path)) %>%
  dplyr::filter(Sample %in% sample_select)

setdiff(sample_select, unique(file_df$Sample))

readr::write_delim(file_df, "fastq_list_for_Hannah_230420.txt")
