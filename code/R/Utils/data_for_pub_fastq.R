# Author: Anna Lyubetskaya. Date: 20-02-21
# Gather all FASTQ data for a publication
# Standard FASTQ naming structure: [SAMPLE_NAME]_SX_L00Y_[I1/I2/R1/R2]_001.fastq.gz

# From terminal, code to split an RA FASTQ file into R1 and R2 from https://www.biostars.org/p/141256/: 
# cd /XXXX
# tar -xzf BMS_Phase_2.tar --wildcards --no-anchored *.fastq*
# cd /XXXX
# zcat 1126819/read-RA_si-CTCAGCGGGA_lane-001-chunk-001.fastq.gz | paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > 1126819/HumanPanc_1284133B_ROI1_FFPE_s4_A_S1_L001_R1_001.fastq) | cut -f 5-8 | tr "\t" "\n" > 1126819/HumanPanc_1284133B_ROI1_FFPE_s4_A_S1_L001_R2_001.fastq
# zcat 1126819/read-RA_si-CTCAGCGGGA_lane-002-chunk-001.fastq.gz | paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > 1126819/HumanPanc_1284133B_ROI1_FFPE_s4_A_S1_L002_R1_001.fastq) | cut -f 5-8 | tr "\t" "\n" > 1126819/HumanPanc_1284133B_ROI1_FFPE_s4_A_S1_L002_R2_001.fastq
# zcat 1126821/read-RA_si-TAAACCCTAG_lane-001-chunk-001.fastq.gz | paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > 1126821/HumanPanc_S001890_FFPE_s5_C_S2_L001_R1_001.fastq) | cut -f 5-8 | tr "\t" "\n" > 1126821/HumanPanc_S001890_FFPE_s5_C_S2_L001_R2_001.fastq
# zcat 1126821/read-RA_si-TAAACCCTAG_lane-002-chunk-001.fastq.gz | paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > 1126821/HumanPanc_S001890_FFPE_s5_C_S2_L002_R1_001.fastq) | cut -f 5-8 | tr "\t" "\n" > 1126821/HumanPanc_S001890_FFPE_s5_C_S2_L002_R2_001.fastq
# zcat 1126822/read-RA_si-GTCTAGTTTG_lane-001-chunk-001.fastq.gz | paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > 1126822/HumanPanc_S001891_FFPE_s3_D_S3_L001_R1_001.fastq) | cut -f 5-8 | tr "\t" "\n" > 1126822/HumanPanc_S001891_FFPE_s3_D_S3_L001_R2_001.fastq
# zcat 1126822/read-RA_si-GTCTAGTTTG_lane-002-chunk-001.fastq.gz | paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > 1126822/HumanPanc_S001891_FFPE_s3_D_S3_L002_R1_001.fastq) | cut -f 5-8 | tr "\t" "\n" > 1126822/HumanPanc_S001891_FFPE_s3_D_S3_L002_R2_001.fastq
# gzip */*.fastq


## ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_s3.R")
source("code/R/Utils/utils_pathology_spots.R")


## PARAMETERS ----


# The work folder

# S3 path information
s3_bucket <- "s3://XXXX"

## NGS360 project list
project_ids <- c("P-20200107-0002/200107_NB552115_0021_AHVJFYBGXC/", 
                 "P-20200617-0002/200616_NB552115_0038_AHVH7FBGXC/", 
                 "P-20200625-0003/200625_NB552115_0040_AHJ3NTBGXF/",
                 "P-20201203-0009/201203_NB501257_0337_AHYCYCAFXY/", 
                 "P-20201203-0009/201203_NB552115_0050_AH3GVYAFX2/", 
                 "P-20201203-0009/201221_NB501257_0343_AH5HCNAFX2/", 
                 "P-20201203-0009/201221_NB552115_0056_AH3FFHAFX2/", 
                 "P-20200612-0001/200612_NB552115_0037_AH2KKHBGXF/", 
                 "P-20200625-0002/200624_NB552115_0039_AHJ2VKBGXF/",
                 "P-20201216-0002/201215_NB501257_0340_AHC3HCBGXH/", 
                 "P-20201216-0002/201215_NB552115_0053_AHC57HBGXH/", 
                 "P-20201216-0002/201218_NB552115_0055_AHC3CLBGXH/", 
                 "P-20200928-0001/200928_NB552115_0043_AHFML3BGXF/", 
                 "P-20201104-0003/201113_A00267_0272_BHVML7DMXX/", 
                 "P-20201104-0003/201216_NB501257_0341_AHC3CGBGXH/", 
                 "P-20201104-0003/201216_NB552115_0054_AHC3TVBGXH/", 
                 "P-20201209-0001/201211_NB501257_0338_AHCHCGBGXH/",
                 "P-20201209-0001/201211_NB552115_0051_AHCHC2BGXH/",
                 "P-20201209-0001/201214_NB501257_0339_AHCJW3BGXH/",
                 "P-20201209-0001/201214_NB552115_0052_AHCJM7BGXH/",
                 "P-20201209-0001/201218_NB501257_0342_AH3CWTAFX2/",
                 "P-20201209-0001/201218_NB552115_0055_AHC3CLBGXH/",
                 "P-20201209-0001/201221_NB501257_0343_AH5HCNAFX2/",
                 "P-20201209-0001/201221_NB552115_0056_AH3FFHAFX2/"
                 )

## NGS360 project list
project_image_ids <- c("P-20200107-0002/", "P-20200617-0002/", "P-20200625-0003/",
                       "P-20201203-0009/", "P-20200612-0001/", "P-20200625-0002/",
                       "P-20201216-0002/", "P-20200928-0001/", "P-20201104-0003/", "P-20201209-0001/")

# Extensions to grab from the S3 bucket
extensions <- c("fastq.gz")

# Specific files to grab from Stash
pub_files <- c("filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5",
               "spatial/tissue_hires_image.png", "spatial/tissue_positions_list.csv", "spatial/scalefactors_json.json")

# File with meta data on Stash to rename FASTQ files for the publication
meta_file <- "XXXX"

# Stash path information
stash_folder <- "XXXX"

stash_folder1 <- paste0(stash_folder, "Split/")
dir.create(stash_folder1, showWarnings=FALSE)

stash_folder2 <- paste0(stash_folder, "MergedRenamed/")
dir.create(stash_folder2, showWarnings=FALSE)

stash_folder3 <- paste0(stash_folder, "Processed/")
dir.create(stash_folder3, showWarnings=FALSE)

stash_folder4 <- paste0(stash_folder, "Images/")
dir.create(stash_folder4, showWarnings=FALSE)


## INGEST META DATA ----


# Read sample file, filter for samples that pass QC
meta_df <- readr::read_delim(file=meta_file, delim="\t")


## S3 TO STASH ----


# Investigate the bucket before copying it
file_list_s3 <- c()
for(p in project_ids){
  file_list_s3 <- c(file_list_s3, ping_bucket_my(s3_bucket, s3_folder=p, extensions=extensions))
}

print(file_list_s3)

# For each sample, write all FASTQ files from S3 to Stash
# Parallelize the process or it takes too long; needs good Internet and Stash connectivity too

library(foreach)

cl <- length(project_ids)  # parallel::makeCluster(parallel::detectCores()) / 2
doParallel::registerDoParallel(cl)

foreach::foreach(p = project_ids) %dopar% {
  write_from_s3_to_stash_my(p, s3_bucket=s3_bucket, stash_path=paste0(stash_folder1, gsub("/", "_", p)), extensions=extensions)
};

parallel::stopCluster(cl)


## MERGE FASTQ LANES ----


# Find sets of FASTQ files to merge
file_list <- dir(stash_folder1, pattern="fastq.gz", full.names=FALSE, recursive=FALSE)

length(file_list_s3) == length(file_list)

group_list <- unique(gsub("_S\\d+_L00\\d_", "\\.\\*_", file_list))
group_list <- unique(gsub(paste0(paste(gsub("/", "_", project_ids), collapse="|")), "", group_list))

length(group_list) / 4

# Merge FASTQ files from different lanes but same sample
file_counter <- 0
for(g in group_list){
  sample_id <- gsub("\\.\\*_[IR][12]_001.fastq.gz", "", g)
  sample_suffix <- gsub(".*(_[IR][12]_001.fastq.gz)", "\\1", g)
  
  indx <- file_list[which(grepl(g, file_list))]
  
  meta_sample_df <- meta_df %>%
    dplyr::filter(grepl(sample_id, FullPath))
  
  if(nrow(meta_sample_df) == 1){
    command_string <- paste("cat", paste0(paste0(stash_folder1, indx), collapse=" "), ">", paste0(stash_folder2, meta_sample_df[1, "Sample_Name_Pub"], sample_suffix))
    
    cat(paste0(stash_folder2, meta_sample_df[1, "Sample_Name_Pub"], sample_suffix), " - ", length(indx), "\n")
    file_counter <- file_counter + length(indx)
    #print(command_string)
    
    system(command_string)
  } else{
    #cat("Struggling to find meta data for sample", g, " - ", sample_id, " - ", nrow(meta_sample_df), "\n")
  }
}

print(file_counter)
print(length(file_list))


## COPY AND RENAME SPACERANGER FILES ----


# Find all unique sample IDs from file names
sample_id_list <- gsub("\\.\\*_[IR][12]_001.fastq.gz", "", group_list)

# Find data corresponding to the sample ID on Stash
# Copy-over and rename specific files needed for the publication
for(s in sample_id_list){
  meta_sample_df <- meta_df %>%
    dplyr::filter(grepl(s, FullPath))
  
  if(nrow(meta_sample_df) == 1){
    for(p in pub_files){
      file_to_copy <- paste0(meta_sample_df[1, "FullPath"], p)
      command_string <- paste("cp", file_to_copy, paste0(stash_folder3, meta_sample_df[1, "Sample_Name_Pub"], "_", gsub("spatial/", "", p)))
      
      if(file.exists(file_to_copy)){
        # print(command_string)
        system(command_string)
      } else{
        cat("Struggling to find file to copy", file_to_copy, "\n")
      }
    }
  } else{
    cat("Struggling to find meta data for sample", s, "\n")
  }
}


## COPY AND RENAME IMAGES ----


# Investigate the bucket before copying it
file_list_s3 <- c()
for(p in project_image_ids){
  file_list_s3 <- c(file_list_s3, ping_bucket_my(s3_bucket, s3_folder=paste0(p, "TIFFs"), extensions=c(".tif")))
}
