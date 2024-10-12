# Author: Anna Lyubetskaya. Date: 20-02-21
# Sync 10X data from NGS360 to Stash based on a factor sheet

# Relevant meta data columns:
# (Required) Sample_ID = name of the sample (SpaceRanger folder) exactly as it appears in NGS360
# (Required) Project_ID = NGS360 Project ID for the sample
# (Optional) Sample_Name = a new name for the sample (folder) that user wishes to assign; otherwise, will be the same as Sample_ID
# (Output) Pipeline.FullPath = this column, absent or present in the file initially, will be filled out with the location of the newly synced data

# Data will be deposited into the folder defined by the user in a standard structure


## ENVIRONMENT ----

print('--- Running 10X_ngs_s3_to_stash.R')

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_s3.R")


## PARAMETERS ----


# Only files with these extensions will be synced
extensions <- c("csv", "h5", "html", "jpg", "png", "json", "tsv", "csv", "res")

# File to use to check the integrity of the data
# For example, for a standard 10x pipeline: web_summary.html or cloupe.cloupe
integrity_check <- "cloupe.cloupe"


## PATHS ----


# S3 path information
s3_bucket <- "s3://XXXX"

# Factor sheet based on which to gather data
default_meta_file <- "XXXX"
# Place to write all the data on Stash
default_output_path <- "XXXX"


# Parse the arguments passed in from the cli (sets to defaults if none provided)
option_list <- list(
  optparse::make_option(c("-m", "--meta_file"), action="store", 
                        default=default_meta_file, help="Parse the path to the meta file/factor sheet."),
  optparse::make_option(c("-o", "--output_path"), action="store", 
                        default=default_output_path, help="Parse the path to the output location."))

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list)) 

# Meta data file listing and annotating all samples for processing
if (opt$meta_file == 'default'){
  meta_file <- default_meta_file
} else {
  meta_file <- opt$meta_file    
}


# Check if provided meta file exists
if (!file.exists(meta_file)){
  stop(paste0("ERROR: The following file does not exists: ", meta_file))
}

# Output location for Seurat objects
output_path <- opt$output_path
if (!endsWith(output_path, '/')){
  output_path <- paste0(output_path, '/')
}

# Remove double dashes from the output path if they happened
output_path <- gsub("//", "/", output_path)
output_path <- gsub("/+$", "", output_path)


# # Delete output folder if already exists
# if(dir.exists(output_path)){
#   unlink(output_path, recursive=TRUE)
# }

# Create output root folder if necessary
if(!file.exists(output_path)){
  dir.create(output_path, recursive=TRUE)
}

# Check if user has write permissions to output folder
if (file.access(output_path, mode = 2) != 0){
  stop(paste0("ERROR: Airflow pipeline does have write permissions to the specified output path where data is to be copied to from s3: ", output_path))
}

# Create subfolder for data
output_path <- paste0(output_path, "/1_SpaceRanger_outputs/")
dir.create(output_path, recursive=TRUE)

# Initiate a log file
log_file <- paste0(output_path, "/log.txt")
write("Files copied", log_file)


## READ META DATA ----


# Read meta data file
meta_df <- readr::read_delim(meta_file, delim="\t")

if(!"Sample_ID" %in% colnames(meta_df) || !"Project_ID" %in% colnames(meta_df)){
  stop(paste0("ERROR: Please check that your factor sheet has both Sample_ID and Project_ID columns ", meta_file))
}


# List of projects to sync
project_list <- unique(meta_df$Project_ID)


## INVESTIGATE NGS360 PROJECTS AND SYNC ----


# Iterate over requested projects
sample_ngs360_pooled <- list()
for(project in project_list){
  
  
  ## Investigate NGS360 project
  
  
  # Find all web_summaries in the project - expected of every successful 10x run
  file_list <- ping_bucket_my(s3_bucket, s3_folder=project, extensions=c(integrity_check))
  
  # Find the list of samples requested from this project
  sample_list <- meta_df %>%
    dplyr::filter(Project_ID == project) %>%
    dplyr::pull(Sample_ID)
  
  # Extract sample names from NGS360
  sample_ngs360_pooled[[project]] <- tibble::tibble(Sample_ID = gsub(".*/", "", gsub(paste0("/outs/", integrity_check), "", file_list)),
                                                    NGS360.FullPath = gsub(paste0("/outs/.*"), "/outs/", file_list))
  
  # Count unique sample folders identified on NGS360
  sample_ngs360_list <- sample_ngs360_pooled[[project]]$Sample_ID
  sample_ngs360_count <- table(sample_ngs360_list)
  
  # Check that files requested are present in the correct NGS360 project
  sample_check <- sapply(sample_list, function(x) x %in% sample_ngs360_list)
  
  print(project)
  print(sample_check)
  
  if(FALSE %in% sample_check){
    stop("ERROR: One of the samples requested for this project was not found!")
    print(sample_ngs360_count)
  }
  
  if(max(sample_ngs360_count) > 1){
    stop("ERROR: Some samples have multiple folders named exactly the same!")
    print(sample_ngs360_count)
  }
  
}


## Sync NGS360 project


# Add the NGS360 path to meta data
meta_df <- meta_df %>%
  dplyr::inner_join(dplyr::bind_rows(sample_ngs360_pooled), by="Sample_ID")

# Populate Sample_Name column from Sample_ID column if not defined by user
if(!"Sample_Name" %in% colnames(meta_df)){
  meta_df[["Sample_Name"]] <- meta_df[["Sample_ID"]]
}

# For each sample, write all but BAM/BAI files from S3 to Stash
for(i in 1:nrow(meta_df)){
  stash_path <- paste0(output_path, "/", meta_df$Project_ID[[i]], "_", meta_df$Sample_Name[[i]], "/")
  bucket_folder_list <- write_from_s3_to_stash_my(meta_df$NGS360.FullPath[[i]], 
                                                  s3_bucket=s3_bucket, 
                                                  stash_path=stash_path,
                                                  extensions=extensions)
  
  meta_df[[i, "Pipeline.FullPath"]] <- stash_path
  write(bucket_folder_list, log_file, append=TRUE)
}

# Write the updated factor sheet to the output folder
readr::write_delim(meta_df, paste0(output_path, "/meta_file.txt"), delim="\t")
