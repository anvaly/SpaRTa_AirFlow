# Author: Anna Lyubetskaya. Date: 20-02-04

# Your BMS AWS credentials should be saved as user variables: AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY
# https://domino.web.bms.com/account#user-variables
# Check your environment is setup correctly: Sys.getenv("AWS_ACCESS_KEY_ID"), Sys.getenv("AWS_SECRET_ACCESS_KEY")

# Utilities to connect to S3
# aws.s3: https://github.com/cloudyr/aws.s3

# Important BMS internal NGS buckets:
# - "s3://XXXX" for BAMs and counts
# - "s3://XXXX" for FASTQs


analyze_s3_content_my <- function(bucket_content, extensions=NULL){
  ## Break up bucket content into a list of file names, their sizes, and extensions

  # List of all file names in a given location of an S3 bucket
  bucket_folder_list <- as.vector(unlist(lapply(bucket_content, function(x) x$Key)))
  # List of all file sizes in a given location of an S3 bucket
  bucket_size_list <- as.vector(unlist(lapply(bucket_content, function(x) x$Size)))
  
  # Select files with select extensions
  if(!is.null(extensions)){
    index <- grepl(paste(paste0(".", extensions, "$"), collapse="|"), bucket_folder_list)
    
    # Subset file names and sizes
    bucket_folder_list <- bucket_folder_list[index]
    bucket_size_list <- bucket_size_list[index]
  }
  
  # Total size of all files in a given location of an S3 bucket
  bucket_size_gb <- sum(bucket_size_list) / 1e9
  # A list of all file extensions detected
  bucket_extensions <- unique(gsub("^.+\\.", "", bucket_folder_list))
  
  cat("File number =", length(bucket_folder_list), "\n", 
      "File size (GB) =", bucket_size_gb, "\n", 
      "File extensions =", bucket_extensions, "\n")
  
  return(bucket_folder_list)
}


write_from_s3_to_stash_my <- function(s3_folder, s3_bucket="s3://XXXX", stash_path=NULL,
                                      extensions=NULL){
  ## Find a file on S3 and write it to Stash using aws.s3
  ## Important: please, skip the following extensions to presever Stash space: "gz", "bam", "bai"
  
  ## Examples: 
  ## extensions <- c("csv", "h5", "html", "jpg", "png", "json", "tsv", "csv", "res")
  ## extensions <- c("filtered_feature_bc_matrix.h5", "metrics_summary.csv", "clusters.csv", "differential_expression.csv")
  
  # Find all objects in a given location of an S3 bucket
  bucket_content <- aws.s3::get_bucket(bucket=s3_bucket, prefix=s3_folder)
  
  # List of files with select extensions in a given location of an S3 bucket
  bucket_folder_list <- analyze_s3_content_my(bucket_content, extensions=extensions)

  # Loop through files and save them to a Stash location
  for(file_name in bucket_folder_list){
    s3_file_path <- paste0(s3_bucket, file_name)
    check_object <- aws.s3::head_object(s3_file_path)
    
    if(check_object == TRUE && !is.null(stash_path)){
      # s3_file_object <- aws.s3::get_object(s3_file_path)
      
      stash_file_path <- paste0(stash_path, gsub(paste0(s3_bucket, s3_folder), "", s3_file_path))
      aws.s3::save_object(s3_file_path, file=stash_file_path)
    }
  }
  
  return(bucket_folder_list)
}


ping_bucket_my <- function(s3_bucket="s3://XXXX", extensions=NULL){
  ## Test that you can reach a BMS internal bucket, e.g. s3://XXXX
  
  cat("Bucket =", s3_bucket, " ", "Folder =", s3_folder, "\n")
  bucket_content <- aws.s3::get_bucket(bucket=s3_bucket, prefix=s3_folder,
                                       max=10000000)
  
  file_list <- analyze_s3_content_my(bucket_content, extensions=extensions)
  
  return(file_list)
}
