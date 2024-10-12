# Author: Anna Lyubetskaya. Date: 21-01-18

# Update meta data information from file in a set of Seurat RDS objects

# Gather a list of Seurat RDS objects, update meta data from file, output to RDS
# Assist script for updating and augmenting data


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`


## PARAMETERS ----


# Sample / Cohort name
cohort_name <- "PDAC"
# List of data_seurat@misc[["user.X"]] to update
parameter_list <- c("pt.size.factor")


## PATHS ----


# Input file with sample meta data
meta_file <- "XXXX"

# Input folder
input_path <- "XXXX"


## INGEST DATA ----


# Read sample file; select only cohort-relevant samples
meta_df <- readr::read_delim(file=meta_file, delim="\t") %>%
  # dplyr::filter(Sample_Name %in% cohort_name)
  dplyr::filter(grepl(cohort_name, Sample_Name))


## WRANGLE DATA ----


for(i in 1:nrow(meta_df)){
  
  # Sample name
  sample_name <- meta_df[[i, "Sample_Name"]]

  # Find the RDS object
  file_list <- dir(input_path, pattern=sample_name, full.names=TRUE, recursive=TRUE)
  
  if(length(file_list) == 1){
    
    print(file_list[1])
    
    # Open a connection to the RDS object
    con <- gzfile(file_list[1])
    
    # Ingest the Seurat object
    data_seurat <- readRDS(con)
    
    # Close the connection to be able to overwrite
    close(con)

    for(field in parameter_list){
      data_seurat@misc[[paste0("user.", field)]] <- meta_df[[i, field]]
    }
    
    # Write the updated Seurat object
    filename <- gsub(".rds", "_upd.rds", file_list[1])
    saveRDS(data_seurat, file=filename)

    # Remove the old Seurat object
    unlink(file_list[1])
    # Rename the new Seurat object to have the name of the old Seurat object
    file.rename(filename, file_list[1])

  } else{
    warning(paste0("More than one file found!\n", file_list))
  }
}
