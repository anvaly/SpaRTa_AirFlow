# Author: Anna Lyubetskaya. Date: 19-09-20
# These are functions that assist with data inputs and outputs

# Limit the number of ingested lines from a file. E.g. for test runs.
n_max = 100


report_col_vals_df_my <- function(df, n_vals=10){
  ## List all columns and their unique values
  
  col_vals_unique <- lapply(df, function(x) unique(x))
  
  col_vals_unique_part <- col_vals_unique
  col_vals_unique_part[lengths(col_vals_unique) >= n_vals] <- NULL
  
  return(col_vals_unique)
}


report_col_num_my <- function(df){
  ## Report number of columns in the tibble
  
  # Report loaded columns
  cols <- colnames(df)
  col_num <- length(cols)
  # Adjust number of columns to print
  if(col_num > 20){
    col_num <- 20
  }
  
  cat("Columns, total =", length(cols), ":", cols[1:col_num], "\n")
}


read_file2df_my <- function(in_file, delim="\t", in_cols=NULL, verbose=FALSE, skip=0){
  ## Read a file into a data frame

  df <- readr::read_delim(file=in_file, delim=delim, col_types=in_cols, skip=skip)  # n_max=n_max
  
  # Report on loaded tibble
  if(verbose == TRUE){
    cat("Reading file =", in_file, "\n")
    report_col_num_my(df)
  }

  return(df)
}


read_dir2file_my <- function(in_path, in_regex=NULL, skip_row=0, in_cols=NULL, delim="\t", do_search=FALSE, verbose=FALSE){
  ## Read all files in the folder into a single data frame
  ## https://stackoverflow.com/questions/46299777/add-filename-column-to-table-as-multiple-files-are-read-and-bound
  
  # If provided with a string, assume it's a folder path to scan
  if(do_search == TRUE){
    file_list <- dir(in_path, pattern=in_regex, full.names=TRUE, recursive=TRUE)
    # If provided with a list, assume it's a list of file paths to read
  } else{
    file_list <- in_path
  }

  df <- tibble::tibble(File = file_list) %>%
    # tidyr::extract(File, "VRUNID", "/([A-Z]{2}[0-9]{6})-[A-Z]{2}[0-9]{6}", remove = FALSE) %>%
    dplyr::mutate(Data = lapply(File, function(x) 
      readr::read_delim(file=x, delim=delim, skip=skip_row, col_types=in_cols))) %>%
    tidyr::unnest(Data)
  
  # Report on loaded tibble
  if(verbose == TRUE){
    cat('Files found =', length(file_list), "\n")
    
    report_col_num_my(df)
    
    #Report rows per file
    print("Files and number of rows uploaded per file")
    print(df %>% 
            dplyr::group_by(File) %>%
            dplyr::tally())
  }
  
  return(df) #%>% dplyr::select(-File)
}


create_output_subfolders_my <- function(global_out, folder_labels){
  ## In a given location, create a set of folders if they don't exist already
  
  dir.create(global_out, showWarnings = FALSE)
  
  output_folders <- list()
  for(f in folder_labels){
    output_folders[[f]] <- paste0(global_out, f, "/")
    dir.create(output_folders[[f]], showWarnings = FALSE)
  }
  
  return(output_folders)
}


read_h5_my <- function(input_file){
  ## Read a file in H5 format
  ## rhdf5 library for H5 manual: https://www.bioconductor.org/packages/release/bioc/vignettes/rhdf5/inst/doc/rhdf5.html
  
  # Look at the data structure
  h5_structure <- rhdf5::h5ls(input_file)
  # Open the H5 file connection
  h5_data <- rhdf5::h5read(input_file, "/")
  
  return(h5_data)
}
