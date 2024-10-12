# Author: Anna Lyubetskaya. Date: 19-11-11
# Functions associated with reading and standardizing ARCHS4 dataset

# ARCHS4: https://amp.pharm.mssm.edu/archs4/download.html
# ARCHS4 data is available on gene and transcript level.
# Gene level data is available as RDA while transcript level data is available as H5.
# The structure of gene RDA data and transcript H5 data is not the same.
# Reading transcript level ARCHS4 data requires a large AWS instance while reading gene level ARCHS4 data requires a medium AWS instance.


#source("code/utils/utils_archs4.R")
source("code/utils/utils_in_out.R")


setup_read_archs4 <- function(cell_line_list=c("MC38", "B16"), group="all"){
  ## Setup and read ARCHS4 data

  # Location of ARCHS4 data (end-of-2019 version)  
  archs4_path <- "XXXX"
  archs4_path_der <- paste0(archs4_path, "Input_ARCHS4_derivative/")

  # File translating cell line names to GEO IDs
  cell_line_file <- paste0(archs4_path, "line2ID_", group, ".txt")

  archs4_files <- list()
  data_types <- c("expression", "samples", "genes")
  # A relevant slice of ARCHS4 expression, sample, and gene data data is stored here: Generated on the first read of full data
  for(d in data_types){
    archs4_files[[d]] <- paste0(archs4_path_der, "mouse_select_", group, "_", d, ".txt")
  }

  # Load a list of cell lines and their GSM IDs
  line2ID_df <- read_file2df_my(cell_line_file) %>% 
    tidyr::separate_rows(GEO_ID, sep=",") %>%
    dplyr::filter(cell_line %in% cell_line_list)
  
  # Read only relevant ARCHS4 data: expression data, samples, and genes
  data_archs4 <- read_processed_archs4_my(NULL, line2ID_df$GEO_ID, archs4_files)
  
  # Update sample meta data
  data_archs4$samples <- data_archs4$samples %>%
    dplyr::inner_join(line2ID_df, by="GEO_ID")
  
  
  return(data_archs4)
}


read_processed_archs4_my <- function(file_data_archs4, geo_list, output_files){
  ## Load full ARCHS4 data if necessary, subset to a set of cell lines and save result to files
  ## Read 3 files of ARCHS4 data corresponding to expression, sample, and gene data of a select set of cell lines
  
  # Load and slice full ARCHS4 data
  # Only runs if no aready processed data is available
  if(!is.null(file_data_archs4)){
    process_raw_archs4_my(file_data_archs4, geo_list, output_files)
  } else{
    archs4_files <- output_files
  }
  
  data <- list()
  
  # Load a slice of ARCHS4 expression data selected on the first load of full data; see the else statement
  data[["expression"]] <- readr::read_delim(archs4_files$expression, delim="\t") %>%
    dplyr::mutate_if(is.numeric, as.integer) %>%
    dplyr::mutate(Symbol = toupper(Symbol)) %>%
    dplyr::arrange(Symbol)
  
  data[["expression"]] <- data[["expression"]][, c("Symbol", intersect(geo_list, colnames(data[["expression"]])))]
  
  # Load a slice of ARCHS4 meta data selected on the first load of full data; see the else statement
  data[["samples"]] <- read_file2df_my(archs4_files$samples, delim="\t")
  
  # Load a slice of ARCHS4 meta data selected on the first load of full data; see the else statement
  data[["genes"]] <- read_file2df_my(archs4_files$genes, delim="\t") %>%
    dplyr::rename("Symbol" = "value") %>%
    dplyr::mutate(Symbol = toupper(Symbol))
  
  return(data)
}


process_raw_archs4_my <- function(file_data_archs4=NULL, geo_list=NULL, output_files=NULL, group="all"){
  ## If no processed data found, upload all ARCHS4 data and put relevant slices into files
  
  # archs4_path <- "XXXX"
  # archs4_path_der <- paste0(archs4_path, "Input_ARCHS4_derivative/")
  
  # Find ARCHS4 RDA file
  if(is.null(file_data_archs4)){
    file_data_archs4 <- paste0(archs4_path, "mouse_matrix_v11.rda")
  }
  
  # Find a list of GEO IDs of interest
  if(is.null(geo_list)){
    geo_list <- read_mouse_cell_line_list_my(group=group)$GEO_ID
  }
  
  if(is.null(output_files)){
    data_types <- c("expression", "samples", "genes")
    output_files <- list()
    
    # A relevant slice of ARCHS4 expression, sample, and gene data data is stored here: Generated on the first read of full data
    for(d in data_types){
      output_files[[d]] <- paste0(archs4_path_der, "mouse_select_", group, "_", d, ".txt")
    }
  }
  
  if(!file.exists(output_files$expression)){
    dir.create(archs4_path_der, showWarnings = FALSE)
    
    # Load full ARCHS4 dataset
    if(file.exists(file_data_archs4)){
      read_slice_archs4_rda_my(file_data_archs4, geo_list, output_files)
    } else{  # Expect H5 format
      read_slice_archs4_h5_my(gsub(".rda", ".h5", file_data_archs4), geo_list, output_files)
    }
  }
  
}


read_mouse_cell_line_list_my <- function(cell_line_file=NULL, group="all"){
  ## Read a list of cell lines of interest
  
  # Read the file translating cell line names to GEO IDs  
  if(is.null(cell_line_file)){
    cell_line_file <- paste0(archs4_path, "line2ID_", group, ".txt")
  }
  
  # Load a list of cell lines and their GSM IDs
  line2ID_df <- read_file2df_my(cell_line_file) %>% 
    tidyr::separate_rows(GEO_ID, sep=",")
  
  return(line2ID_df)
}


read_slice_archs4_rda_my <- function(file_data_archs4, geo_list, output_files){
  ## Load and slice full ARCHS4 dataset from RDA format
  ## ARCHS4 dataset: "expression" matrix and "meta" named list
  
  ## Important: Reading gene level ARCHS4 data requires a medium AWS instance
  
  load(file = file_data_archs4)
  
  # Select relevant meta data
  expression_select <- subset(expression, select = geo_list)
  write.table(expression_select, file = output_files$expression, sep = "\t", col.names = NA, row.names = TRUE)
  
  # Find relevant GEO IDs in meta data
  sample_index_select <- meta$Sample_geo_accession %in% geo_list
  
  # Select important meta fields from 37 ARCHS4 meta fields
  # Other interesting fields: "Sample_extract_protocol_ch1", "Sample_data_processing"
  meta_fields <- c("Sample_geo_accession", "Sample_title", "Sample_source_name_ch1", "reads_aligned", "total_reads")  
  
  # Select sample data by meta field and by GEO accession
  sample_meta_data_select <- list()
  for(field in meta_fields){  # Each of these fields should be 209395 long
    sample_meta_data_select[[field]] <- meta[[field]][sample_index_select]
  }
  
  # Select gene meta data
  gene_meta_data_select <- list()
  for(field in names( meta[lengths(meta) == 32544] )){
    gene_meta_data_select[[field]] <- meta[[field]]
  }
  
  # Write selected sample meta data to file
  readr::write_delim(tibble::as_tibble(sample_meta_data_select) %>%
                       dplyr::rename("GEO_ID" = "Sample_geo_accession"), 
                     path = output_files$samples, delim = "\t")
  
  # Write selected sample gene data to file
  readr::write_delim(tibble::as_tibble(gene_meta_data_select), path = output_files$genes, delim = "\t")
}


read_slice_archs4_h5_my <- function(file_data_archs4, geo_list, output_files){
  ## Load and slice full ARCHS4 dataset from H5 format
  ## ARCHS4 dataset: "expression" matrix and "meta" named list
  ## rhdf5 library for H5 manual: https://www.bioconductor.org/packages/release/bioc/vignettes/rhdf5/inst/doc/rhdf5.html
  ## Look at the data structure rhdf5::h5ls(file_data_archs4)
  
  ## Important: Reading transcript level ARCHS4 data requires a large AWS instance.
  
  # Open the H5 file connection
  h5f <- rhdf5::H5Fopen(file_data_archs4)
  # Read expression structure
  expression <- h5f$data
  # Read sample structure
  meta <- h5f$meta
  # Break the connection to the H5 file
  rhdf5::h5closeAll()
  
  # Find all sample IDs
  sample_list <- meta$sample$geo_accession
  # Find all gene symbols
  gene_list <- meta$genes$gene_symbol
  
  # Find relevant GEO IDs in meta data
  sample_index_select <- which(meta$sample$geo_accession %in% geo_list & grepl("RNA", meta$sample$molecule_ch1))
  
  dim(expression$expression)
  print(c(length(gene_list), length(sample_list)))
  
  # Select relevant data
  expression_select <- t(expression$expression[sample_index_select, ])
  dim(expression_select)
  
  # Add column and row names
  colnames(expression_select) <- sample_list[sample_index_select]
  rownames(expression_select) <- gene_list
  
  # Write expression data to file
  readr::write_delim(tibble::as_tibble(expression_select, rownames="Symbol"), 
                     file=output_files$expression, delim="\t")
  
  # Select important meta fields from 37 ARCHS4 meta fields
  meta_fields <- c("geo_accession", "title", "readsaligned", "source_name_ch1", "molecule_ch1", "characteristics_ch1")
  
  # Select sample data by meta field and by GEO accession
  sample_meta_data_select <- list()
  for(field in meta_fields){
    sample_meta_data_select[[field]] <- meta$samples[[field]][sample_index_select]
  }
  
  # Meta data tibble
  meta_df <- tibble::as_tibble(sample_meta_data_select) %>%
    dplyr::rename("GEO_ID" = "geo_accession")
  
  # Write selected sample meta data to file
  readr::write_delim(meta_df, path = output_files$samples, delim = "\t")
  
  # Write selected sample gene data to file
  readr::write_delim(tibble::as_tibble(meta$genes$gene_symbol), file=output_files$genes, delim="\t")
}
