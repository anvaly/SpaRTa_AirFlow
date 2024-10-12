# Author: Anna Lyubetskaya. Date: 21-09-10

# Create a diet Seurat object


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

library(Seurat)


## PARAMETERS ----


# Seurat RDS files are tagged as follows
cohort_name <- "MSKCC_lung_chestwall"

# Split the original object into subsets
do_split_by <- "user.Sample_Name"  # "user.Sample_Name" or NULL

# Other user-defined columns to add to the diet output
user_cols <- c()  # "CellCounts"


## PATHS ----


# Input folder
input_path <- "XXXX"

# Output folder
output_path <- "XXXX"

# Create output folders
dir.create(output_path, showWarnings = FALSE)

# Output meta data column lists
col_meta_output <- paste0(output_path, "col_meta_list.txt")
# Start the meta data column list file
write(paste(c("Sample_Name", "Meta column list"), collapse="\t"), file=col_meta_output, append=TRUE)


## INGEST DATA ----


# Find all files following the regex
file_list <- dir(input_path, pattern=paste0(cohort_name, ".*.rds"), full.names=TRUE)


## WRANGLE DATA ----


diet_to_file <- function(data_seurat, filename_output){
  ## Diet and object and write to file
  
  Seurat::DefaultAssay(data_seurat) <- "SCT"
  
  # Strip away all unnecessary parts
  data_seurat_diet <- Seurat::DietSeurat(data_seurat, assays = "SCT", counts = FALSE)
  
  # Save diet Seurat object
  saveRDS(data_seurat_diet, file = filename_output)
  
  return(data_seurat_diet)
}


for(rds_file in file_list){
  
  # Load data
  data_seurat <- readRDS(rds_file)
  
  # Check if important meta data fields are present and if not, dummy them
  col_names <- colnames(data_seurat@meta.data)
  
  if(!("user.Clustering" %in% names(data_seurat@misc))){
    data_seurat@misc$user.Clustering <- "temp_snn_res"
    data_seurat@meta.data[["temp_snn_res"]] <- 0
  }
  
  if(!("Pathology.Group" %in% col_names)){
    data_seurat@meta.data[["Pathology.Group"]] <- "Tissue"
  }
  
  # Find pathology columns
  path_cols <- sort(col_names[grep("Pathology.*percent", col_names)])
  
  # Find clustering columns
  clust_cols <- sort(col_names[grep("_snn_res", col_names)])
  
  # The following columns will remain in the output
  col_names_retain <- c("nCount_Spatial", "nFeature_Spatial", "nCount_SCT", "nFeature_SCT", "Pathology.Group",
                        user_cols, path_cols, clust_cols)
  
  if(is.null(do_split_by)){
    
    # Ouptut path
    filename_output <- gsub(".rds", "_diet.rds", paste0(output_path, "/", gsub("^.+/", "", rds_file)))
    print(filename_output)
    
    # Remove extra meta data
    data_seurat@meta.data <- data_seurat@meta.data %>%
      dplyr::select(dplyr::all_of(col_names_retain))
    
    # Diet and write to file
    data_seurat_diet <- diet_to_file(data_seurat, filename_output)
    cat(paste0("File Created:", filename_output, "\n"))
    
    # Write the list of all columns to a file
    write(paste(c(data_seurat_diet@misc$user.Sample_Name, 
                  paste(colnames(data_seurat_diet@meta.data), collapse=";")), collapse="\t"),
          file=col_meta_output, append=TRUE)
    
  } else{
    
    # Split data by a factor before dieting
    split_factor <- unique(data_seurat@meta.data[[do_split_by]])
    
    # Cycle through every factor value and split the object by barcodes
    for(s in split_factor){
      # Find a subset of barcodes by selected meta feature
      barcode_list <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinates") %>%
        dplyr::filter(!!rlang::sym(do_split_by) == s) %>%
        dplyr::pull(Coordinates)
      
      # Output path
      filename_output <- paste0(output_path, "/", s, "_diet.rds")
      print(filename_output)
      
      # Diet and write to file
      data_seurat_subset <- subset(data_seurat, cells=barcode_list)
      
      # Remove dashes because R
      s <- gsub("-", ".", s)
      
      if(s %in% names(data_seurat_subset@images)){
        # Remove extra images
        data_seurat_subset@images <- data_seurat_subset@images[s]
        
        # Remove extra meta data
        data_seurat_subset@meta.data <- data_seurat_subset@meta.data %>%
          dplyr::select(dplyr::all_of(col_names_retain))
        
        # Diet and write to file
        data_seurat_diet <- diet_to_file(data_seurat_subset, filename_output)
        
        # Write the list of all columns to a file
        write(paste(c(s, paste(colnames(data_seurat_diet@meta.data), collapse=";")), collapse="\t"),
              file=col_meta_output, append=TRUE)
        
      } else{
        warning("Image ", s, " not found in ", rds_file, "\n")
      }
    }
    
  }
}
