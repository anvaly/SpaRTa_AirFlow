# Author: Anna Lyubetskaya. Date: 23-08-29
# Annotate each barcode with an integrated pathology class


## ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`


## PARAMETERS ----


## PATHS ----


# Cluster labels assigned to each cluster across the cohorts
cohort_annotation_file <- "XXXX"

# The list of cohorts to process
cohort_rds_list <- paste0("XXXX",
                          c("PDAC108_path14_parenchyma_4K_reftop20_filt2K_cca_sct.rds",
                            "PDAC108_path14_stroma_4K_reftop20_filt2K_cca_sct.rds"))

# Cohort labels
cohort_labels <- c("parenchyma", "stroma")

# Output folder
output_file <- "XXXX"


## WRANGLE ----


# Read the cohort annotation file
cohort_annotation <- readr::read_delim(cohort_annotation_file, delim="\t") %>%
  dplyr::mutate(cluster = as.factor(cluster))

# Iterate over cohorts
barcode_list <- list()
for(i in 1:length(cohort_rds_list)){
  
  # Read in the Seurat RDS
  data_seurat <- readRDS(cohort_rds_list[i])
  
  # Cohort annotation
  annotation_df <- cohort_annotation %>%
    dplyr::filter(cohort == cohort_labels[i]) %>%
    dplyr::select(-cohort)
  
  # Name the cluster column with the resolution
  colnames(annotation_df) <- gsub("cluster", annotation_df[1, "resolution"], colnames(annotation_df))
  
  # Extract and annotate barcodes
  barcode_list[[i]] <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate") %>%
    dplyr::select(Coordinate, rlang::sym(annotation_df[[1, "resolution"]])) %>%
    dplyr::left_join(annotation_df %>%
                       dplyr::select(-resolution), by=annotation_df[[1, "resolution"]])
  
  # Cleanup
  rm(data_seurat)
  gc()
}


coordinate_list <- unique(dplyr::bind_rows(barcode_list)$Coordinate)

# Combine all barcode annotations
barcode_df <- dplyr::bind_rows(barcode_list) %>%
  dplyr::group_by(Coordinate) %>%
  dplyr::summarise(n = dplyr::n_distinct(integrated_snn_res.0.3),
                   n_call = dplyr::n_distinct(call),
                   n_split = dplyr::n_distinct(split),
                   call = paste(sort(unique(call)), collapse=";"),
                   split = paste(sort(unique(split)), collapse=";"))

duplicates <- barcode_df %>% 
  dplyr::filter(n > 1)

duplicates_v2 <- barcode_df %>% 
  dplyr::filter(n > 1 & n_call > 1)

sort(table(duplicates_v2$call))
sort(table(duplicates_v2$split))

# Resolve ambiguous cohort assignments
barcode_df <- barcode_df %>%
  dplyr::mutate(label = gsub(";.+", "", split))

# Write the result to file
readr::write_delim(barcode_df, paste0(output_file, ".txt"), delim="\t")


# Print subcohorts for future integration
for(l in unique(barcode_df$label)){
  # Resolve ambiguous cohort assignments
  barcode_loc_df <- barcode_df %>%
    dplyr::filter(label == l)
  
  # Write the result to file
  readr::write_delim(barcode_loc_df, paste0(output_file, "_", l, ".txt"), delim="\t")
}
