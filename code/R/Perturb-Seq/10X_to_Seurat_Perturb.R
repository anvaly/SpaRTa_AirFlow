# Author: Anna Lyubetskaya. Date: 20-11-27
# Create a normalized, filtered Seurat object from Perturb-Seq 10X folder data


## ENVIRONMENT ----


# Allow piping throughout the package
`%>%` <- magrittr::`%>%`


source("code/utils/utils_tibble.R")
source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_in_out_perturb.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_10X_matrix.R")


## PARAMETERS ----


guide_umi_top_threshold <- 50
guide_umi_top_ratio <- 5

mt_threshold <- 20
ribo_threshold <- 30

cohort_regex <- "MC38_Glides"

## PATHS ----


# Input location containing 10X folders
input_path <- "XXXX"

# Output location for Seurat objects
output_path <- "XXXX"

# Create output folder
dir.create(output_path, showWarnings=FALSE)

# Input file with sample meta data
meta_file <- paste0(input_path, "PerturbSeq_Metadata_Sheet.txt")


## INGEST, WRANGLE, OUTPUT ----


# Read sample file
meta_df <- readr::read_delim(file=meta_file, delim="\t") %>%
  tidyr::drop_na() %>%
  dplyr::filter(grepl(cohort_regex, Cohort))


for(i in 1:nrow(meta_df)){
  
  # Meta data parameters
  filename_in <- meta_df[[i, "FullPath"]]
  filename_gex_in <- paste0(filename_in, "filtered_feature_bc_matrix/matrix_gex.mtx")
  filename_crispr_in <- paste0(filename_in, "filtered_feature_bc_matrix/matrix_crispr.mtx")
  sample_name <- meta_df[[i, "Sample_Name"]]
  
  filename_out1 <- paste0(output_path, sample_name, "_all.rds")
  filename_out2 <- paste0(output_path, sample_name, "_crispr.txt")
  filename_out3 <- paste0(output_path, sample_name, "_crispr_summary.txt")

  
  ## Split data ----
  
  # If it hasn't been done before, split the Perturb-Seq H5 object into GEX and CRISPR parts
  if(!file.exists(filename_gex_in)){
    read_split_10X_matrix_my(filename_in)
  }
  
  
  # if(!file.exists(filename_out1)){
    print(filename_in)
    
    
    ## Load data ----
    
    # Read a 10X Spatial folder
    data_seurat <- matrix_to_Seurat_my(filename_gex_in)
    
    # Read the CRISPR data from H5
    data_seurat_crispr <- matrix_to_Seurat_my(filename_crispr_in)

    
    ## Process the guide data ----
    
    # Guide data as a wide tibble
    crispr_wide_df <- tibble::as_tibble(as.matrix(Seurat::GetAssayData(data_seurat_crispr)), rownames="Guides")
                                   
    # Write guide data to a file
    readr::write_delim(crispr_wide_df, file=filename_out2, delim="\t", col_names=TRUE)

    # Long tibble of guide UMIs
    crispr_df <- crispr_wide_df %>%
      df_wide2long_my(key="Coordinate", val="Guide_UMI") 
    
    # Find well cells with well covered top guide
    coordinate_top_guide_df <- crispr_df %>% 
      dplyr::group_by(Coordinate) %>% 
      dplyr::arrange(desc(Guide_UMI)) %>% 
      dplyr::slice(1) %>%
      dplyr::select(Coordinate, Guides) %>%
      dplyr::rename(TopGuide = Guides)
      
    # Find top 2 guides for each coordinate
    crispr_stat_df <- crispr_df %>%
      dplyr::group_by(Coordinate) %>%
      dplyr::slice_max(n=2, Guide_UMI, with_ties=FALSE) %>%
      dplyr::arrange(Coordinate, desc(Guide_UMI))
    
    # Add a Guide Status column
    crispr_stat_df[["GuideStatus"]] <- rep(c("Guide1", "Guide2"), nrow(crispr_stat_df)/2)
    
    # Make the stat table wide and calculate the ratio between the 2 top guides
    crispr_stat_wide_df <- crispr_stat_df %>%
      df_long2wide_my(rows="Coordinate", cols="GuideStatus", value="Guide_UMI") %>%
      dplyr::inner_join(coordinate_top_guide_df, by="Coordinate") %>%
      dplyr::mutate(Ratio = round(Guide1 / (Guide2+1), 1),
                    Pass = Guide1 >= guide_umi_top_threshold & Ratio >= guide_umi_top_ratio,
                    Good = ifelse(Guide1 >= guide_umi_top_threshold & Ratio >= guide_umi_top_ratio, TopGuide, ""),
                    Background = ifelse(Guide1 == 0, "background", ""),
                    Category = paste0(Good, Background))
    
    # Summary of cells that had good guide data
    table(crispr_stat_wide_df$Pass)

    # Write the result of guide analysis to file
    readr::write_delim(crispr_stat_wide_df, filename_out3, delim="\t")
    
    # Ensure that the annotation vector order matches Seurat meta data rowname order
    crispr_stat_df <- crispr_stat_df[match(rownames(data_seurat@meta.data), crispr_stat_df$Coordinate),]
    # Check that meta data matches the annotation tibble
    table(rownames(data_seurat@meta.data) == crispr_stat_df$Coordinate)
    
    # Add Guide category to the meta data
    data_seurat@meta.data[["CRISPRGuideCategory"]] <- crispr_stat_wide_df$Category
    data_seurat@meta.data[["CRISPRGeneCategory"]] <- gsub("\\.\\d+$", "", crispr_stat_wide_df$Category)
    
    # List of coordinates with definitive guides
    coordinates_select <- crispr_stat_wide_df$Coordinate[which(crispr_stat_wide_df$Category != "")]
    # Remove all cells with ambigious guides
    # data_seurat <- subset(data_seurat, cells=coordinates_select)
    
    
    ## Calculate ribosomal and mitochondrial content ----
    
    # Calculate mitochondrial and ribosomal content
    data_seurat <- calculate_mt_ribo_my(data_seurat)
    
    # Select coordinates with good MT and ribo content
    coordinates_select <- rownames(data_seurat@meta.data %>% 
                                     dplyr::filter(mito_percent <= mt_threshold & ribo_percent <= ribo_threshold))
    # Remove all cells with too high MT and ribo content
    data_seurat <- subset(data_seurat, cells=coordinates_select)

    
    
    ## SCTransform ----
    
    # Perform SCTransform normalization
    data_seurat <- Seurat::SCTransform(data_seurat, assay="RNA",vst.flavor = "v2", vars.to.regress = c("mito_percent", "ribo_percent"),
                                       return.only.var.genes = FALSE, verbose = FALSE)
    
    
    ## Calculate mean QC metrics on final objects ----
    
    # List of meta data fields to parse
    field_list <- c("nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT", "mito_percent", "ribo_percent")
    
    # Add mean QC metrics to the misc slot
    for(f in field_list){
      val <- mean(data_seurat@meta.data[[f]])
      data_seurat <- add_misc_to_seurat_object_my(data_seurat, field=paste0("qc.mean.", f), dict_value=val)
    }
    # Add user-defined meta-data to a Seurat object misc slot tagged with "user." prefix
    for(v in colnames(meta_df)){
      data_seurat <- add_misc_to_seurat_object_my(data_seurat, field=paste0("user.", v), dict_value=meta_df[[i, v]])
    }
    
    ## Write data ----
    
    # Save the full dataset to RDS
    saveRDS(data_seurat, file = filename_out1)

  # }
}
