# Author: Andrew Fisher. Date: 23-09-14
# Integrates pathology distance computations with merged seurat object


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
library(Seurat)
library(tidyverse)
library(httr)
httr::set_config(httr::config(http_version = 2))

source("code/R/Utils/utils_10X_signatures.R")
source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_matrix.R")

source("code/utils/utils_signatures.R")

# apply scaling factor to distances (px -> um)
fetch_px_res <- function(conc_id){
  conc_user <- Sys.getenv("concentriq_api_user")
  conc_pass <- Sys.getenv("concentriq_api_pass")
  
  base <- 'XXXX'
  endpoint <- 'images/'
  
  # First query Concentriq to get image id
  call1   <- paste0(base,endpoint,conc_id)
  get_img <- httr::GET(call1, authenticate(conc_user, conc_pass, type='basic'), config = httr::config(ssl_verifypeer = FALSE))
  api_return <- content(get_img)
  
  # pull metadata
  px_res <- api_return$data$mppx
  
  return(px_res)
}


## PARAMETERS ----

# Run name
run_name <- "PDAC_TxNaive_Trichrome"

# Sample / Cohort name
cohort_name <- "PDAC"

# Write the final RDS object
do_rds_output <- TRUE

# List of specific samples to select (or NULL)
sample_list <- NULL # c("PDAC_E1265_ROI1", "PDAC_E1265_ROI2", "PDAC_E1265_ROI3", "PDAC_E1265_ROI4")

# Gene abundance filters
sct_threshold <- 0.5
spot_threshold <- 5

# Signatures to select
sig_select <- c("PDAC.collisson.classical", "PDAC.moffitt.basal", "PDAC.moffitt.activatedstroma", "PDAC.moffitt.normalstroma",
                "PDAC.P19.Bcell", "PDAC.P19.Fibroblast", "PDAC.P19.Macrophage", "PDAC.P19.Tcell",
                "BMS.Pathway.TNFa", "BMS.Pathway.TGFB", "BMS.Pathway.IFNg", "BMS.Pathway.IFNa",
                "BMS.PDAC.Hypoxia.CL", "PDAC.U.Immune.Tcell.exhausted"
)

# Path to a signature file
sig_path <- "data/import/Signatures/signatures_220210.txt"


## PATHS ----

# Input file with sample meta data
meta_file <- "/XXXX/meta_data_pdac_probes.txt"

# Input folder
input_path <- "/XXXX/PDAC108_path14_merge.rds"

# Input folder for computed distances
distances_path <- "/XXXX/PDAC_TxNaive_PathDistances"
tissue_path <- "/XXXX/PDAC_TxNaive_TissueDistances"

# Output folder for RDS objects
output_base_path <- "/XXXX/h_Distance_models/"
output_path <- file.path(output_base_path, run_name)

# Output folder for images
output_figs <- file.path(output_path, 'figs')

# Create the output folder for all figures
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
dir.create(output_figs, recursive = TRUE, showWarnings = FALSE)


## INGEST DATA ----

# Read sample file; select only cohort-relevant samples
meta_df <- readr::read_delim(file=meta_file, delim="\t") %>%
  dplyr::filter(grepl(cohort_name, Sample_Name) &
                  Best_Looking == TRUE)

# Read distance files
dist_ls <- list.files(distances_path, pattern = "\\.csv$", full.names = TRUE)
distances_df_ls <- sapply(dist_ls, read.csv, USE.NAMES = TRUE)
names(distances_df_ls) <- basename(names(distances_df_ls))
tiss_ls <- list.files(tissue_path, pattern = "\\.csv$", full.names = TRUE)
tissue_df_ls <- sapply(tiss_ls, read.csv, USE.NAMES = TRUE, simplify = FALSE)
names(tissue_df_ls) <- basename(names(tissue_df_ls))

# merge into single df
distances_df <- dplyr::bind_rows(distances_df_ls) %>%
  dplyr::left_join(dplyr::bind_rows(tissue_df_ls), by = "X") %>%
  tibble::column_to_rownames("X")

# ensure no duplicate row.names in table
if(nrow(distances_df) > length(unique(row.names(distances_df)))){
  stop("Duplicate probe names found in distances table")   
}

# Filter meta data for samples of interest if requested
if(!is.null(sample_list)){
  meta_df <- meta_df %>%
    dplyr::filter(Sample_Name %in% sample_list)
}

# add value to metadata
meta_df$um_per_px <- sapply(meta_df$Concentriq_Image_ID, fetch_px_res)
# add value to distance dataframe matching sample ID
scaled_distances_df <- distances_df %>%
  rownames_to_column(var = "Barcode") %>%
  mutate(Sample_Name = gsub("(.*)\\:(.*)", "\\1", Barcode)) %>%
  left_join(dplyr::select(meta_df, c("Sample_Name", "um_per_px"))) %>%
  column_to_rownames("Barcode") %>%
  dplyr::select(-Sample_Name) %>%
  mutate_at(vars(-one_of('um_per_px')), ~. * !!sym('um_per_px'))


## RUN SUMMARY ----

# Question about how many spots fit certain conditions and modeling constraints
# For example, what is the breakdown of tumor distances when we start censoring spots for adjacency to other biological structures?
param.edge <- 200 # distance in um from tissue border to exclude from analysis
param.other.dist <- 100 # minimum distance in um from other tissue structures for inclusion
param.other.class <- c("Pathology.Distance.TLSImmature.unfilled.invertFALSE",
                       "Pathology.Distance.TLSMature.unfilled.invertFALSE",
                       "Pathology.Distance.Vessel.unfilled.invertFALSE",
                       "Pathology.Distance.BenignEpi.unfilled.invertFALSE",
                       "Pathology.Distance.ExoEndo.unfilled.invertFALSE",
                       "Pathology.Distance.LuminalNec.unfilled.invertFALSE",
                       "Pathology.Distance.NormalAdj.unfilled.invertFALSE")
param.binw <- 50 # bin width (um length) for histogram
param.max <- 1000 # max distance to consider in um

# let's do a param sweep of other.dist to see how many spots fit the criteria
filtered_stromal_spots_df <- function(in_param.other.dist, in_param.edge, in_param.max, in_param.other.class, in_df){
  filt_df <- in_df %>%
    # swap NA with Inf so cases without a particular features aren't filtered out
    dplyr::mutate_all(~replace(., is.na(.), Inf)) %>%
    # apply other class minimum distance filter
    dplyr::filter_at(vars(all_of(in_param.other.class)), all_vars(. > in_param.other.dist)) %>%
    # apply minimum distance from tissue sample edge
    dplyr::filter(Pathology.Distance.Tissue.filled.invertTRUE > in_param.edge) %>%
    # 'stroma' selection based on distances (without loading class from Seurat)
    dplyr::filter(Pathology.Distance.NonEpi.unfilled.invertFALSE < 50) %>%
    # apply maximal distance to consider
    dplyr::filter(Pathology.Distance.Tumor.unfilled.invertFALSE < in_param.max)
  return(filt_df)
}
count_filtered_spots <- function(in_param.other.dist, in_param.edge, in_param.max, in_param.other.class, in_df){
  spot_df <- filtered_stromal_spots_df(in_param.other.dist, in_param.edge, in_param.max, in_param.other.class, in_df)
  return(nrow(spot_df))
}
# function to determine how many individual samples remain in filtered spots
count_filtered_samples <- function(in_param.other.dist, in_param.edge, in_param.max, in_param.other.class, in_df){
  filt_df <- filtered_stromal_spots_df(in_param.other.dist, in_param.edge, in_param.max, in_param.other.class, in_df)
  # get sample ID from row names
  sample_df <- filt_df %>%
    rownames_to_column(var = "Barcode") %>%
    mutate(Sample_Name = gsub("(.*)\\:(.*)", "\\1", Barcode)) %>%
    group_by(Sample_Name) %>%
    summarise(spot_count=n())
  return(sample_df)
}

# apply the spot count function to all combinations of parameters
param_search_df <- expand.grid(param.edge=seq(0,2000,50), param.max=seq(500,2000,100), param.other.dist=seq(0,1000,50))
f = mapply(count_filtered_spots, param_search_df$param.other.dist, param_search_df$param.edge, param_search_df$param.max, MoreArgs = list(param.other.class, scaled_distances_df))
param_search_df <- cbind(param_search_df, f)

# plot the resulting spot counts
my_breaks = c(0, 1, 100, 1000, 10000, 50000, 100000, 250000)
param_search_df %>%
  dplyr::filter(param.max == 2000) %>%
  ggplot(aes(x=param.edge, y=param.other.dist, fill=f)) +
  geom_raster(hjust = 1, vjust=0) +
  ggtitle("Count of Spots With Various Distance Constraints") +
  xlab("Minimum Distance from Tissue Edge (um)") +
  ylab("Minimum Distance from Other Classes (um)") +
  scale_fill_viridis_c(breaks = my_breaks, labels = my_breaks,
                       trans = scales::pseudo_log_trans(sigma = 0.001),
                       guide = guide_colorbar(title = "Spots (n)")) +
  theme_bw(16) +
  theme(legend.position = "right", legend.key.height = unit(2.5, "cm"))
  
# previous plot suggests 500um as a upper limit to the other class distance parameter
param_search_df %>%
  dplyr::filter(param.max == 2000) %>%
  ggplot(aes(x=param.edge, y=param.other.dist, fill=f)) +
  geom_raster(hjust = 1, vjust=0) +
  ylim(c(0, 500)) +
  ggtitle("Count of Spots With Various Distance Constraints") +
  xlab("Minimum Distance from Tissue Edge (um)") +
  ylab("Minimum Distance from Other Classes (um)") +
  scale_fill_viridis_c(breaks = my_breaks[1:(length(my_breaks)-2)], labels = my_breaks[1:(length(my_breaks)-2)],
                       trans = scales::pseudo_log_trans(sigma = 0.001),
                       guide = guide_colorbar(title = "Spots (n)")) +
  theme_bw(16) +
  theme(legend.position = "right", legend.key.height = unit(2.5, "cm"))

# previous plot suggests 500um as a upper limit to the other class distance parameter (non-log)
param_search_df %>%
  dplyr::filter(param.max == 2000) %>%
  ggplot(aes(x=param.edge, y=param.other.dist, fill=f)) +
  geom_raster(hjust = 1, vjust=0) +
  ylim(c(0, 200)) +
  ggtitle("Count of Spots With Various Distance Constraints") +
  xlab("Minimum Distance from Tissue Edge (um)") +
  ylab("Minimum Distance from Other Classes (um)") +
  scale_fill_viridis_c(breaks = c(0,50000,100000), limits = c(0, 100000), labels = c('0','50k','100k'), guide = guide_colorbar(title = "Spots (n)")) +
  theme_bw(16) +
  theme(legend.position = "right", legend.key.height = unit(2.5, "cm"))


# the plot shows lack of sensitivity to tissue edge distance
hist(scaled_distances_df$Pathology.Distance.Tissue.filled.invertTRUE)



# summarize sample representation for a given spot count
spots_by_sample_df <- count_filtered_samples(100, 200, 1000, param.other.class, scaled_distances_df)
# quick plot
spots_by_sample_df %>%
  ggplot() +
  geom_col(aes(x = reorder(Sample_Name, spot_count), y = spot_count)) +
  xlab("Sample") +
  ylab("Spot Count (after filtering)") +
  ggtitle(paste0("Count of Stroma Spots \n(",sum(spots_by_sample_df$spot_count)," spots, edge offset: 200um, class offset: 100um, max distance: 1mm)")) +
  theme_bw(12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





plt_df <- scaled_distances_df %>%
  dplyr::mutate_all(~replace(., is.na(.), Inf)) %>%
  dplyr::filter_at(vars(param.other.class), all_vars(. > param.other.dist)) %>%
  dplyr::filter(Pathology.Distance.Tissue.filled.invertTRUE > param.edge) %>%
  # 'stroma' selection based on distances (without loading class from Seurat)
  dplyr::filter(Pathology.Distance.NonEpi.unfilled.invertFALSE < 50)

ggplot(plt_df) +
  geom_histogram(aes(x = Pathology.Distance.Tumor.unfilled.invertFALSE), breaks = seq(0, param.max, param.binw)) +
  xlab("NonEpi Spot Distance to Tumor (um)") +
  ylab("Spot Count") +
  ggtitle(paste0("Spot Distances of NonEpi vs Tumor\n(",nrow(plt_df)," spots, edge offset: ",param.edge,"um, class offset: ",param.other.dist,"um)")) +
  theme_bw(14)















## INGEST SEURAT OBJECT ----

# Open a connection to the RDS object
con <- gzfile(input_path)

# Ingest the Seurat object
data_seurat <- readRDS(con)

# Close the connection to be able to overwrite
close(con)


# Merge the distance values into the seurat features table per sample
# add to meta.data table
data_seurat@meta.data <- merge(data_seurat@meta.data, distances_df, by="row.names", all=FALSE, all.x=TRUE, all.y=FALSE)



## PROCESS SEURAT OBJECT ----
## Calculate signature scores and add them to Seurat meta data ----


# Find genes abundant in this sample
#gene_list <- seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold,
#                                             assay="SCT", slot="data")
# above function didn't work, just use the 1000 list already contained in SCT
gene_list <- data_seurat@assays$SCT@var.features

# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select, 
                                            sig_length_min=1, sig_length_max=1000, ratio_threshold=0)

# Signature names in the object
sig_names <- names(signature_list)

# Invert signature list
sig_invert_list <- invert_list_my(signature_list)

# Add signature scores to a seurat object
data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.", assay="SCT")

# Seurat renames column names, this step finds new signature names
sig_names_upd <- colnames(data_seurat@meta.data)[grep("sig.", colnames(data_seurat@meta.data))]







## SAVE SEURAT FILE ----

if(do_rds_output == TRUE){
  # Write the updated Seurat object
  filename <- file.path(output_path, basename(input_path))
  saveRDS(data_seurat, file=filename)
}













for (row_i in 1:nrow(meta_df)){
  print(row_i)
  
  # Sample name
  sample_name <- meta_df[[row_i, "Sample_Name"]]
  # Trichrome file to integrate
  trichrome_job_id <- dplyr::filter(halo_summary_df, Analysis.Region==sample_name)[["Job.Id"]]
  trichrome_file <- dir(file.path(trichrome_path, "ObjectData"), pattern = paste0('job',trichrome_job_id), full.names = TRUE, recursive = TRUE)
  
  # Find the RDS object
  file_list <- dir(input_path, pattern=paste0(sample_name), full.names=TRUE, recursive=TRUE)
  
  # Ensure singular matching of found files
  if(length(file_list) == 1 & length(trichrome_file) == 1){
    
    print(file_list[1])
    print(trichrome_file[1])
    
    ## INGEST SEURAT DATA ----
    
    
    # Open a connection to the RDS object
    con <- gzfile(file_list[1])
    
    # Ingest the Seurat object
    data_seurat <- readRDS(con)
    
    # Close the connection to be able to overwrite
    close(con)
    
    
    ## EXTRACT SEURAT SPOT DATA ----
    
    
    # Seurat image object
    image_structure <- data_seurat@images[[1]]
    
    # Seurat spot center coordinates in full resolution: Y = imagerow, X = imagecol
    X <- image_structure@coordinates[["imagecol"]]
    Y <- image_structure@coordinates[["imagerow"]]
    
    
    ## LOAD TRICHROME DATA ----
    
    
    # precalc the center of detected analysis objects
    halo_obj_df <- read.csv(trichrome_file) %>%
      dplyr::mutate(Visium.Registered.X = (Visium.Registered.XMin + Visium.Registered.XMax)/2) %>%
      dplyr::mutate(Visium.Registered.Y = (Visium.Registered.YMin + Visium.Registered.YMax)/2) %>%
      dplyr::mutate(X = mean(c(XMin, XMax))) %>%
      dplyr::mutate(Y = mean(c(YMin, YMax)))
      
    
    ## MAP TRICHROME OBJECTS TO SPOTS ----
    
    
    d <- raster::pointDistance(halo_obj_df[, c("Visium.Registered.X","Visium.Registered.Y")], cbind(X,Y), lonlat = FALSE)
    r <- apply(d, 1, which.min)
    # add spot call to the halo table
    halo_obj_df$Visium.Spot <- row.names(data_seurat@meta.data)[r]
    
    # now that we have spot labels, summarize halo table per spot
    # note collagen (blue) is stain 1 in results
    # note Tissue Area seems to be filtered for background (Region Area is not)
    var_list <- c("Tissue.Area.Analyzed..um.2.",
                  "Stain.1.Area..um.2.",
                  "Stain.1.Weak.Area..um.2.",
                  "Stain.1.Moderate.Area..um.2.",
                  "Stain.1.Strong.Area..um.2.")
    
    halo_spot_df <- halo_obj_df %>%
      dplyr::group_by(Visium.Spot, Classifier.Label) %>%
      dplyr::summarize_at(dplyr::all_of(var_list), sum)
    
    # Now how to treat Visium spots with no trichrome data?
    # true zeros need to be separate from spots lacking coreg tissue
    # HALO analysis modified to provide non-stromal objects (glass, tumor)
    # will need to zero out measures in tumor and glass (not in alg scope)
    # however keep Tumor area for collagen % over tissue measures
    halo_spot_df[halo_spot_df$Classifier.Label=='Tumor', 4:ncol(halo_spot_df)] <- 0
    halo_spot_df[halo_spot_df$Classifier.Label=='Glass', 3:ncol(halo_spot_df)] <- 0
    halo_spot_df$not_empty <- ifelse(halo_spot_df$Classifier.Label=="Glass", 0, 1)
    
    # compute prop stromal value (note tumor and glass will remain zero)
    halo_spot_df <- halo_spot_df %>%
      dplyr::mutate(Proportion.Stroma.Collagen.Positive = Stain.1.Area..um.2./Tissue.Area.Analyzed..um.2.) %>%
      dplyr::mutate(Proportion.Stroma.Collagen.Weak.Positive = Stain.1.Weak.Area..um.2./Tissue.Area.Analyzed..um.2.) %>%
      dplyr::mutate(Proportion.Stroma.Collagen.Moderate.Positive = Stain.1.Moderate.Area..um.2./Tissue.Area.Analyzed..um.2.) %>%
      dplyr::mutate(Proportion.Stroma.Collagen.Strong.Positive = Stain.1.Strong.Area..um.2./Tissue.Area.Analyzed..um.2.)
    
    # now Tumor and Glass have zero measures, drop classifier label and condense again
    # ignore NaN so glass values don't cancel out stromal prop values
    # use the not_empty column to drop spots that were entirely glass (no coreg tissue underneath)
    halo_spot_df <- halo_spot_df %>%
      dplyr::select(-Classifier.Label) %>%
      dplyr::group_by(Visium.Spot) %>%
      dplyr::summarize_if(is.numeric, ~sum(., na.rm = TRUE)) %>%
      dplyr::filter(not_empty > 0) %>%
      dplyr::select(-not_empty)
    
    # compute prop tissue (tumor+stroma) value
    halo_spot_df <- halo_spot_df %>%
      dplyr::mutate(Proportion.Tissue.Collagen.Positive = Stain.1.Area..um.2./Tissue.Area.Analyzed..um.2.) %>%
      dplyr::mutate(Proportion.Tissue.Collagen.Weak.Positive = Stain.1.Weak.Area..um.2./Tissue.Area.Analyzed..um.2.) %>%
      dplyr::mutate(Proportion.Tissue.Collagen.Moderate.Positive = Stain.1.Moderate.Area..um.2./Tissue.Area.Analyzed..um.2.) %>%
      dplyr::mutate(Proportion.Tissue.Collagen.Strong.Positive = Stain.1.Strong.Area..um.2./Tissue.Area.Analyzed..um.2.)
    
    # Add distance to the pathology feature as a meta data factor of the Seurat object
    trichrome_df <- halo_spot_df %>%
      dplyr::select(-Tissue.Area.Analyzed..um.2.) %>%
      tibble::column_to_rownames("Visium.Spot") %>%
      dplyr::rename_all(function(x) paste0("Trichrome.",x))
    
    # missing spots means no coreg tissue existed, should be NA in seurat
    data_seurat@meta.data <- merge(data_seurat@meta.data, trichrome_df, by="row.names", all.x = TRUE) %>%
      tibble::column_to_rownames("Row.names")
    
    
    ## VISUALIZE ----
    pathology_var_name <- "Trichrome.Proportion.Tissue.Collagen.Positive"
    p1 <- spatial_feature_plot_my(data_seurat, pathology_var_name, min.cutoff=0, max.cutoff=1, crop=TRUE)
    
    # Save distance visualization to file
    filename <- file.path(output_figs, paste0("Trichrome_", sample_name, "_", pathology_var_name))
    write_plot2file_my(p1, filename, num_row=1, num_col=1)
    
    
    ## OVERWRITE SEURAT FILE* ----
    # * or write to csv
    
    if(do_rds_output == TRUE){
      # Write the updated Seurat object
      filename <- file.path(output_path, paste0(sample_name, "_trichrome_pathology.rds"))
      saveRDS(data_seurat, file=filename)
    }
    if(do_csv_output == TRUE){
      # Write the distance data to a csv (only distance data)
      filename <- file.path(output_path, paste0(sample_name, "_trichrome_pathology.csv"))
      data_seurat@meta.data %>%
        dplyr::select(contains("Trichrome.")) %>%
        write.csv(file = filename, row.names = TRUE)
    }
    
  } else{
    warning("More than one file found!\n", file_list)
  } # end conditional check for singular file per sample
  
  
} # end looping through samples





