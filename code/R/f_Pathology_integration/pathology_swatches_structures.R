# Author: Anna Lyubetskaya. Date: 21-04-22
# Select specific spots based on meta data and create H&E "swatches"

# Important concern: The swatches extracted from different images are of a different size - matches the original microscope seetting of the image


## ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_image.R")
source("code/R/Utils/utils_10X_signatures.R")
source("code/R/Utils/utils_pathology_spots.R")

source("code/R/b_Cluster_Biomarkers/Seurat_cluster_utils.R")


## PARAMETERS ----


# Run name
run_name <- "Nerve"

# Path to processed Seurat data
sample_name <- "PDAC108_path14_5K_harmony"
sample_exclude <- NULL

# Name for the signature group being plotted
sig_select <- "PDAC.U.Nervous"  # e.g., NULL or "PDAC.collisson.classical", "PDAC.moffitt.activatedstroma"
sig_empirical <- FALSE  # If want to swap scores for bootstrapped p-values
sig_pval_threshold1 <- 0.2  # score, e.g. 0.5
sig_pval_threshold2 <- 0  # score, e.g. 0

# Name of the pathology field to filter by groups
pathology_select <- "Pathology.Nerve.percent"  # e.g., "Tumor"
pathology_threshold1 <- 50  # upper, 75
pathology_threshold2 <- 0  # lower, 25

# SCT and spot number thresholds for gene filtering
sct_threshold <- 0.5
spot_threshold <- 5

# Assay and slot to use to calculate signature score
assay <- "SCT"
slot <- "data"

# User clustering resolution
# E.g., integrated_snn_res.0.4 or user.Sample_Name
resolution_user <- "user.Sample_Name"

# Group swatches by
swatch_group <- resolution_user  # e.g., user.Sample_Name
patient_group <- "user.Sample_Name"  # currently only "user.Sample_Name" works

# Number of spots to return
topn_spots <- 10000
topn_structures <- 2
str_size_max <- 100
combine_str_vis <- FALSE

# Name of the classifier to use in the Seurat object meta data
classifier_name <- paste0("Swatches_", run_name)

# Subset to specific clusters for modeling if needed
cluster_select <- NULL  # c(3, 10, 0, 4, 6, 9)

# Remove spots with a lot of white space or blur
remove_bad_spots <- FALSE

# Coefficient to expand the intended spot patch
expand_coef <- 2

# Randomize spot selection after thresholding?
do_random <- FALSE

# Visualizes combined patches in particular configuration
user_nrow <- 1  # default = NULL

# Select most relevant section per patient
# default == NULL, otherwise user.Block_ID or user.Region_ID
select_rep_section <- NULL  # "user.Block_ID"

# Check spots for contiguity?
do_contiguity <- TRUE

# Explicit boolean requesting rendering of swatches
draw_swatches <- TRUE


## PATHS ----


# Path to a signature file
sig_path <- "data/import/Signatures/signatures_pdac_20230720.txt"

# Location of pre-processed data
input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, classifier_name, "_", sample_name, "_", run_name, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Seurat data
data_seurat <- readRDS(input_path)

# Remove samples if necessary
if(!is.null(sample_exclude)){
  barcode_list <- rownames(data_seurat@meta.data[which(!data_seurat@meta.data[["user.Sample_Name"]] %in% sample_exclude),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# Remove clusters if necessary
if(!is.null(cluster_select)){
  barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data[[resolution_user]] %in% cluster_select),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# Remove bad spots by pathology if necessary
if(remove_bad_spots == TRUE){
  barcode_list <- rownames(data_seurat@meta.data)[which(data_seurat@meta.data$Pathology.White_Space.percent <= 10 & data_seurat@meta.data$Pathology.Blur_Tissue.percent <= 10)]
  if(length(barcode_list) > 0){
    data_seurat <- subset(data_seurat, cells=barcode_list)
  }
}

# Select a single representative section for every patient
if(!is.null(select_rep_section)){
  section_select <- data_seurat@meta.data %>% 
    tibble::rownames_to_column("Coordinate")
  
  if("user.Sample_Name" != swatch_group){
    section_select <- section_select %>% 
      dplyr::group_by(user.Sample_Name, user.Region_ID, user.Block_ID, !!rlang::sym(swatch_group))
  } else{
    section_select <- section_select %>% 
      dplyr::group_by(user.Sample_Name, user.Region_ID, user.Block_ID)
  }
  
  section_select <- section_select%>% 
    dplyr::summarize(CountSpots = dplyr::n_distinct(Coordinate)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(user.Sample_Name, !!rlang::sym(select_rep_section)) %>%
    dplyr::summarise(CountGroups = dplyr::n_distinct(!!rlang::sym(swatch_group)),
                     CountSpots = sum(CountSpots)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(CountGroups), desc(CountSpots)) %>%
    dplyr::group_by(!!rlang::sym(select_rep_section)) %>%
    dplyr::top_n(1)
  
  barcode_list <- rownames(data_seurat@meta.data)[which(data_seurat@meta.data$user.Sample_Name %in% section_select$user.Sample_Name)]
  if(length(barcode_list) > 0){
    data_seurat <- subset(data_seurat, cells=barcode_list)
  }
  
}


## CALCULATE TARGET SIGNATURE SCORES ----


## Calculate signature scores and add them to Seurat meta data

# Find genes abundant in this sample
gene_list <- seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold, assay=assay, slot=slot)

# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=c(sig_select), sig_length_min=1, sig_length_max=1000, ratio_threshold=0)
cat("Target signature length =", length(signature_list[[1]]))

# Add signature scores to a seurat object
data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.", assay=assay)

# Find column names in the Seurat object - Seurat can rename original signatures
sig_names <- colnames(data_seurat@meta.data[grep("sig.", colnames(data_seurat@meta.data))])

# Column name for the signature empirical p-value if available
if(sig_empirical == FALSE){
  sig_name_pval <- sig_names
}else {
  sig_name_pval <- paste0(sig_names, ".EPvalue")  
}


## CREATE A RANDOM SIGNATURE SCORE BACKGROUND ----


if(sig_empirical == TRUE){
  # Write/read random gene expression scores to file
  filename <- paste0(output_path, "/sig_random_", sample_name, "_", sig_names, "_", pathology_select, ".txt")
  
  # Create a random signature background distribution and score the target signature against it
  data_seurat <- signature_empirical_pvalue_my(data_seurat, signature_list[[1]], output_path, sample_name, sig_names[1], sig_random_filename=filename,
                                               n_simulations=1000, sct_threshold=1, spot_threshold=10, assay="SCT", slot="data", output_col_name=sig_name_pval)
}


## WRANGLE META DATA ----


# Extract coordinate meta data
meta_init_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

# Add user.Sample_Name to meta data if absent: happens in stand-alone samples
if(!"user.Sample_Name" %in% colnames(meta_init_df)){
  meta_init_df["user.Sample_Name"] = data_seurat@misc$user.Sample_Name
}

# Select relevant meta data
meta_df <- meta_init_df  %>%
  dplyr::select(dplyr::all_of(unique(c("Coordinate", "user.Sample_Name", patient_group, sig_name_pval, pathology_select, resolution_user, swatch_group))))

# Anonymize two main parameters by which spots are selected
colnames(meta_df) <- gsub(sig_name_pval, "Signature", colnames(meta_df))
colnames(meta_df) <- gsub(pathology_select, "Pathology", colnames(meta_df))


## VISUALIZE DISTRIBUTIONS ----


# Find mean signature and pathology score by cluster
if("user.Sample_Name" == swatch_group){
  clust_df <- meta_df %>%
    dplyr::group_by(!!rlang::sym(swatch_group)) %>%
    dplyr::summarise(SignatureMean = mean(Signature),
                     PathologyMean = mean(Pathology))
} else{
  clust_df <- meta_df %>%
    dplyr::group_by(!!rlang::sym(patient_group), !!rlang::sym(swatch_group)) %>%
    dplyr::summarise(SignatureMean = mean(Signature),
                     PathologyMean = mean(Pathology))
}

# Define custom colors
cols <- define_cols_my(n=length(unique(clust_df[[swatch_group]])))
names(cols) <- sort(unique(clust_df[[swatch_group]]))

if(resolution_user %in% colnames(clust_df)){
  # Plot cluster mean signature score v cluster mean pathology score  
  filename <- paste0(output_path, "scatter_clust_", sample_name, "_", sig_select, "_", pathology_select)
  p <- create_scatter_plot_my(clust_df, x_label="PathologyMean", y_label="SignatureMean", 
                              fill_label=resolution_user, facet_var=c(patient_group, "free_y"), 
                              filename=filename, size=5, labels=c(pathology_select, sig_select, sample_name), cols=cols)
}

# Plot spot signature score v spot pathology score
filename <- paste0(output_path, "scatter_spot_", sample_name, "_", sig_select, "_", pathology_select)
p <- create_scatter_plot_my(meta_df, x_label="Pathology", y_label="Signature", 
                            fill_label=resolution_user, facet_var=c(patient_group, "fixed"), 
                            filename=filename, size=0.5, labels=c(pathology_select, sig_select, sample_name), stroke=0, cols=cols)


## SELECT SPOTS ----


extract_topn_my <- function(meta_df, topn_spots=5, category="", do_random=FALSE){
  # Extract top N hits in a tibble based on two parameters
  
  if(do_random == TRUE){
    meta_df <- meta_df[sample(1:nrow(meta_df)),]
  }
  
  meta_df <- meta_df %>%
    dplyr::slice_head(n = topn_spots) %>%
    dplyr::mutate(Category = category)
  
  cat(category, nrow(spot_lists[[category]]), "\n")
  
  return(meta_df)
  
}


spot_lists <- list()

# Threshold each sample individually
for(s in sort(unique(meta_df[[swatch_group]]))){
  
  meta_loc_df <- meta_df %>%
    dplyr::filter(!!rlang::sym(swatch_group) == s)
  
  # Pathology high, signature high
  category <- "Hi-Hi"
  
  spot_lists[[paste(category, s)]] <- meta_loc_df %>%
    dplyr::filter(Pathology >= pathology_threshold1 & Signature >= sig_pval_threshold1) %>%
    dplyr::arrange(desc(Pathology), desc(round(Signature, 1)), !!rlang::sym(patient_group)) %>%
    extract_topn_my(topn_spots=topn_spots, category=category, do_random=do_random)
  
  
  # Pathology low, signature high
  category <- "Lo-Hi"
  spot_lists[[paste(category, s)]] <- meta_loc_df %>%
    dplyr::filter(Pathology < pathology_threshold2 & Signature >= sig_pval_threshold1) %>%
    dplyr::arrange(Pathology, desc(round(Signature, 1)), !!rlang::sym(patient_group)) %>%
    extract_topn_my(topn_spots=topn_spots, category=category, do_random=do_random)
  
  # Pathology high, signature low
  category <- "Hi-Lo"
  spot_lists[[paste(category, s)]] <- meta_loc_df %>%
    dplyr::filter(Pathology >= pathology_threshold1 & Signature < sig_pval_threshold2) %>%
    dplyr::arrange(desc(Pathology), round(Signature, 1), !!rlang::sym(patient_group)) %>%
    extract_topn_my(topn_spots=topn_spots, category=category, do_random=do_random)
  
  # Pathology low, signature low
  category <- "Lo-Lo"
  spot_lists[[paste(category, s)]] <- meta_loc_df %>%
    dplyr::filter(Pathology < pathology_threshold2 & Signature < sig_pval_threshold2) %>%
    dplyr::arrange(Pathology, round(Signature, 1), !!rlang::sym(patient_group)) %>%
    extract_topn_my(topn_spots=topn_spots, category=category, do_random=do_random)
}


# Tibble of selected spots
spot_df <- dplyr::bind_rows(spot_lists)


## VISUALIZE SPOTS ----


# Add categories to the Seurat object
data_seurat@meta.data[classifier_name] <- "-"
for(spot_name in names(spot_lists)){
  if(nrow(spot_lists[[spot_name]]) > 1){
    data_seurat@meta.data[spot_lists[[spot_name]]$Coordinate, classifier_name] <- gsub(" .+", "", spot_name)
  }
}
table(data_seurat@meta.data[classifier_name])

# Assign these colors to the categories
cols <- list("Hi-Hi" = "red", "Hi-Lo" = "orange", "Lo-Hi" = "darkblue", "Lo-Lo" = "grey", "-" = "white")

# Plot selected spots
# For a given variable, plot PCA, UMAP, and spatial distributions
filename <- paste0(output_path, paste(c("spatial", sig_select, pathology_select), collapse="_"))
Seurat_pca_umap_spatial_my(data_seurat, classifier_name, filename, cols=cols)


## CHECK SWATCHES FOR CONTIGUITY ----


# Setup a tibble to collect spot contiguity
cluster_id <- spot_df %>%
  dplyr::mutate(Cluster_ID = "") %>%
  dplyr::select(Coordinate, Cluster_ID) %>%
  unique() %>%
  tidyr::drop_na() %>%
  tibble::column_to_rownames("Coordinate")

# Go through samples and check swatch groups for neighbors
if(do_contiguity == TRUE){
  for(group in sort(unique(spot_df[[swatch_group]]))){
    
    for(category in unique(spot_df$Category)){
      
      # Grab swatches in the group
      spot_loc_df <- spot_df %>%
        dplyr::filter(Category == category & !!rlang::sym(swatch_group) == group)
      
      
      for(sample in unique(spot_loc_df[[patient_group]])){
        
        # Extract spot information for a specific group and specific sample
        coordinate_list <- spot_loc_df %>%
          dplyr::filter(!!rlang::sym(patient_group) == sample) %>%
          dplyr::pull(Coordinate)
        
        
        ## A version of seurat_spot_neigbors_my function in utils_10X_image.R
        
        # Extract Seurat coordinates
        xy <- data_seurat@images[[sample]]@coordinates[coordinate_list, c("row", "col")]
        barcode_list <- rownames(xy)
        xy <- mapply(xy, FUN=as.numeric)
        
        # Find the matrix with the indices of points belonging to the set of the k nearest neighbours of each other
        if(length(barcode_list) > 1){
          spot_nb <- spdep::dnearneigh(xy, 0, 2, row.names=barcode_list)
          link_check <- max(sapply(spot_nb, function(x) max(x)))
          
          
          if(link_check > 1){
            
            ## A version of seurat_image_to_clusters_my function in utils_10X_image.R
            
            # Transform neighborhoods into a weighted list
            spot_nb_list <- as(spdep::nb2mat(spot_nb, style="B", zero.policy=TRUE), "CsparseMatrix")
            
            # Calculate graph adjacency
            ## set add.rownames = NULL to pass the cell ID to the output dgcMatrix
            spot_graph <- igraph::graph.adjacency(spot_nb_list, mode="undirected", add.rownames = NULL)
            
            # Identify graph subclusters
            spot_clusters <- igraph::clusters(spot_graph)
            
            # Select barcodes corresponding to the cluster      
            for(cluster in unique(spot_clusters$membership)){
              barcode_list <- names(spot_clusters$membership)[which(spot_clusters$membership == cluster)]
              cluster_id[barcode_list, "Cluster_ID"] <- cluster
            }
          } else{
            last_cluster <- as.numeric(max(cluster_id$Cluster_ID))
            if(is.na(last_cluster)){
              last_cluster <- 0
            }
            cluster_id[barcode_list, "Cluster_ID"] <- as.character((last_cluster+1):(last_cluster+length(barcode_list)))
          }
        } else{
          cluster_id[barcode_list, "Cluster_ID"] <- "1"
        }
        
      }
    }
  }
} else{
  cluster_id[["Cluster_ID"]] <- as.character(1:nrow(cluster_id))
}

# Add contiguity information to the spot tibble
spot_df <- spot_df %>%
  dplyr::inner_join(cluster_id %>%
                      tibble::rownames_to_column("Coordinate"), by="Coordinate")

# Write to file
filename <- paste0(output_path, paste(c("table", sig_select, pathology_select), collapse="_"), ".txt")
readr::write_delim(spot_df, filename, delim="\t")


# Count the number of structures per sample
sample_df <- spot_df %>%
  dplyr::select(user.Sample_Name, Category, Cluster_ID) %>%
  unique() %>%
  dplyr::group_by(user.Sample_Name, Category) %>%
  dplyr::summarise(StructureCount = dplyr::n_distinct(Cluster_ID)) %>%
  dplyr::ungroup()

# Write to file
filename <- paste0(output_path, paste(c("table_sample", sig_select, pathology_select), collapse="_"), ".txt")
readr::write_delim(sample_df, filename, delim="\t")


# Count the number of structures per sample
patient_df <- sample_df %>%
  dplyr::inner_join(meta_init_df %>%
                      dplyr::select(user.Sample_Name, user.Block_ID) %>%
                      unique(), by="user.Sample_Name") %>%
  dplyr::group_by(user.Block_ID, Category) %>%
  dplyr::summarise(StructureCount = sum(as.numeric(StructureCount)))

# Write to file
filename <- paste0(output_path, paste(c("table_patient", sig_select, pathology_select), collapse="_"), ".txt")
readr::write_delim(patient_df, filename, delim="\t")


## CREATE SWATCHES ----


# Concentriq image information
image_df <- meta_init_df %>% 
  dplyr::select(dplyr::all_of(unique(c("user.Sample_Name", patient_group, "user.Concentriq_Image_ID", "user.Concentriq_Repo_ID")))) %>% 
  unique()


if(draw_swatches == TRUE){
  # Keep log of swatches
  log_file <- paste0(output_path, "log.txt")
  write("Group\tPatch Category\tSample\tCluster\tRow (Y)\tCol (X)\tRadius", log_file, append=FALSE)
  
  # Go through samples and extract swatches
  for(group in sort(unique(spot_df[[swatch_group]]))){  # Iterate over a group (clusters or pathology)
    
    for(category in unique(spot_df$Category)){  # Iterate over 4 corners
      
      swatch_combo_list <- list()
      
      for(sample in sort(unique(spot_df[[patient_group]]))){ # Iterate over samples
        
        # Grab swatches in the group
        spot_loc_df <- spot_df %>%
          dplyr::filter(Category == category & 
                          !!rlang::sym(swatch_group) == group & 
                          !!rlang::sym(patient_group) == sample)
        
        # List of structures to visualize starting from largest
        clust_list <- spot_loc_df %>%
          dplyr::group_by(Cluster_ID) %>%
          dplyr::summarize(n = dplyr::n_distinct(Coordinate)) %>%
          dplyr::filter(n <= str_size_max) %>%
          dplyr::arrange(desc(n)) %>%
          dplyr::slice_head(n = topn_structures) %>%
          dplyr::pull(Cluster_ID)
        
        for(cl in clust_list){
          
          cat(group, "\t", category, "\t", sample, "\t", cl, "\n")
          
          swatch_list <- list()
          
          # Extract spot information for a specific group and specific sample
          coordinate_list <- spot_loc_df %>%
            dplyr::filter(!!rlang::sym(patient_group) == sample & Cluster_ID == cl) %>%
            dplyr::pull(Coordinate)
          
          # Extract sample information for a specific group and specific sample / patient
          section_name <- spot_loc_df %>%
            dplyr::filter(!!rlang::sym(patient_group) == sample & Cluster_ID == cl) %>%
            dplyr::pull(user.Sample_Name) %>%
            unique()
          
          
          # Concentriq image ID
          concetriq_id <- image_df %>%
            dplyr::filter(user.Sample_Name == section_name) %>%
            dplyr::pull(user.Concentriq_Image_ID)
          
          if(length(coordinate_list) >= 1){
            
            # Collect a patch for every coordinate
            # Notice the reversal of column order! col = X; row = Y
            xy_list <- data_seurat@images[[sample]]@coordinates[coordinate_list, c("imagecol", "imagerow")]
            rad <- data_seurat@images[[sample]]@scale.factors$spot_diameter_fullres / 2
            
            write(paste(c(group, category, sample, cl, xy_list, rad), collapse="\t"), log_file, append=TRUE)
            
            # Create patches
            swatch_list[[cl]] <- fetch_image_patch(in_concentriq_id=concetriq_id, 
                                                   xy_list=xy_list, in_radius_px=rad,
                                                   expand_val=rad*expand_coef,
                                                   vis_type="structure")
          }
          
          
          # Visualize structure
          if(length(swatch_list) > 0 && combine_str_vis == FALSE){
            title <- paste0(group, ".  ", run_name, ".  ", category, ".  ", sample)
            filename <- paste0(output_path, paste(c("swatch", group, run_name, category, sample, cl), collapse="_"), ".png")
            
            patch_image_list <- combine_and_vis_patches_my(swatch_list, title=title, filename=filename,
                                                           do_resize=FALSE)
            
          } else{
            swatch_combo_list <- c(swatch_combo_list, swatch_list)
          }
          
        }
      }
      
      # If all swatches are a single-spot patches, combine them together
      if(length(swatch_combo_list) > 0 && combine_str_vis == TRUE){
        title <- paste0(group, ".  ", run_name, ".  ", category)
        filename <- paste0(output_path, paste(c("swatch", group, run_name, category), collapse="_"), ".png")
        
        if(is.null(user_nrow)){
          user_nrow <- floor(log2(length(swatch_combo_list)))
        }
        
        patch_image_list <- combine_and_vis_patches_my(swatch_combo_list, title=title, filename=filename, 
                                                       do_resize=TRUE, nrow=user_nrow)
        
      }
      
    }
  }
}
