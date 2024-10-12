# Author: Anna Lyubetskaya. Date: 20-08-06

# Plot a set of signatures for a set of Seurat RDS objects
# Produce Spatial plots in various configurations (PCA / UMAP are also available)
# Allows grouping plots by sample or signature as well as picking various scaling options

# Important: Make sure that signatures will not be renamed on ingestion (their names are R compliant)

# Bug to fix: plotting_types list is not currently subsetted with signatures


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_signatures.R")

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_signatures.R")
source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_10X_vis.R")

library(Seurat)


## PARAMETERS ----


# Run ID
run_name <- "exhaust"

# The cohort of interest: regex and name
cohort_name <- "PDAC"
cohort_regex <- "PDAC108_path14_5K_harmony"

# What to plot on one plot
# - all signatures ("signatures"), samples split into individual plots
# - all samples ("samples"), signatures split into individual plots
# - nothing ("nothing"), plot everything individually but with a PCA and a UMAP
plot_together <- "signatures"

# Do scale plot
do_scale_plot <- FALSE

# Create PCA/UMAP
do_dim_plot <- FALSE

# Quartile parameters
min_cutoff <- "q1"
max_cutoff <- "q95"

# Gene abundance filters
sct_threshold <- 0.5
spot_threshold <- 5

# From which assay and slot to plot data
assay <- "SCT"
slot <- "data"

# Signature filter: filter out signature that don't ever reach this intensity
sig_threshold <- 0.1

# Specific signatures to use
sig_select <- c("PDAC.P19.Bcell", "PDAC.P19.Tcell")

# Or single genes to plot
gene_interest_list <- c("CD38", "CD4", "CD8A", "CTLA4", "ENTPD1", "FOXP3", 
                        "HAVCR2", "ICOS", "IL2RA", "LAG3", "PDCD1", "TIGIT", "TOX")

# Or a meta data field
meta_cont_list <- NULL
meta_cat_list <- c("Pathology.Group", "integrated_snn_res.0.1")

# In case the input is a single integrated object and user wants to break it up into individual samples for more flexibility
breakup_object <- TRUE

# Create an underlying tissue plot as a first plot
plot_HE <- TRUE


## PATHS ----


# Input folder
# Location of pre-processed data
input_paths <- c("XXXX")

# Output path
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_name, "_", run_name, "/")

# Create output folder
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)

# Path to a signature file
sig_path <- "data/import/Signatures/signatures_220819.txt"


## INGEST DATA ----


# Find all files following the regex
file_list <- unname(unlist(sapply(input_paths, function(p) dir(p, pattern=paste0(cohort_regex, ".*.rds"), full.names=TRUE))))

# Ingest a set of RDS Seurat objects
seurat_list <- read_rds_list_simple_my(file_list)

# In case the input is a single integrated object and user wants to break it up into individual samples for more flexibility
if(breakup_object == TRUE && length(seurat_list) == 1){
  
  # List to temporary store subset objects
  seurat_subset_list <- list()
  
  # Find samples and corresponding barcodes to subset
  subset_df <- tibble::as_tibble(seurat_list[[1]]@meta.data, rownames="Coordinate") %>%
    dplyr::select(Coordinate, user.Sample_Name)
  
  for(s in unique(subset_df$user.Sample_Name)){
    # Barcodes for the select samples
    cell_list <- subset_df %>%
      dplyr::filter(user.Sample_Name == s) %>%
      dplyr::pull(Coordinate)
    
    # Subset object
    data_subset_seurat <- subset(seurat_list[[1]], cells=cell_list)
    
    # Only leave the relevant image
    data_subset_seurat@images <- data_subset_seurat@images[unique(data_subset_seurat@meta.data$user.Sample_Name)]
    
    seurat_subset_list[[s]] <- data_subset_seurat
  }
  
  seurat_list <- seurat_subset_list
  rm(seurat_subset_list)
  gc()
}

# Define feature set to plot
feature_set <- c(meta_cat_list, paste0("sig.", sig_select), meta_cont_list, gene_interest_list)
# Define plotting types from user inputs
plotting_types <- c(rep("d", length(meta_cat_list)), rep("f", length(sig_select)), rep("f", length(meta_cont_list)), rep("f", length(gene_interest_list)))


## CALCULATE SIGNATURE SCORE ----


# For each dataset, update the counts slot of SCT assay to be the same as the counts slot of Spatial assay
if(slot == "counts"){
  for(dataset in names(seurat_list)){
    seurat_list[[dataset]]@active.assay <- "Spatial"
  }
}

# Calculate signature scores and add them to Seurat meta data
if(!is.null(sig_path)){
  for(dataset in names(seurat_list)){
    # Find genes abundant in this sample
    gene_list <- seurat_select_abundant_genes_my(seurat_list[[dataset]], sct_threshold=sct_threshold, spot_threshold=spot_threshold)
    
    # Load signatures and filter them down to only well represented genes
    signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select,
                                                sig_length_min=1, sig_length_max=1000, ratio_threshold=0)
    
    # Add signature scores to a seurat object
    seurat_list[[dataset]] <- add_signature_scores_my(seurat_list[[dataset]], signature_list, prefix="sig.",
                                                      assay=assay)
  }
}


## VISUALIZE DATA ----


## Set a common color scheme for all categorical values across all samples
col_list <- list()
for(i in meta_cat_list){
  val_list <- sort(unique(unlist(unname(sapply(seurat_list, function(x) unique(x@meta.data[[i]]))))))
  
  col_list[[i]] <- define_cols_my(n=length(val_list))
  
  names(col_list[[i]]) <- val_list
}


# For each dataset, select two coordinates and put miminum and maximum values in there
# This doesn't work if the min_cutoff and max_cutoff is anything but q0 and q100
if(do_scale_plot == TRUE){
  # Column names across the set of samples
  col_names <- unique(c(unlist(unname(sapply(seurat_list, function(x) colnames(x@meta.data))))))
  # Grab signatures present in the set of samples
  sig_names_upd <- col_names[grepl("sig.", col_names)]
  
  
  for(sig in sig_names_upd){
    sig_min <- max(0, round(min( sapply(seurat_list, function(x) min(x@meta.data[[sig]])) ), 1))
    sig_max <- round(max( sapply(seurat_list, function(x) quantile(x@meta.data[[sig]], 0.85)) ), 1)
    
    for(dataset in names(seurat_list)){
      
      ## Cutoff signature outliers to global dataset min and max
      
      indx_min <- which(seurat_list[[dataset]]@meta.data[[sig]] < sig_min)
      indx_max <- which(seurat_list[[dataset]]@meta.data[[sig]] > sig_max)
      
      # Set minimum value
      seurat_list[[dataset]]@meta.data[indx_min, sig] <- sig_min
      # Set maximum value
      seurat_list[[dataset]]@meta.data[indx_max, sig] <- sig_max
      
      
      ## Set a single spot to a specific min and a max value to set the colorbar
      
      seurat_list[[dataset]] <- set_spatial_min_max_my(seurat_list[[dataset]], c(sig), min_val=sig_min, max_val=sig_max)
      
    }
  }
}


# Create spatial plots
if(plot_together == "signatures"){
  
  ## Plot all signatures for a single dataset on one plot
  
  # Go through samples
  for(dataset in names(seurat_list)){
    
    # Ensure default
    # Seurat::DefaultAssay(object = seurat_list[[dataset]]) <- "SCT"
    
    # Column names across the sample
    col_names <- c(colnames(seurat_list[[dataset]]@meta.data), rownames(seurat_list[[dataset]]))
    # Grab signatures present in the sample
    sig_names_upd <- intersect(feature_set, col_names)
    
    output_file <- paste0(output_path, "spatial_", dataset)
    if(!file.exists(output_file)){
      batch_spatial_feature_plot_my(seurat_list[dataset], sig_names_upd, output_file=output_file, plot_HE=plot_HE, 
                                    min.cutoff=min_cutoff, max.cutoff=max_cutoff, plot_type=plotting_types, col_list=col_list)
    }
    
  }
} else if(plot_together == "samples"){
  
  ## Plot all datasets for a single signature on one plot
  
  # Column names across the set of samples
  col_names <- unique(unlist(unname(sapply(seurat_list, function(x) colnames(x@meta.data)))))
  # Grab signatures present in the set of samples
  sig_names_upd <- intersect(feature_set, col_names)
  
  # Go through signatures
  for(sig in sig_names_upd){
    
    if("sig." %in% sig){
      title <- paste0(sig, "\n", paste0(signature_list[[gsub("sig.", "", sig)]], collapse="; "))
    } else{
      title <- sig
    }
    output_file <- paste0(output_path, "spatial_", sig)
    
    batch_spatial_feature_plot_my(seurat_list, sig, output_file=output_file, title=title,
                                  min.cutoff=min_cutoff, max.cutoff=max_cutoff, slot="counts")
    
  }
  
} else if(plot_together == "nothing"){
  
  ## Create individual plots for each marker to visualize
  
  for(dataset in names(seurat_list)){
    # Select a dataset
    data_seurat <- seurat_list[[dataset]]
    
    # Create a subfolder for organization
    dataset_name <- data_seurat@misc$user.Sample_Name
    out_path <- paste0(output_path, dataset_name,"/")
    dir.create(out_path, showWarnings = FALSE)
    
    # Seurat renames column names, this step finds new signature names
    # Column names across the sample
    col_names <- colnames(seurat_list[[dataset]]@meta.data)
    # Grab signatures present in the sample
    sig_names_upd <- intersect(feature_set, col_names)
    
    # Cycle through features to plot
    for(sig in sig_names_upd){
      # Proceed if the signature reaches the expression threshold
      if(max(data_seurat@meta.data[[sig]]) >= sig_threshold){
        
        # List of signature genes to add to the plot as caption
        caption <- paste0(sig, "\n", paste(signature_list[[gsub("sig.", "", sig)]], collapse=";"))
        
        if(do_dim_plot == TRUE){
          # Plot a classic two-dimensional projection: PCA
          p1 <- feature_plot_my(data_seurat, sig, reduction="pca")
          
          # Plot a classic two-dimensional projection: UMAP
          p2 <- feature_plot_my(data_seurat, sig, reduction="umap")
        }
        
        # Plot signature on the original tissue slice
        p3 <- spatial_feature_plot_my(data_seurat, sig, min.cutoff=min_cutoff, max.cutoff=max_cutoff, 
                                      name=paste(dataset_name, "\n", sig))
        
        # Write combo plot of PCA and UMAP to file          
        if(do_dim_plot == TRUE){
          filename <- paste0(output_path, "combo_", dataset_name, "_", sig)
          write_plot2file_my(patchwork::wrap_plots(list(p1, p2), nrow=1), filename, num_row=1, num_col=2)
        }
        
        # Write combo plot of tissue expression to file
        filename <- paste0(out_path, dataset_name, "_", sig, "_tissue")
        write_plot2file_my(p3 + 
                             patchwork::plot_annotation(caption = caption), 
                           filename, num_row=1, num_col=length(names(data_seurat@images)))
        
      }
    }
  }
}
