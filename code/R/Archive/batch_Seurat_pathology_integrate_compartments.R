#!/opt/tbio/domino_202108/binaries/R-4.0.5/bin/Rscript
# Author: Anna Lyubetskaya. Date: 21-01-18
# Add pathology compartments to a Seurat object


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(optparse))

source("/repos/P02567_TBIO-3021_10X_pilot/code/utils/utils_ggplot.R")
source("/repos/P02567_TBIO-3021_10X_pilot/code/R/Utils/utils_pathology_halo.R")
source("/repos/P02567_TBIO-3021_10X_pilot/code/R/Utils/utils_10X_image.R")
source("/repos/P02567_TBIO-3021_10X_pilot/code/R/b_Cluster_Biomarkers/Seurat_cluster_utils.R")


## PARAMETERS ----

parameter_list <- list(
  make_option(c("-c", "--cohort_name"),action="store", type="character", default="PDAC",
              help = "Cohort Regex"),
  make_option(c("-i", "--iterator"),action="store", type="integer", default=1,
              help = "Row number of metadata file to use")
)

# parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
# parser <- OptionParser(option_list = parameter_list)
opt <- parse_args(OptionParser(option_list=parameter_list))

# Sample / Cohort name
# cohort_name <- "PDAC"
cohort_name <- opt$c

# Baseline class
reference_class <- "Tissue"   # "Tissue" or NULL

# Rat colon
# remove_class_completely <- c("Tissue", "Connective_Debris")  # c("Tissue")
# remove_class_from_percent <- NULL  # c("Tissue")
# PDAC
# remove_class_completely <- c(NULL)
# remove_class_from_percent <- c("Tissue", "HumanPanc_ROI1_FFPE_B_Dec20", "HumanPanc_ROI1_FFPE_A_Apr21", 
#                                "HumanPanc_ROI2_FFPE_B_Dec20", "HumanPanc_ROI2_FFPE_B_Apr21",
#                                "Tumor epi", "Nontumor epi", "Viable_Tissue", "Blur_Tissue")
remove_class_completely <- NULL
# others <- c("learning islets", "learning tumor", "learning normal ducts", "learning tumor", "learning necrosis",
#             "learning stroma", "learning tumor", "learning exocrine", "learning - tumor", "learning - exocrine", 
#             "learning - benign", "learning - iselts", "learning - blood vessel", "learning - muscle-like from adjacent tissue (muscle)", 
#             "learning - benign (ducts)", "learning - benign (exocrine)", "learning iselts")
remove_class_from_percent <- c("Tissue", "Viable_Tissue", "Blur_Tissue", "Focused_Tissue","Adipose", "TLS-Mature", "TLS-Immature", "TLS-Aggregate", "Lymph_Node")
expected_classes <- c("Tissue", "Focused_Tissue", "Blue_Tissue", "Adipose", "TLS-Mature", "TLS-Immature", "TLS-Aggregate", "Lymph_Node", "Exocrine", "Benign Glands", "Blood", "Luminal Debris", "Blood Vessels", "Stroma", "Epithelium")
# Syng
# remove_class_completely <- c("Capsule", "OLD_B16", "OLD_MC38", "Tumor_Periphery", "Layer 3")
# remove_class_from_percent <- c("Tissue", "Tumor_Periphery", "Tumor_Core", "Melanin")
# HNSCC
# remove_class_completely <- c(NULL)
# remove_class_from_percent <- c("Tissue", "Viable_Tissue", "Blur_Tissue")

# Perform a subset of the Seurat object to only those spots that have a pathology annotation
# Subset necessitates re-clustering
do_subset <- TRUE
subset_threshold <- 10000

# Prefix to add to column names
col_prefix <- "Pathology."

# Don't write the final RDS object
no_rds_output <- FALSE

# Perform DEA or not
do_findmarkers <- FALSE

# User defined class order for plotting
class_order <- NULL
# RatColon
#class_order <- c("Outer_Muscular_Coat", "Myenteric_Plexus", "Inner_Muscular_Coat", "Submucosa", 
#                 "Mucosa_Muscularis",  "Crypt", "Lamina_Propria", "Colonocyte_Border")

# User defined colors when it has to be fully prescriptive
cols <- NULL
# PDAC
#cols <- list("Epithelium"="#c9415e", "NonEpi1_smooth"="#d3d3d3", "NonEpi2_rough"="#e5e5aa",
#             "Blood"="#992525", "Exocrine"="#4cb0e1", "Islets"="#153c65", "Luminal Debris"="#939393", "Stroma"="#e5e5aa")
# Syngeneics
#cols <- list("MC38" = "#db5a7c", "B16" = "#056db5", "Necrosis" = "#fae573")
# RatColon
#cols <- list("Outer_Muscular_Coat" = "#4dafdf", "Myenteric_Plexus" = "#fcb353", "Inner_Muscular_Coat" = "#153c65", "Submucosa" = "#b0e0ea", 
#             "Mucosa_Muscularis" = "#fae573",  "Crypt" = "#ff8086", "Lamina_Propria" = "#d3d3d3", "Colonocyte_Border" = "#96257d")

# Set seed for clustering
set.seed <- 531


## PATHS ----


# Input file with sample meta data
meta_file <- "XXXX"

# Input folder
input_path <- "XXXX"

# Output folder
output_path <- "XXXX"

# Output folder
output_init_figs <- "XXXX"

# Create output folders
dir.create(output_path, showWarnings = FALSE)
dir.create(output_init_figs, showWarnings = FALSE)


## INGEST DATA ----


# Read sample file; select only cohort-relevant samples
meta_df <- readr::read_delim(file=meta_file, delim="\t") %>%
  tidyr::drop_na() %>%
  dplyr::filter(grepl(cohort_name, Sample_Name) &
                  Pathology_Annotation_Name == "Compartments" &
                  Quality_Status == TRUE)


## WRANGLE DATA ----

i <- opt$i
# for(i in 1:nrow(meta_df)){
  # Sample name
  sample_name <- meta_df[[i, "Sample_Name"]]
  # Pathology file to integrate
  path_file <- meta_df[[i, "Pathology_Full_Path"]]
  # Any image rotation to apply
  pi_coef <- as.numeric(meta_df[[i, "Rotate"]])
  
  # Find the RDS object
  file_list <- dir(input_path, pattern=paste0(sample_name), full.names=TRUE, recursive=TRUE)
  
  if(length(file_list) == 1){
    
    # Create the output folder for all figures
    output_figs <- paste0(output_init_figs, sample_name, "/")
    dir.create(output_figs, showWarnings = FALSE)
    
    print(file_list[1])
    print(path_file)
    
    
    ## INGEST SEURAT DATA ----
    
    
    # Open a connection to the RDS object
    con <- gzfile(file_list[1])
    
    # Ingest the Seurat object
    data_seurat <- readRDS(con)
    
    # Close the connection to be able to overwrite
    close(con)
    
    
    ## EXTRACT SEURAT SPOT DATA ----
    
    
    # Seurat image object
    image_structure <- data_seurat@images$slice1
    
    # Seurat spot center coordinates in full resolution: Y = imagerow, X = imagecol
    X <- image_structure@coordinates[["imagecol"]]
    Y <- image_structure@coordinates[["imagerow"]]
    
    # Spot diameter at full resolution
    spot_radius <- image_structure@scale.factors$spot_diameter_fullres / 2
    
    
    ## ADD PATHOLOGY DATA TO SEURAT OBJECT ----
    
    
    # Load the annotations as a list of owin-s
    polys <- halo_annot_to_poly_my(path_file, region_list = expected_classes)
    pathology_classes_init <- sort(names(polys))
    
    print(pathology_classes_init)
    
    # Remove certain classes from consideration entirely: legacy classes or areas of tissue we want to filter out
    # pathology_classes_init <- setdiff(pathology_classes_init, c(reference_class, remove_class_completely))
    pathology_classes_init <- intersect(expected_classes, pathology_classes_init)
    
    # Rotate image by the user-defined amount declared in the input meta data    
    if(!is.na(pi_coef)){
      polys <- lapply(polys, function(x) spatstat::rotate(x, angle=pi*pi_coef, centre="centroid"))
    }
    
    # Compute overlap of each spot with each pathology annotation region
    spot_data_list <- list()
    for(cl in pathology_classes_init){
      spot_data_list[[cl]] <- spot_annotation_overlap_my(X, Y, spot_radius, polys[[cl]])
    }
    
    # Add annotations to the Seurat object
    for(cl in pathology_classes_init){
      col_name <- paste0(col_prefix, cl)
      data_seurat@meta.data[[col_name]] <- round(spot_data_list[[cl]])
    }
    
    # Remove certain pathology classes that shouldn't be counted towards the total
    # Data will remain in the Seurat object but won't be counted towards the Pathology_Group
    pathology_classes <- setdiff(pathology_classes_init, c(reference_class, remove_class_from_percent))
    
    # Subset pathology input to classes of interest for subsequent plotting
    polys <- polys[c("Tissue", pathology_classes)]
    
    # Sum across non-reference classes; this number will be used for subsetting to include spots that received no pathology annotation
    annotation_total_col <- paste0(col_prefix, "AnnotationTotal")
    col_names <- paste0(col_prefix, pathology_classes)
    data_seurat@meta.data[[annotation_total_col]] <- unname(rowSums(data_seurat@meta.data[col_names]))
    
    # Choose the reference for the spot surface in pixels - either from a pathology region (e.g. Tissue) or from sum of all pathology regions included in the analysis
    total_col <- paste0(col_prefix, "Total")
    if(!is.null(reference_class)){
      data_seurat@meta.data[[total_col]] <- round(spot_annotation_overlap_my(X, Y, spot_radius, reference_class))
    } else{
      data_seurat@meta.data[[total_col]] <- data_seurat@meta.data[[annotation_total_col]]
    }
    
    # Calculate the ratio of each class relative to the spot surface in pixels
    for(cl in pathology_classes){
      col_name <- paste0(col_prefix, cl, ".percent")
      data_seurat@meta.data[[col_name]] <- round(spot_data_list[[cl]] / data_seurat@meta.data[[total_col]] * 100)
    }
    
    
    ## SUBSET SEURAT OBJECT IF NECESSARY ----
    
    
    # Select spots characterized by pathology using a user-defined pixel threshold
    spots_remaining <- rownames(data_seurat@meta.data)[which(data_seurat[[annotation_total_col]] > subset_threshold)]
    
    if(do_subset == TRUE && length(spots_remaining) < ncol(data_seurat)){
      
      
      cat("Subsetting", length(spots_remaining), "spots of", ncol(data_seurat), "\n")
      
      
      # Subset Seurat object
      data_subset_seurat <- subset(data_seurat, cells = spots_remaining)
      
      # Remove any floating tissue spots using the contiguity filter
      data_subset_seurat <- tissue_contiguity_filter_my(data_subset_seurat, spot_num=50)
      
      
      # Meta data column names
      seurat_col_names <- colnames(data_subset_seurat@meta.data)
      
      # Remove old cluster data
      col_names_leave <- seurat_col_names[!grepl("SCT_snn_res", seurat_col_names)]
      data_subset_seurat@meta.data <- data_subset_seurat@meta.data[,col_names_leave]
      
      # Re-cluster data
      params <- cluster_params_my()
      # If no covariates present, default to Wilcoxon for DEA
      params[["test_use"]] <- "wilcox"  # wilcox; "LR", "negbinom", "poisson", or "MAST"
      params[["assay"]] <- "SCT"  # Spatial or SCT
      
      data_subset_seurat <- cluster_analysis_my(data_subset_seurat, params, sample_name, output_figs)
      
      # Reference cluster resolution: Table of cluster names and number of spots in each
      cluster_table_ref <- table(droplevels.data.frame(data_subset_seurat[[gsub(".0$", "", data_subset_seurat@misc$user.Clustering)]]))
      
      # Define new clustering resolution using previous clustering resolution: match the number of clusters as close as possible
      for(clust_res_loc in names(data_subset_seurat@meta.data)[grep("snn_res", names(data_subset_seurat@meta.data))]){
        
        # Table of cluster names and number of spots in each
        cluster_table <- table(droplevels.data.frame(data_subset_seurat[[clust_res_loc]]))
        clust_res <- clust_res_loc
        
        # Find clustering resolution that finds the same number of clusters as the pathology
        if(length(cluster_table) >= length(cluster_table_ref)){
          break
        }
      }
      
      # Add pathology-based clustering resolution to the Seurat object
      data_subset_seurat@misc$user.Clustering <- clust_res
      
    } else{
      
      # Keep the Seurat object as is
      data_subset_seurat <- data_seurat
      
      # Preferred clustering resolution
      clust_res <- gsub(".0$", "", data_subset_seurat@misc$user.Clustering)
      
    }
    
    
    ## ADD PATHOLOGY DOMINANT CLASS TO SEURAT OBJECT ----
    
    
    # Seurat meta-data column names storing pathology class percentages
    path_cols <- paste0(col_prefix, pathology_classes, ".percent")
    
    # Seurat meta data
    meta_seurat_df <- tibble::as_tibble(data_subset_seurat@meta.data, rownames="Coordinate")
    
    # Add a column representing "no pathology classification found for a given spot"
    col_extra <- paste0(col_prefix, "None.percent")
    meta_seurat_df[[col_extra]] <- 100 - unname(rowSums(meta_seurat_df[path_cols]))
    meta_seurat_df[which(meta_seurat_df[[col_extra]] < 0), col_extra] <- 0
    
    # Extract pathology classes and Seurat clustering data
    clust_df <- meta_seurat_df %>%
      dplyr::select(dplyr::all_of(c("Coordinate", clust_res, path_cols, col_extra))) %>%
      df_wide2long_my(key="Pathology_Group", val="Percent", start_col=3) %>%
      dplyr::mutate(Pathology_Group = gsub(paste0(col_prefix, "|.percent"), "", Pathology_Group))
    
    # Find the most abundant pathology class for each spot
    clust_max_df <- clust_df %>%
      dplyr::group_by(Coordinate) %>%
      dplyr::arrange(desc(Percent)) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(Pathology_Group, Percent)
    
    # Add the dominant cluster tag to the meta data tibble
    meta_seurat_df <- meta_seurat_df %>%
      dplyr::left_join(clust_max_df %>%
                         dplyr::select(Coordinate, Pathology_Group), by="Coordinate") %>%
      tibble::column_to_rownames("Coordinate") #%>%
    #tidyr::replace_na(list(Pathology_Group = "None"))
    
    # Add the dominant cluster tag to the Seurat object
    data_subset_seurat@meta.data <- meta_seurat_df
    colnames(data_subset_seurat@meta.data) <- gsub("Pathology_Group", paste0(col_prefix, "Group"), 
                                                   colnames(data_subset_seurat@meta.data))
    
    
    ## AGGREGATE PATHOLOGY DATA ----
    
    
    # Factorize coordinates to keep them in a specific order when plotting
    clust_df$Coordinate <- factor(clust_df$Coordinate, levels=clust_max_df$Coordinate)
    
    # Aggregate pathology classes by cluster
    clust_sum_df <- clust_df %>%
      dplyr::group_by(Pathology_Group, !!rlang::sym(clust_res)) %>%
      dplyr::summarise(Percent_Sum = sum(Percent)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(!!rlang::sym(clust_res)) %>%
      dplyr::mutate(Cluster_Sum = sum(Percent_Sum)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Percent = round(Percent_Sum / Cluster_Sum * 100))
    
    
    ## COLOR SCHEMA ----
    
    
    # Establish a color schema
    class_list <- pathology_classes
    
    if(is.null(cols)){
      num_classes <- length(pathology_classes)
      cols <- define_cols_my(n=num_classes, col_type="jet")
      
      # Establish either user-defined or alphabetic order
      if(!is.null(class_order)){
        names(cols) <- class_order
      } else{
        names(cols) <- sort(class_list)      
      }
    }
    
    # Add special classes: Tissue and None
    cols[["Tissue"]] <- "black"
    cols[["None"]] <- "white"
    
    
    ## PLOT INITIAL PATHOLOGY VIEWS ----
    
    
    # Pathology classification image reflecting full resolution of the annotation
    # Visually check that the classification and spots make sense
    filename <- paste0(output_figs, "pathology_", sample_name, "_spatial_pathology_overlay.png")
    plot_pathology_ingestion_check_my(polys, X, Y, spot_radius, filename, main_class="Tissue", title=sample_name, cols=cols, class_list=class_list)
    
    if(do_subset == TRUE){
      # Plot spot coverage by pathology assignments
      p1 <- spatial_feature_plot_my(data_subset_seurat, feature=annotation_total_col, min.cutoff="q0", max.cutoff="q100", name="")
      
      # Plot de novo clusters
      p2 <- spatial_dim_plot_my(data_subset_seurat, group.by=clust_res, title=paste0(sample_name, "\n", clust_res))
      
      # Plot de novo clusters of spots with pathology annotation
      p3 <- spatial_dim_plot_my(data_subset_seurat, group.by=clust_res, title=paste0(sample_name, "\n", clust_res))
      
      filename <- paste0(output_figs, "pathology_", sample_name, "_spatial_pathology_coverage")
      write_plot2file_my(patchwork::wrap_plots(list(p1, p2, p3), nrow=1), filename, num_row=1, num_col=3)
    }
    
    
    ## PLOT PATHOLOGY CLASSES INDEPENDENTLY ----
    
    
    # Spatial representation of each pathology layer by spot
    filename <- paste0(output_figs, "pathology_", sample_name, "_individual_classes")
    p <- batch_spatial_feature_plot_my(list(sample_name = data_seurat), paste0(col_prefix, pathology_classes_init), 
                                       output_file=filename, min.cutoff="q0", max.cutoff="q100")
    
    
    ## PLOT CLUSTERS VS PATHOLOGY ----
    
    
    # Plot pathology classes by spot in Seurat clusters
    p <- create_bar_plot_my(clust_df, x_label="Coordinate", y_label="Percent", fill_label="Pathology_Group", 
                            facet_var=c(clust_res, "free"), position="stack", filename=NULL, 
                            labels=c("Spot", "Pathology Class, %", "")) +
      scale_fill_manual(values=cols[tidyr::drop_na(clust_df) %>% 
                                      dplyr::pull(Pathology_Group) %>% 
                                      unique %>%
                                      intersect(names(cols))])
    
    filename <- paste0(output_figs, "pathology_", sample_name, "_bar_spot_clust_by_path")
    write_plot2file_my(p, filename)
    
    
    # Plot pathology classes by Seurat clusters
    p <- create_bar_plot_my(clust_sum_df, x_label=clust_res, y_label="Percent", fill_label="Pathology_Group", 
                            position="stack", filename=NULL, 
                            labels=c("Cluster", "Pathology Class, %", "")) +
      scale_fill_manual(values=cols[tidyr::drop_na(clust_df) %>% 
                                      dplyr::pull(Pathology_Group) %>% 
                                      unique %>%
                                      intersect(names(cols))])
    
    filename <- paste0(output_figs, "pathology_", sample_name, "_bar_clust_by_path")
    write_plot2file_my(p, filename)
    
    
    ## OVERWRITE SEURAT FILE ----
    
    
    if(no_rds_output == FALSE){    
      # Write the updated Seurat object
      filename <- paste0(output_path, sample_name, "_annotated_pathology.rds")
      saveRDS(data_subset_seurat, file=filename)
    }
  } else{
    warning("More than one file found!")
  }
# }
