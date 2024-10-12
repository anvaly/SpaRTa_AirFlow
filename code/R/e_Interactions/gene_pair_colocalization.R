# Author: Anna Lyubetskaya. Date: 21-11-04
# Measure gene pair co-localization


## SETUP ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_specialized_plots.R")

source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_vis.R")


## PARAMETERS ----


# Cohort name
cohort_name <- "Syng_B16-MC38_cca_sct"
# Select a specific sample within a cohort, e.g. HumanPanc_ROI1_FFPE_A_Apr21
cohort_subset <- NULL

# SCT threshold for abundant genes
sct_threshold <- 1
# Spot threshold for abundant genes
spot_threshold <- 10

# Seurat assay and slot to extract
assay <- "SCT"
slot <- "scale.data"

# Create gene pair individual plots
do_individual_plots <- FALSE

# Gene list
gene_interest_list <- c("PDCD1", "CD274", # "LAG3"
                        "TAP2", "TAP1", "B2M", 
                        "STAT1", #"JAK1", 
                        "CXCL9", "CXCL10", "CXCL11", #"IFNGR1",
                        "CD8A", "CD8B1", "PRF1", "GZMB", 
                        #"FOXP3",
                        "CD68", "CSF1R", 
                        #"XCR1", "ITGAE", "BATF3",
                        "COL1A1", "COL1A2", "ACTA2", "VIM", #"FAP",
                        "MIF",
                        "SHISAL2B", "RHOX5",
                        "TYRP1", "DCT", "PMEL", "MLANA"
)

gene_categories <- c("Other", "Other", # "Other",
                     "AP", "AP", "AP",
                     "IFN", #"Other",
                     "IFN", "IFN", "IFN", #"IFNg",
                     "CD8", "CD8", "CD8", "CD8",
                     #"T-reg",
                     "Mac", "Mac",
                     #"DC", "DC", "DC",
                     "Fib", "Fib", "Fib", "Fib", #"Fib",
                     "Other",
                     "Tumor", "Tumor",
                     "Tumor", "Tumor", "Tumor", "Tumor"
)

# Create a gene annotation tibble
gene_interest_df <- tibble::tibble(Symbol = gene_interest_list,Category = gene_categories)
gene_interest_df["Symbol"] <- factor(gene_interest_df[["Symbol"]], levels=gene_interest_list)


## PATHS ----


# Location of annotated Seurat objects
input_path <- paste0("XXXX")
#input_path <- paste0("XXXX")

# Output information
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_name, "/")

# Create output folder
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Read the Seurat object
data_seurat <- readRDS(input_path)

# Subset an integrated cohort to a specific sample if necessary
if(!is.null(cohort_subset)){
  data_seurat <- subset(data_seurat, cells=rownames(data_seurat@meta.data %>% 
                                                      dplyr::filter(user.Sample_Name == cohort_subset)))
}


## WRANGLE SEURAT DATA ----


# Find user defined resolution
resolution <- data_seurat@misc$user.Clustering

# List of abundant genes
gene_list <- seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, 
                                             spot_threshold=spot_threshold, assay="SCT", slot="data")

# Extract meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

# Expression data in non-scaled SCT values
data_sct_df <- seurat_expression_to_long_tibble_my(data_seurat, assay=assay, slot="data") %>%
  dplyr::filter(Symbol %in% gene_interest_list)

# Add levels to symbols for future plotting order
data_sct_df["Symbol"] <- factor(data_sct_df[["Symbol"]], levels=gene_interest_list)

# Expression data to long tibble
data_df <- seurat_expression_to_long_tibble_my(data_seurat, assay=assay, slot=slot) %>%
  dplyr::filter(Symbol %in% gene_interest_list) %>%
  dplyr::inner_join(meta_df %>% 
                      dplyr::select(Coordinate, !!rlang::sym(resolution), user.Sample_Name),
                    by="Coordinate")

cat("Found", length(unique(data_df$Symbol)), "genes of", length(gene_interest_list))
setdiff(gene_interest_list, unique(data_df$Symbol))

# Calculate mean expression of each gene of interest in each cluster
data_clust_df <- data_df %>% 
  dplyr::group_by(!!rlang::sym(resolution), Symbol) %>% 
  dplyr::summarise(ScoreMean = mean(SCT_scale.data), ScoreSD = sd(SCT_scale.data)) %>% 
  dplyr::ungroup()

#  Write to file
filename <- paste0(output_path, cohort_name, "_score_by_clust.txt")
readr::write_delim(data_clust_df, filename, delim="\t")

# Find the cluster with the highest expression of the gene
data_clust_max_df <- data_clust_df %>% 
  dplyr::group_by(Symbol) %>% 
  dplyr::top_n(1, ScoreMean)


## CALCULATE DISTANCE BETWEEN GENE PAIRS ----


dist_list <- list()

# Calculate correlation and mutual information between receptor and ligand signals
for(gene_x in gene_interest_list){
  for(gene_y in gene_interest_list){
    
    # Extract spot SCT data for the pair
    data_loc_df <- data_df %>%
      dplyr::filter(Symbol == gene_x) %>%
      dplyr::inner_join(data_df %>%
                          dplyr::filter(Symbol == gene_y), 
                        by=c("Coordinate", resolution, "user.Sample_Name"))
    
    # Find receptor and ligand thresholds
    x_mean <- mean(data_loc_df$SCT_scale.data.x)
    y_mean <- mean(data_loc_df$SCT_scale.data.y)
    
    cat("Gene means:", x_mean, y_mean, "\n")
    
    # Find status of receptor and ligand in each spot
    data_loc_df <- data_loc_df %>% 
      dplyr::mutate(GeneX.Status = SCT_scale.data.x > x_mean,
                    GeneY.Status = SCT_scale.data.y > y_mean)
    
    # Status counts for a receptor-ligand pair
    pair_table <- table(data_loc_df[c("GeneX.Status", "GeneY.Status")])
    
    # Save all distance information to a list
    if(ncol(pair_table) == 2 && nrow(pair_table) == 2){
      dist_list[[paste(gene_x, gene_y)]] <- c(gene_x, gene_y,
                                              round(cor(data_loc_df$SCT_scale.data.x, data_loc_df$SCT_scale.data.y, method="pearson"), 3),
                                              round(entropy::mi.plugin(pair_table), 3),
                                              pair_table[1,1], pair_table[2,1], pair_table[1,2], pair_table[2,2])
    }
  }
}

# Tranform list of lists to a tibble
dist_df <- tibble::as_tibble(as.data.frame(do.call(rbind, dist_list)), .name_repair="unique") %>%
  dplyr::rename(GeneX = V1, GeneY = V2, R2 = V3, MI = V4, 
                NONE = V5, GeneXOnly = V6, GeneYOnly = V7, Both = V8) %>%
  dplyr::mutate(R2 = as.numeric(R2), MI = as.numeric(MI),
                NONE = as.numeric(NONE), GeneXOnly = as.numeric(GeneXOnly), GeneYOnly = as.numeric(GeneYOnly), Both = as.numeric(Both),
                PercentOccupancy = round(Both / (GeneXOnly + GeneYOnly + Both) * 100)) %>%
  dplyr::arrange(desc(R2), desc(MI)) %>%
  dplyr::left_join(data_clust_max_df, by=c("GeneX" = "Symbol")) %>%
  dplyr::rename(GeneXClusterMax = !!rlang::sym(resolution), GeneXMaxMean = ScoreMean, GeneXMaxSD = ScoreSD) %>%
  dplyr::left_join(data_clust_max_df, by=c("GeneY" = "Symbol")) %>%
  dplyr::rename(GeneYClusterMax = !!rlang::sym(resolution), GeneYMaxMean = ScoreMean, GeneYMaxSD = ScoreSD) %>%
  dplyr::left_join(gene_interest_df, by=c("GeneX" = "Symbol")) %>%
  dplyr::left_join(gene_interest_df, by=c("GeneY" = "Symbol"))

#  Write to file
filename <- paste0(output_path, cohort_name, "_pair_dist.txt")
readr::write_delim(dist_df, filename, delim="\t")

# Plot distance distributions
hist(dist_df$R2)
hist(dist_df$MI)
hist(dist_df$PercentOccupancy)

# Establish distance threshold - correlation or MI
cor_threshold <- mean(dist_df$R2) + sd(dist_df$R2)
mi_threshold <- mean(dist_df$MI) + sd(dist_df$MI)

# Filter pairs by their expression distance - correlation or MI
dist_loc_df <- dist_df %>%
  dplyr::filter((R2 >= cor_threshold | MI >= mi_threshold) &
                  GeneX != GeneY)

#  Write to file
filename <- paste0(output_path, cohort_name, "_pair_dist_filt.txt")
readr::write_delim(dist_df, filename, delim="\t")


## HEATMAPS ----


# Make a boxplot showing gene expression levels for reference

p <- create_box_plot_my(data_sct_df %>%
                     dplyr::inner_join(gene_interest_df, by="Symbol"), 
                   x_label="Symbol", y_label="SCT_data", fill_label="Category", filename=NULL, 
                   labels=c("", "Expression, SCT", "")) +
scale_fill_manual(values=define_cols_my(n=length(unique(gene_categories)), col_type="jet"))

# Write to file
filename <- paste0(output_path, cohort_name, "_box_sct")
write_plot2file_my(p, filename)

# Setup heatmap parameters
params <- list(
  cell_value = "R2",
  row_label = "GeneX", 
  col_label = "GeneY", 
  distance = "euclidean",
  row_annotation = c("Category.x", "GeneXClusterMax"),
  col_annotation = NULL,  # c("Category.y", "GeneYClusterMax"),
  range = c(-0.25, 0, 0.25),
  colors = c("darkblue", "white", "darkred"),
  column_order = intersect(gene_interest_list, dist_df$GeneY),
  row_order = intersect(gene_interest_list, dist_df$GeneX),
  cluster_rows=FALSE, cluster_cols=FALSE)

# Create a heatmap of DE biomarkers
filename <- paste0(output_path, cohort_name, "_hm_corr.png")
# Plot the heatmap
hm <- create_heatmap_my(dist_df, params, row_meta_df=dist_df, col_meta_df=dist_df, filename=filename, width=8, height=8)


# Setup heatmap parameters
params <- list(
  cell_value = "MI",
  row_label = "GeneX", 
  col_label = "GeneY", 
  distance = "euclidean",
  row_annotation = c("Category.x", "GeneXClusterMax"),
  col_annotation = NULL,  # c("Category.y", "GeneYClusterMax"),
  range = c(0, 0.1, 0.5),
  colors = c("white", "pink", "darkred"),
  column_order = intersect(gene_interest_list, dist_df$GeneY),
  row_order = intersect(gene_interest_list, dist_df$GeneX),
  cluster_rows=FALSE, cluster_cols=FALSE)

# Create a heatmap of DE biomarkers
filename <- paste0(output_path, cohort_name, "_hm_mi.png")
# Plot the heatmap
hm <- create_heatmap_my(dist_df, params, row_meta_df=dist_df, col_meta_df=dist_df, filename=filename, width=10, height=10)


# Setup heatmap parameters
params <- list(
  cell_value = "PercentOccupancy",
  row_label = "GeneX", 
  col_label = "GeneY", 
  distance = "euclidean",
  row_annotation = c("Category.x", "GeneXClusterMax"),
  col_annotation = NULL,  # c("Category.y", "GeneYClusterMax"),
  range = c(0, 20, 50),
  colors = c("white", "pink", "darkred"),
  column_order = intersect(gene_interest_list, dist_df$GeneY),
  row_order = intersect(gene_interest_list, dist_df$GeneX),
  cluster_rows=FALSE, cluster_cols=FALSE)

# Create a heatmap of DE biomarkers
filename <- paste0(output_path, cohort_name, "_hm_occupancy.png")
# Plot the heatmap
hm <- create_heatmap_my(dist_df, params, row_meta_df=dist_df, col_meta_df=dist_df, filename=filename, width=8, height=8)


## VISUALIZE DATA ----


if(nrow(dist_loc_df) > 0 && do_individual_plots == TRUE){
  # Define colors
  cols <- c(GeneXGeneY = "blue4", GeneX = "springgreen", GeneY = "deepskyblue", None = "grey95")
  
  # Plot the highly correlated gene pairs
  for(i in 1:nrow(dist_loc_df)){
    
    #Initialize variables
    gene_x <- dist_loc_df[[i, "GeneX"]]
    gene_y <- dist_loc_df[[i, "GeneY"]]
    
    # Extract spot SCT data for the gene pair
    data_loc_df <- data_df %>%
      dplyr::filter(Symbol == gene_x) %>%
      dplyr::inner_join(data_df %>%
                          dplyr::filter(Symbol == gene_y), 
                        by=c("Coordinate", resolution, "user.Sample_Name"))
    
    # Find gene specific thresholds
    x_mean <- mean(data_loc_df$SCT_scale.data.x)
    y_mean <- mean(data_loc_df$SCT_scale.data.y)
    
    cat("Gene means:", x_mean, y_mean, "\n")
    
    # Find status of both genes in each spot
    data_loc_df <- data_loc_df %>% 
      dplyr::mutate(GeneX.Status = ifelse(SCT_scale.data.x > x_mean, "GeneX", ""),
                    GeneY.Status = ifelse(SCT_scale.data.y > y_mean, "GeneY", ""),
                    Status = paste0(GeneX.Status, GeneY.Status),
                    Status = ifelse(Status == "", "None", Status))
    
    # Set level order for the status
    data_seurat@meta.data["Status"] <- data_loc_df[match(rownames(data_seurat@meta.data), data_loc_df$Coordinate), "Status"]
    data_seurat@meta.data["Status"] <- factor(data_seurat@meta.data$Status,
                                              levels=c("GeneXGeneY", "GeneX", "GeneY", "None"))
    
    
    # Initiate output files
    filename1 <- paste0(output_path, "aux_", cohort_name, "_", gene_x, "_", gene_y, "_spatial")
    filename2 <- paste0(output_path, "aux_", cohort_name, "_", gene_x, "_", gene_y, "_box")
    filename3 <- paste0(output_path, "aux_", cohort_name, "_", gene_x, "_", gene_y, "_spatial_binary")
    
    # Plot RL levels by cluster
    p <- create_box_plot_my(data_df %>%
                              dplyr::filter(Symbol %in% c(gene_x, gene_y)) %>%
                              dplyr::mutate(Symbol = ifelse(Symbol == gene_x, "GeneX", Symbol),
                                            Symbol = ifelse(Symbol == gene_y, "GeneY", Symbol)),
                            x_label=resolution, y_label="SCT_scale.data", fill_label="Symbol",
                            filename=NULL, labels=c("Clusters", "Normalized Expression",
                                                    paste("GeneX =", gene_x, "GeneY =", gene_y))) +
      scale_fill_manual(values=cols)
    
    # Write to file
    write_plot2file_my(p, filename2)
    
    # Plot signature on the original tissue slice
    p1 <- spatial_feature_plot_my(data_seurat, gene_x, name=gene_x)
    
    # Plot signature on the original tissue slice
    p2 <- spatial_feature_plot_my(data_seurat, gene_y, name=gene_y)
    
    # Spatial binary plot of RL expression status
    p3 <- spatial_dim_plot_my(data_seurat, group.by="Status", cols=cols)
    
    # Write combo plot to file
    title <- paste(cohort_name, "\nGeneX =", gene_x, " GeneY =", gene_y,
                   "\nR2 =", dist_loc_df[[i, "R2"]], " MI =", dist_loc_df[[i, "MI"]],
                   "\nNA =", dist_loc_df[[i, "None"]], "XOnly =", dist_loc_df[[i, "GeneXOnly"]],
                   "YOnly =", dist_loc_df[[i, "GeneXOnly"]], "Both =", dist_loc_df[[i, "Both"]])
    
    write_plot2file_my(patchwork::wrap_plots(list(p1, p2), nrow=2) +
                         patchwork::plot_annotation(title = title),
                       filename1, num_row=2, num_col=length(data_seurat@images))
    
    write_plot2file_my(p3 + patchwork::plot_annotation(title = title),
                       filename3, num_row=1, num_col=length(data_seurat@images))
    
  }
}
