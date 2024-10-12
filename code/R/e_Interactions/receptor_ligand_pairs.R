# Author: Anna Lyubetskaya. Date: 21-03-11

# This script measures co-expression of receptor-ligand gene pairs in ST data
# The source of RL pairs: FANTOM5
# Co-expression is measured by calculating Pearson calculation and mutual information between scaled normalized gene profiles


## SETUP ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_vis.R")


## PARAMETERS ----


# Run name
run_name <- "FANTOM5"

# Cohort name
cohort_name <- "PDAC108_path14_harmonyepi_rpca_sct"
# Select a specific sample within a cohort, e.g. HumanPanc_ROI1_FFPE_A_Apr21
cohort_subset <- NULL

# SCT threshold for abundant genes
sct_threshold <- 1
# Spot threshold for abundant genes
spot_threshold <- 1000

# Seurat assay and slot to extract
assay <- "SCT"
slot <- "scale.data"

# Subset to specific clusters for modeling if needed and instill that order for vis
cluster_order <- c(3, 10, 0, 4, 6, 9)


## PATHS ----


# Location of annotated Seurat objects
input_path <- paste0("XXXX")

# Receptor-ligand pairs
receptor_ligand_files <- c("XXXX")

# Output information
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_name, "_", run_name, "/")

# Create output folder
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Read all files with RL pairs
rl_list <- list()
for(f in receptor_ligand_files){
  rl_list[[f]] <- readr::read_delim(f, delim="\t") %>%
    dplyr::select(Pair.Name, Ligand.ApprovedSymbol, Receptor.ApprovedSymbol)
}

# Read in receptor-ligand pairs
rl_df <- dplyr::bind_rows(rl_list)

# Read the Seurat object
data_seurat <- readRDS(input_path)


# Find user defined resolution
resolution <- data_seurat@misc$user.Clustering

# Remove clusters if necessary
if(!is.null(cluster_order)){
  barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data[[resolution]] %in% cluster_order),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# Subset an integrated cohort to a specific sample if necessary
if(!is.null(cohort_subset)){
  data_seurat <- subset(data_seurat, cells=rownames(data_seurat@meta.data %>% 
                                                      dplyr::filter(user.Sample_Name == cohort_subset)))
  
  data_seurat@images <- data_seurat@images[cohort_subset]
}


## WRANGLE SEURAT DATA ----


# List of abundant genes
gene_list <- seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, 
                                             spot_threshold=spot_threshold, assay="SCT", slot="data")

# Extract meta data
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")


# If user didn't define the cluster order, just use the default cluster order
if(is.null(cluster_order)){
  cluster_order <- unique(sort(meta_df[[resolution]]))
}


# Expression data to long tibble
data_df <- seurat_expression_to_long_tibble_my(data_seurat, assay=assay, slot=slot) %>%
  dplyr::filter(Symbol %in% c(rl_df$Ligand.ApprovedSymbol, rl_df$Receptor.ApprovedSymbol) &
                  Symbol %in% gene_list) %>%
  dplyr::inner_join(meta_df %>% 
                      dplyr::select(Coordinate, !!rlang::sym(resolution)),
                    by="Coordinate")

cat("Found", length(unique(data_df$Symbol)), "genes of", 
    length(unique(c(rl_df$Ligand.ApprovedSymbol, rl_df$Receptor.ApprovedSymbol))))

# Calculate mean expression of each gene of interest in each cluster
data_clust_df <- data_df %>% 
  dplyr::group_by(!!rlang::sym(resolution), Symbol) %>% 
  dplyr::summarise(ScoreMean = mean(SCT_scale.data), 
                   ScoreSD = sd(SCT_scale.data)) %>% 
  dplyr::ungroup()

#  Write to file
filename <- paste0(output_path, cohort_name, "_score_by_clust.txt")
readr::write_delim(data_clust_df, filename, delim="\t")

# Find the cluster with the highest expression of the gene
data_clust_max_df <- data_clust_df %>% 
  dplyr::group_by(Symbol) %>% 
  dplyr::top_n(1, ScoreMean)


## PROCESS RECEPTOR-LIGAND PAIRS BY SPOT ----


# Map each receptor-ligand pair to filtered ST expression data
rl_spot_df <- rl_df %>%
  dplyr::select(Receptor.ApprovedSymbol, Ligand.ApprovedSymbol) %>%
  dplyr::inner_join(data_df, by=c("Receptor.ApprovedSymbol" = "Symbol")) %>%
  dplyr::inner_join(data_df, by=c("Ligand.ApprovedSymbol" = "Symbol", "Coordinate" = "Coordinate")) %>%
  dplyr::rename(Receptor.SCT = SCT_scale.data.x, Ligand.SCT = SCT_scale.data.y) %>%
  dplyr::mutate(SampleName = gsub(":.+$", "", Coordinate))

# A list of receptor-ligand pairs with robust expression data
rl_pairs <- rl_spot_df %>% 
  dplyr::select(Receptor.ApprovedSymbol, Ligand.ApprovedSymbol) %>% 
  unique()

# Sample names
sample_list <- unique(rl_spot_df$SampleName)

cat("Number of pairs =", nrow(rl_pairs))


## CALCULATE DISTANCE BETWEEN RECEPTOR AND LIGAND ----


dist_list <- list()

# Calculate correlation and mutual information between receptor and ligand signals
for(i in 1:nrow(rl_pairs)){
  # Extract spot SCT data for the pair
  data_loc_df <- rl_spot_df %>% 
    dplyr::filter(Receptor.ApprovedSymbol == rl_pairs[[i, 1]] & 
                    Ligand.ApprovedSymbol == rl_pairs[[i, 2]])
  
  # Find receptor and ligand thresholds
  r_mean <- mean(data_loc_df$Receptor.SCT)
  l_mean <- mean(data_loc_df$Ligand.SCT)
  
  # Find status of receptor and ligand in each spot
  data_loc_df <- data_loc_df %>% 
    dplyr::mutate(Receptor.Status = Receptor.SCT > r_mean,
                  Ligand.Status = Ligand.SCT > l_mean)
  
  # Status counts for a receptor-ligand pair
  rl_table <- table(data_loc_df[c("Receptor.Status", "Ligand.Status")])
  
  # Save all distance information to a list
  dist_list[[i]] <- c(rl_pairs[[i, 1]], rl_pairs[[i, 2]],
                      round(cor(data_loc_df$Receptor.SCT, data_loc_df$Ligand.SCT, method="pearson"), 3),
                      round(entropy::mi.plugin(rl_table), 3),
                      rl_table[1,1], rl_table[2,1], rl_table[1,2], rl_table[2,2])
}

# Transform list of lists to a tibble
dist_df <- tibble::as_tibble(as.data.frame(do.call(rbind, dist_list)), .name_repair="unique") %>%
  dplyr::rename(Receptor = V1, Ligand = V2, R2 = V3, MI = V4, 
                NONE = V5, ReceptorOnly = V6, LigandOnly = V7, ReceptorLigand = V8) %>%
  dplyr::mutate(R2 = as.numeric(R2), MI = as.numeric(MI),
                NONE = as.numeric(NONE), ReceptorOnly = as.numeric(ReceptorOnly),
                LigandOnly = as.numeric(LigandOnly), ReceptorLigand = as.numeric(ReceptorLigand)) %>%
  dplyr::arrange(desc(R2), desc(MI)) %>%
  dplyr::inner_join(rl_df, by=c("Receptor" = "Receptor.ApprovedSymbol", "Ligand" = "Ligand.ApprovedSymbol")) %>%
  dplyr::left_join(data_clust_max_df, by=c("Receptor" = "Symbol")) %>%
  dplyr::rename(ReceptorClusterMax = !!rlang::sym(resolution), ReceptorMaxMean = ScoreMean, ReceptorMaxSD = ScoreSD) %>%
  dplyr::left_join(data_clust_max_df, by=c("Ligand" = "Symbol")) %>%
  dplyr::rename(LigandClusterMax = !!rlang::sym(resolution), LigandMaxMean = ScoreMean, LigandMaxSD = ScoreSD)

#  Write to file
filename <- paste0(output_path, cohort_name, "_rl_dist.txt")
readr::write_delim(dist_df, filename, delim="\t")

# Establish distance threshold - correlation or MI
cor_threshold <- mean(dist_df$R2) + sd(dist_df$R2) * 2
mi_threshold <- mean(dist_df$MI) + sd(dist_df$MI) * 2

# Select RL pairs that pass the threshold
dist_df <- dist_df %>%
  dplyr::mutate(Associated = R2 >= cor_threshold | MI >= mi_threshold)

# Plot distance distributions: correlation
p1 <- create_hist_plot_my(dist_df, x_label="R2", fill_label="Associated", 
                          intercept=c(-1, 0, cor_threshold, 1), binwidth=0.05, filename=NULL, labels=c("R2", "Number of RL pairs", ""))

# Plot distance distributions: mutual information
p2 <- create_hist_plot_my(dist_df, x_label="MI", fill_label="Associated", 
                          intercept=c(0, mi_threshold), binwidth=0.005, filename=NULL, labels=c("R2", "Number of RL pairs", ""))

filename <- paste0(output_path, cohort_name, "_dist_hist")
write_plot2file_my(patchwork::wrap_plots(list(p1, p2), nrow=2), filename, num_row=2)

# Filter pairs by their expression distance - correlation or MI
dist_loc_df <- dist_df %>%
  dplyr::filter(Associated == TRUE)

#  Write to file
filename <- paste0(output_path, cohort_name, "_rl_dist_filt.txt")
readr::write_delim(dist_df, filename, delim="\t")


## VISUALIZE DATA ----


data_df[[resolution]] <- factor(data_df[[resolution]], levels=cluster_order)


if(nrow(dist_loc_df) > 0){
  
  # Plot the receptor ligand pairs
  for(i in 1:nrow(dist_loc_df)){
    
    if(length(intersect(c(dist_loc_df[[i, 1]], dist_loc_df[[i, 2]]), rownames(data_seurat))) == 2){
      
      # Extract spot SCT data for the pair
      data_loc_df <- rl_spot_df %>% 
        dplyr::filter(Receptor.ApprovedSymbol == dist_loc_df[[i, 1]] & 
                        Ligand.ApprovedSymbol == dist_loc_df[[i, 2]])
      
      r_mean <- mean(data_loc_df$Receptor.SCT)
      l_mean <- mean(data_loc_df$Ligand.SCT)
      
      # Find status of receptor and ligand in each spot
      data_loc_df <- data_loc_df %>% 
        dplyr::mutate(Receptor.Status = ifelse(Receptor.SCT > r_mean, "Receptor", ""),
                      Ligand.Status = ifelse(Ligand.SCT > l_mean, "Ligand", ""),
                      Status = paste0(Receptor.Status, Ligand.Status),
                      Status = ifelse(Status == "", "None", Status))
      
      data_seurat@meta.data["Status"] <- data_loc_df[match(rownames(data_seurat@meta.data), data_loc_df$Coordinate), "Status"]
      data_seurat@meta.data["Status"] <- factor(data_seurat@meta.data$Status, 
                                                levels=c("ReceptorLigand", "Receptor", "Ligand", "None"))
      
      filename1 <- paste0(output_path, cohort_name, "_", dist_loc_df[[i, 1]], "_", dist_loc_df[[i, 2]], "_spatial")
      filename2 <- paste0(output_path, cohort_name, "_", dist_loc_df[[i, 1]], "_", dist_loc_df[[i, 2]], "_box")
      filename3 <- paste0(output_path, cohort_name, "_", dist_loc_df[[i, 1]], "_", dist_loc_df[[i, 2]], "_spatial_binary")
      
      if(!file.exists(paste0(filename1, ".png"))){
        # Define colors
        cols <- c(ReceptorLigand = "blue4", Receptor = "springgreen", Ligand = "deepskyblue", None = "grey95")
        
        # Write combo plot to file
        title <- paste(cohort_name, "\nRECEPTOR =", dist_loc_df[[i, 1]], " LIGAND =", dist_loc_df[[i, 2]],
                       "\nR2 =", dist_loc_df[[i, 3]], " MI =", dist_loc_df[[i, 4]], 
                       "\nNA =", dist_loc_df[[i, 5]], "ROnly =", dist_loc_df[[i, 6]], 
                       "LOnly =", dist_loc_df[[i, 7]], "RL =", dist_loc_df[[i, 8]])
        
        # Plot RL levels by cluster
        p <- create_box_plot_my(data_df %>% 
                                  dplyr::filter(Symbol %in% dist_loc_df[i, c(1, 2)]) %>%
                                  dplyr::mutate(Symbol = ifelse(Symbol == dist_loc_df[[i, "Receptor"]], "Receptor", Symbol),
                                                Symbol = ifelse(Symbol == dist_loc_df[[i, "Ligand"]], "Ligand", Symbol)),
                                x_label=resolution, y_label="SCT_scale.data", fill_label="Symbol",
                                filename=NULL, labels=c("Clusters", "Normalized Expression", title,
                                                        cols=cols))
        
        write_plot2file_my(p, filename2)
        
        # # Plot signature on the original tissue slice
        # p1 <- spatial_feature_plot_my(data_seurat, dist_loc_df[[i, 1]], title=dist_loc_df[[i, 1]])
        # 
        # # Plot signature on the original tissue slice
        # p2 <- spatial_feature_plot_my(data_seurat, dist_loc_df[[i, 2]], title=dist_loc_df[[i, 2]])
        
        # # Spatial binary plot of RL expression status
        # p3 <- spatial_dim_plot_my(data_seurat, group.by="Status", cols=cols)

        # write_plot2file_my(patchwork::wrap_plots(list(p1, p2), nrow=2) +
        #                      patchwork::plot_annotation(title = title), 
        #                    filename1, num_row=2, num_col=length(data_seurat@images))
        
        # write_plot2file_my(p3 + patchwork::plot_annotation(title = title), 
        #                    filename3, num_row=1, num_col=length(data_seurat@images))
      }
      
    }
  }
}
