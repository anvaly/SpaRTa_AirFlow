# Author: Anna Lyubetskaya. Date: 21-03-29
# Merge and compare DEA results


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

library(ggplot2)

source("code/utils/utils_tibble.R")
source("code/utils/utils_specialized_plots.R")


## PARAMETERS ----


# The name of the input sample / cohort
cohort_name <- "HumanPDAC_10K_sct"

# Resolution column name
resolution <- "integrated_snn_res.0.4"  # "Pathology_Group", "integrated_snn_res.0.4"

# PCT threshold for plotting
pct_list <- c(20, 30, 40, 50, 60, 70, 80, 90)
pct_list <- c(50)


## PATHS ----


# Input folder
input_path_table1 <- "XXXX"
input_path_table2 <- "XXXX"

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_name, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Annotated DE genes, file 1
dea1_df <- readr::read_delim(input_path_table1, delim="\t")

# Annotated DE genes, file 2
dea2_df <- readr::read_delim(input_path_table2, delim="\t")

# Gene list shared between both DE gene lists
gene_list <- unique(intersect(dea1_df$Symbol, dea2_df$Symbol))

# List of clusters in the first file
clust1_list <- as.character(unique(dea1_df$cluster))
# List of clusters in the second file
clust2_list <- as.character(unique(dea2_df$cluster))


for(pct_min in pct_list){
  
  ## WRANGLE DATA ----
  
  
  # First set of clusters and markers
  clust1_df <- dea1_df %>%
    dplyr::filter(pct_1 >= pct_min) %>%
    dplyr::filter(Symbol %in% gene_list & direction == "UP") %>%
    dplyr::select(Symbol, cluster) %>%
    unique
  
  # Second set of clusters and markers
  clust2_df <- dea2_df %>%
    dplyr::filter(pct_1 >= pct_min) %>%
    dplyr::filter(Symbol %in% gene_list & direction == "UP") %>%
    dplyr::select(Symbol, cluster) %>%
    unique
  
  
  ## CALCULATE INTERSECTION ----
  
  
  intersect_list <- list()
  
  # Go through every pair of clusters
  for(c1 in unique(clust1_df$cluster)){
    for(c2 in unique(clust2_df$cluster)){
      # Pair name
      name <- paste(c1, c2)
      
      # Genes in the first cluster from the first file
      gene1_list <- clust1_df %>%
        dplyr::filter(cluster == c1) %>%
        dplyr::pull(Symbol)
      
      # Genes in the second cluster from the second file
      gene2_list <- clust2_df %>%
        dplyr::filter(cluster == c2) %>%
        dplyr::pull(Symbol)
      
      # Intersect gene lists between two clusters
      intersect_count <- length(intersect(gene1_list, gene2_list))
      group_min <- min(length(gene1_list), length(gene2_list))
      intersect_perc <- round(intersect_count / group_min * 100)
      
      # Store data in a list
      intersect_list[[name]] <- c(c1, c2, intersect_count, group_min, intersect_perc)
    }
  }
  
  # Tranform list of lists into a long tibble
  dist_df <- tibble::as_tibble(as.data.frame(do.call(rbind, intersect_list)), .name_repair="unique") %>%
    dplyr::rename(cluster1 = V1, cluster2 = V2, IntersectCount = V3, GroupMin = V4, IntersectPercent = V5) %>%
    dplyr::mutate(cluster1 = as.character(cluster1),
                  cluster2 = as.character(cluster2))
  
  
  ## FIND FIT ----
  
  
  # A pair of clusters with the biggest marker intersection
  dist_fit_df <- dist_df %>%
    dplyr::filter(IntersectPercent > 0) %>%
    dplyr::group_by(cluster1) %>%
    dplyr::top_n(1, IntersectPercent)
  
  # Clusters that got lost in the intersection
  clust_missing1 <- setdiff(clust1_list, sort(unique(dist_fit_df$cluster1)))
  clust_missing2 <- setdiff(clust2_list, sort(unique(dist_fit_df$cluster2)))
  
  cat("Cluster:", pct_min, "missing from 1st file:", clust_missing1, "\n")
  cat("Cluster:", pct_min, "missing from 2nd file:", clust_missing2, "\n\n")
  
  
  ## CALCULATE INTERSECTION ----
  
  
  # Distribution stats
  dist_min <- min(dist_df$IntersectPercent)
  dist_max <- max(dist_df$IntersectPercent)
  dist_mean <- round(mean(dist_df$IntersectPercent))
  dist_sd <- round(sd(dist_df$IntersectPercent))
  
}


# Plot the heatmap of the last intersection
ggplot(dist_df, aes(cluster1, cluster2, fill=IntersectPercent)) + 
  geom_tile() +
  scale_fill_gradient2(low="red", mid="white", high="blue", midpoint=dist_mean)
