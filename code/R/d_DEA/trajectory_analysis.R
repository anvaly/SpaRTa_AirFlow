# Author: Anna Lyubetskaya. Date: 22-08-29
# Attempt trajectory analysis


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

source("code/utils/utils_ggplot.R")


## PARAMETERS ----


# The name of the input sample / cohort
cohort_name <- "path_epi_0.2"

# Resolution column name
resolution <- "integrated_snn_res.0.2"

# Use the following clusters
# cluster_list <- c(0, 1, 2, 3, 4, 5, 6, 8)
cluster_list <- c(0, 1, 2, 3)


## PATHS ----


# Input folder
input_path_object <- "XXXX"

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, cohort_name, "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


## INGEST DATA ----


# Seurat data
data_seurat <- readRDS(input_path_object)

# Subset Seurat data to specific clusters
if(!is.null(cluster_list)){
  barcode_list <- data_seurat@meta.data %>%
    dplyr::filter(!!rlang::sym(resolution) %in% cluster_list) %>%
    rownames()
  
  data_seurat <- subset(data_seurat, cells=barcode_list)
}


## TRAJECTORY ANALYSIS ----


cds <- SeuratWrappers::as.cell_data_set(data_seurat)

# Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors(cds)

cds <- monocle3::cluster_cells(cds)

p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)

filename <- paste0(output_path, "/cluster_", cohort_name)
write_plot2file_my(patchwork::wrap_plots(list(p1, p2), nrow=1), filename, num_row=1, num_col=4)


cds <- monocle3::learn_graph(cds)

p3 <- plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

filename <- paste0(output_path, "/cluster_", cohort_name, "_trajectory")
write_plot2file_my(p3, filename, num_row=1, num_col=2)


cds <- order_cells(cds)

p4 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=FALSE, label_leaves=FALSE,
                 label_branch_points=FALSE, graph_label_size=1.5)

filename <- paste0(output_path, "/cluster_", cohort_name, "_pseudotime")
write_plot2file_my(p4, filename, num_row=1, num_col=2)


cds_sub_class <- choose_graph_segments(cds)

cds_sub_class <- preprocess_cds(cds_sub_class, num_dim = 50)
cds_sub_class <- reduce_dimension(cds_sub_class)
cds_sub_class <- monocle3::cluster_cells(cds_sub_class)
cds_sub_class <- monocle3::learn_graph(cds_sub_class)

p5 <- plot_cells(cds_sub_class, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

filename <- paste0(output_path, "/cluster_", cohort_name, "_trajectory_classical")
write_plot2file_my(p5, filename, num_row=1, num_col=2)

cds_test_res <- graph_test(cds_sub_class, neighbor_graph="principal_graph", cores=4, verbose=TRUE)
pr_deg_ids <- tibble::as_tibble(as.matrix(subset(cds_test_res, q_value < 0.05)), rownames="Coordinate")

filename <- paste0(output_path, "/cluster_", cohort_name, "_trajectory_classical.txt")
readr::write_delim(pr_deg_ids, filename, delim="\t")


cds_sub_basal1 <- choose_graph_segments(cds)

cds_sub_basal1 <- preprocess_cds(cds_sub_basal1, num_dim = 50)
cds_sub_basal1 <- reduce_dimension(cds_sub_basal1)
cds_sub_basal1 <- monocle3::cluster_cells(cds_sub_basal1)
cds_sub_basal1 <- monocle3::learn_graph(cds_sub_basal1)

p5 <- plot_cells(cds_sub_basal1, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

filename <- paste0(output_path, "/cluster_", cohort_name, "_trajectory_basal3")
write_plot2file_my(p5, filename, num_row=1, num_col=2)

cds_test_res <- graph_test(cds_sub_basal1, neighbor_graph="principal_graph", cores=4, verbose=TRUE)
pr_deg_ids <- tibble::as_tibble(as.matrix(subset(cds_test_res, q_value < 0.05)), rownames="Coordinate")

filename <- paste0(output_path, "/cluster_", cohort_name, "_trajectory_basal3.txt")
readr::write_delim(pr_deg_ids, filename, delim="\t")
