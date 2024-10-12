
# Allow piping throughout the package
`%>%` <- magrittr::`%>%`


source("code/utils/utils_tibble.R")
source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_10X_in_out.R")
source("code/R/Utils/utils_10X_in_out_perturb.R")
source("code/R/Utils/utils_10X_vis.R")
source("code/R/Utils/utils_10X_dea.R")
source("code/R/Utils/utils_10X_object.R")
source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_perturb.R")


## PARAMETERS ----


guide_umi_top_threshold <- 50
guide_umi_top_ratio <- 5

pca_var_threshold <- 0.2
npcs <- 10

cohort_regex <- "MC38_Glides"

## PATHS ----


# Input location containing 10X folders
input_path <- "XXXX"

# Output location for figures
output_path <- "XXXX"
output_figs <- "XXXX"
# Create output folder
dir.create(output_figs, showWarnings=FALSE)

# Input file with sample meta data
meta_file <- paste0(input_path, "PerturbSeq_Metadata_Sheet.txt")


## INGEST, WRANGLE, OUTPUT ----

# Read sample file
meta_df <- readr::read_delim(file=meta_file, delim="\t") %>%
  tidyr::drop_na() %>%
  dplyr::filter(grepl(cohort_regex, Cohort))



# for(i in 1:nrow(meta_df)){
  i <-1

  
  # Meta data parameters
  filename_in <- meta_df[[i, "FullPath"]]
  filename_gex_in <- paste0(filename_in, "filtered_feature_bc_matrix/matrix_gex.mtx")
  filename_crispr_in <- paste0(filename_in, "filtered_feature_bc_matrix/matrix_crispr.mtx")
  sample_name <- meta_df[[i, "Sample_Name"]]
  
  filename_out1 <- paste0(output_path, sample_name, "_all.rds")
  filename_out2 <- paste0(output_path, sample_name, "_crispr.txt")
  filename_out3 <- paste0(output_path, sample_name, "_crispr_summary.txt")
  
  path_out <- paste0(output_figs, sample_name)
  dir.create(path_out, showWarnings=FALSE)
  ## Load data ----
  
  # Read in GEX data
  data_seurat <- matrix_to_Seurat_my(filename_gex_in)
  
  # Read the CRISPR data
  data_seurat_crispr <- matrix_to_Seurat_my(filename_crispr_in)
  
  # Get names of guides in this sample
  guide_names <- data_seurat_crispr@assays$RNA@meta.features %>% rownames()
  num_guides <- guide_names %>% length
  guide_metadata_df <- tibble::as_tibble(data_seurat_crispr@meta.data)
  
  ## Plot distributions of guide library
  plot_title <- paste0("Guide UMI Count Distribution - ",sample_name)
  p1 <- ggplot(guide_metadata_df, aes(x = nCount_RNA))+
    geom_histogram(bins=50) +
    scale_x_log10()+
    scale_y_continuous(expand = c(0,0))+
    ggtitle(plot_title)+
    theme_classic()

  plot_title <- paste0("Feature Count Distribution - ",sample_name)
  p2 <- ggplot(guide_metadata_df, aes(x = nFeature_RNA))+
    geom_histogram(binwidth = 1) +
    scale_x_continuous(breaks = seq(1,num_guides,3))+
    scale_y_continuous(expand = c(0,0))+
    ggtitle(plot_title)+
    theme_classic()
  filename <- paste0(path_out, "/guide_library_metrics")
  write_plot2file_my(in_plot = patchwork::wrap_plots(p1+p2), filename = filename,  num_col = 2, num_row =1)
  
  ## Guide Library Entropy ----
  # A higher entropy indicates a less confident assignment
  # Normalize the data by relative counts
  data_seurat_crispr <- Seurat::NormalizeData(data_seurat_crispr, normalization.method = "RC", scale.factor = 1) %>% Seurat::ScaleData()
  
  # Calculate Entropy of guide library over UMIs per cell.
  normalized_guide_counts <- data_seurat_crispr@assays$RNA@data # each cell sums to one, so this is a probability
  guide_entropy <- calculate_guide_entropy_my(normalized_guide_counts)
  # Add to Seurat Object Metadata
  data_seurat_crispr@meta.data$guide_entropy <- guide_entropy
  
  # update standalone metadata df for plotting
  guide_metadata_df <- dplyr::mutate(guide_metadata_df, guide_library_entropy =guide_entropy)
  guide_metadata_df <- dplyr::mutate(guide_metadata_df, cell_bc = rownames(data_seurat_crispr@meta.data))
  
  # Plot entropy histogram
  plot_title <- paste0("Guide Library Entropy -", sample_name)
  p1<- ggplot(guide_metadata_df, aes(x =guide_library_entropy))+
    geom_histogram()+
    scale_y_continuous(expand = c(0,0))+
    ggtitle(plot_title)+
    theme_classic()
  
  #plot UMIs vs entropy
  plot_title <- paste0("UMI Count - ",sample_name)
  p2<- ggplot(guide_metadata_df, aes(x = nCount_RNA, y = guide_library_entropy))+
    geom_point(alpha = 0.1)+
    scale_x_log10()+
    ggtitle(plot_title)+
    theme_classic()

  #plot number of detected features vs entropy
  plot_title <- paste0("Features - ",sample_name)
  p3 <- ggplot(guide_metadata_df, aes(x = nFeature_RNA, y = guide_library_entropy))+
    geom_point(alpha = 0.1, position = "jitter")+
    ggtitle(plot_title)+
    theme_classic()
  
  filename <- paste0(path_out,"/entropy_plots")
  write_plot2file_my(in_plot = patchwork::wrap_plots(p1+p2+p3), filename = filename, num_col = 3, num_row =1)
  
  # GEX vs Guide Library Depth
  # Curious if GEX library size is dependent on GEX library size (i.e. does more guide umis mean more gex umis?)
  gex_metadata_df <- tibble::as_tibble(data_seurat@meta.data)
  gex_metadata_df <- dplyr::mutate(gex_metadata_df, cell_bc = rownames(data_seurat@meta.data))
  
  #merge gex and guide metadata dfs
  merged_metadata_df <- dplyr::left_join(gex_metadata_df, guide_metadata_df, by = "cell_bc", suffix = c(".gex", ".guide"))
  
  plot_title <- paste0("Guide v GEX UMI")
  p <- ggplot(merged_metadata_df, aes(x = nCount_RNA.gex, y = nCount_RNA.guide))+
    geom_point(alpha = 0.3)+
    # scale_y_log10()+
    # scale_x_log10()+
    ggtitle(plot_title)+
    theme_classic()
  filename <- paste0(path_out, "/gex_vs_guide_library_umi")
  write_plot2file_my(in_plot = p, filename = filename, width = 4, height = 3)
  
  ## PCA and UMAP on the guide library
  data_seurat_crispr <- seurat_pca_my(data_seurat_crispr, npcs = npcs, var_threshold=pca_var_threshold, features = guide_names)
  data_seurat_crispr <- seurat_umap_nb_my(data_seurat_crispr$data, num_dimensions=data_seurat_crispr$num_pcs)
  
  # 
  p1 <- feature_plot_my(data_seurat_crispr, var = c("nCount_RNA")) + ggtitle("UMIs") + theme(aspect.ratio = 1, axis.title = element_blank())
  p2 <- feature_plot_my(data_seurat_crispr, var = c("nFeature_RNA")) + ggtitle("Features") + theme(aspect.ratio = 1, axis.title = element_blank())
  p3 <- feature_plot_my(data_seurat_crispr, var = c("guide_entropy")) + ggtitle("Entropy") + theme(aspect.ratio = 1, axis.title = element_blank())
  filename <- paste0(path_out, "/guide_UMAP")
  write_plot2file_my(in_plot = patchwork::wrap_plots(p1,p2,p3, ncol = 3), filename = filename, num_col = 3, num_row = 1, width = 10, height = 4)
  # }
  

## Process the guide data ----

# Guide data as a wide tibble
crispr_wide_df <- tibble::as_tibble(as.matrix(Seurat::GetAssayData(data_seurat_crispr,slot = "counts")), rownames="Guides")

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

# plot top guide vs 2nd highest guide per cell
ggplot(crispr_stat_wide_df, aes(x = Guide1, y = Guide2))+
  geom_point(alpha = 0.2)+
  scale_y_continuous(expand = c(.01,0))+
  scale_x_continuous(expand = c(.01,0))+
  theme_classic()

crispr_stat_wide_df <- dplyr::rename(crispr_stat_wide_df, cell_bc = Coordinate)
full_metadata <- dplyr::left_join(crispr_stat_wide_df, merged_metadata_df, by = "cell_bc")

ggplot(full_metadata, aes(x = nFeature_RNA.guide, y = Guide1))+
  geom_jitter(alpha = 0.1)+
  theme_classic()

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
data_seurat_crispr@meta.data[["CRISPRGeneCategory"]] <- crispr_stat_wide_df$Category
data_seurat_crispr@meta.data[["CRISPRGeneCategory"]] <- gsub("\\.\\d+$", "", crispr_stat_wide_df$Category)

p5 <- dim_plot_my(data_seurat_crispr, group.by = "CRISPRGeneCategory") + theme(aspect.ratio = 1)
write_plot2file_my(in_plot = p5, filename = paste0(output_figs, "guide_umap_grouped_by_guide_all_Cells.png"), width = 8, height = 6)
  
coordinates_select <- crispr_stat_wide_df$cell_bc[which(crispr_stat_wide_df$Category != "")]
data_seurat_crispr_subset <- subset(data_seurat_crispr, cells=coordinates_select)
p4 <- dim_plot_my(data_seurat_crispr_subset, group.by = "CRISPRGeneCategory") + theme(aspect.ratio = 1)
write_plot2file_my(in_plot = p4, filename = paste0(output_figs, "guide_umap_grouped_by_guide.png"), width = 8, height = 6)


raw_guide_counts <- data_seurat_crispr@assays$RNA@counts
normcounts_df <- tibble::as_tibble(Matrix::t(normalized_guide_counts))
normcounts_df <- dplyr::mutate(normcounts_df, cell_bc = colnames(normalized_guide_counts), .before =1)

rawcounts_df <- tibble::as_tibble(Matrix::t(raw_guide_counts))
rawcounts_df <- dplyr::mutate(rawcounts_df, cell_bc = colnames(raw_guide_counts), .before =1)
# 

metadata_df <- data_seurat_crispr@meta.data %>% tibble::as_tibble()
metadata_df$cell_bc <- rownames(data_seurat_crispr@meta.data) 

norm_data_df <- dplyr::left_join(normcounts_df,metadata_df, .before = 2) 
raw_data_df <- dplyr::left_join(rawcounts_df, metadata_df, .before = 2)
# 
norm_data_df_long <- tidyr::pivot_longer(norm_data_df, cols = guide_names, names_to = "Guide", values_to = "relative_guide_counts")
raw_data_df_long <- tidyr::pivot_longer(raw_data_df, cols = guide_names, names_to = "Guide", values_to = "raw_guide_counts")
  

  # 
  ggplot(raw_data_df_long, aes(x = raw_guide_counts))+
    geom_histogram()+
    scale_x_log10("Guide Counts")+
    scale_y_continuous("Count (cells)")+
    facet_wrap(.~Guide)+
    coord_cartesian(ylim = c(0, 300))+
    theme_classic()
  # 
  
  #plot heatmap of normalized and raw guide library
  filename <- paste0()
  pheatmap::pheatmap(Matrix::t(data_seurat_crispr@assays$RNA@data))
  raw_heatmap <- pheatmap::pheatmap(Matrix::t(data_seurat_crispr@assays$RNA@counts))
  
  
  #}

