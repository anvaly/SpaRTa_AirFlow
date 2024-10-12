# Authors: Yulong Bai, Hannah Pliner, Anna Lyubetskaya. Date: 23-03-10
# A variety of assist functions for 10X_to_Seurat_ST.R


spot_clean_convertToSeurat_my <- function(slide_obj, image_dir, slice="slice1", filter_matrix=TRUE){
  ## Since Read10X_Image function is used, for CytAssist samples, original convertToSeurat function from SpotClean converts coordinates features to character.
  ## Line 14 added to make sure the coordinates features are integer. 
  
  # create Seurat object
  object <- SeuratObject::CreateSeuratObject(SummarizedExperiment::assay(slide_obj), assay = "Spatial")
  
  # load image and add to Seurat object
  image <- Seurat::Read10X_Image(image.dir = image_dir, filter.matrix = filter_matrix)
  image@coordinates[colnames(image@coordinates)] <- lapply(image@coordinates[colnames(image@coordinates)], as.integer)
  ts_coord <- SeuratObject::GetTissueCoordinates(image)
  
  if(nrow(ts_coord)>ncol(object)){
    stop("The slide object has fewer spots than the image data. Consider setting filter_matrix=TRUE?")
  }
  if(nrow(ts_coord)<ncol(object)){
    object <- object[,rownames(ts_coord)]
  }
  
  image <- image[Seurat::Cells(x = object)]
  Seurat::DefaultAssay(object = image) <- "Spatial"
  object[[slice]] <- image
  
  return(object)
}


find_top_genes_my <- function(orig_seur, decont_seur) {
  ## Hannah's function to identify top changed genes from before and after spotclean Seurat objects
  ## This function can be used to pick up genes for spotclean figures when there is no good markers from prior knowledge
  
  decont_mat <- decont_seur@assays$Spatial@counts
  orig_mat <- orig_seur@assays$Spatial@counts
  both_genes <- intersect(row.names(decont_mat), row.names(orig_mat))
  both_spots <- intersect(colnames(decont_mat), colnames(orig_mat))
  orig_mat <- orig_mat[both_genes,both_spots]
  decont_mat <- decont_mat[both_genes,both_spots]
  diff_mat <- orig_mat - decont_mat
  
  gene_stats <- data.frame(gene = row.names(diff_mat),
                           total_change = Matrix::rowSums(diff_mat),
                           total_abs_change = Matrix::rowSums(abs(diff_mat)),
                           mean_change = Matrix::rowMeans(diff_mat),
                           mean_abs_change = Matrix::rowMeans(abs(diff_mat)),
                           mean_percent_change = Matrix::rowMeans((diff_mat * 100)/(orig_mat + .001)))
  
  gene_stats <- gene_stats[order(gene_stats$total_abs_change, decreasing = T),]
  
  return(gene_stats)
}


apply_spotclean_my <- function(filename_in, filename_out4){
  ## Apply SpotClean to the current Visium object
  
  # Start with raw data including non-tissue spots
  seur_raw <- SpotClean::read10xRawH5(paste0(filename_in, "/raw_feature_bc_matrix.h5"))
  
  # Read SpaceRanger calc-ed positions
  if(file.exists(paste0(filename_in, "/spatial/tissue_positions.csv"))) {
    tissue_filename <- paste0(filename_in, "/spatial/tissue_positions.csv")
  } else if(file.exists(paste0(filename_in, "/spatial/tissue_positions_list.csv"))) {
    tissue_filename <- paste0(filename_in, "/spatial/tissue_positions_list.csv")
  } else {
    stop("Can't find tissue positions file.")
  }
  
  # Start with raw data including non-tissue spots
  slide_info <- SpotClean::read10xSlide(tissue_filename,
                                        paste0(filename_in, "/spatial/tissue_lowres_image.png"),
                                        paste0(filename_in, "/spatial/scalefactors_json.json"))
  
  # Fix the issue with an extra coordinates in the image slot
  cells <- intersect(colnames(seur_raw), slide_info$slide$barcode)
  slide_info$slide <- slide_info$slide %>%
    dplyr::filter(barcode %in% cells)
  
  slide_obj <- SpotClean::createSlide(count_mat = seur_raw, slide_info = slide_info, gene_cutoff = 0)
  
  # Decontaminate raw data
  decont_obj <- SpotClean::spotclean(slide_obj, verbose = TRUE, kernel = "gaussian", gene_cutoff = 0)
  
  # Convert to a Seurat object
  ## convertToSeurat function of SpotClean pacakge uses Seurat Read10X_Image function,
  ## Use modified version to make sure coordinates are integer
  # decont_seur <- SpotClean::convertToSeurat(decont_obj, image_dir = paste0(filename_in, "/spatial/"))
  decont_seur <- spot_clean_convertToSeurat_my(decont_obj, image_dir = paste0(filename_in, "/spatial/"))
  decont_seur$contamination_rate <- decont_obj@metadata$contamination_rate
  
  # List of genes present in both object - before and after SpotClean
  ref_gene_list <- intersect(rownames(data_seurat_filt), rownames(decont_seur))
  
  # Generate a few figs for checking SpotClean
  p_list <- list()
  
  # Spots under and outside of tissue
  p_list[[1]] <- SpotClean::visualizeLabel(slide_obj, "tissue") + 
    ggplot2::ggtitle(paste(sample_name, "tissue")) + 
    ggplot2::scale_fill_manual(name = NULL, labels = c("Empty", "Tissue"), 
                               values = c("#5D2E8C", "#F96E46")) + theme_classic() + 
    ggplot2::theme(legend.position="top") + 
    Seurat::NoAxes() 
  
  # Contamination rates by spot
  p_list[[2]] <- spatial_feature_plot_my(decont_seur, feature="contamination_rate", title=paste(sample_name, "SpotClean contamination rate"))
  
  # Compare before and after SpotClean objects to find most changed genes
  top_gene_mat <- find_top_genes_my(orig_seur = data_seurat_filt, decont_seur = decont_seur)
  top_gene_mat_output <- paste0(output_folders[["2d_Matrices_spotclean_topgene"]], sample_name, "_spotclean_topgene.txt")
  write.table(top_gene_mat, top_gene_mat_output, quote = F, sep = "\t", row.names = F)
  
  # Select and plot genes before and after SpotClean
  if(is.numeric(spotclean_gene_list) & length(spotclean_gene_list)==1){
    topNgene <- spotclean_gene_list
    spotclean_gene_list <- top_gene_mat$gene[1:topNgene]
  }
  genes2plot <- intersect(ref_gene_list, spotclean_gene_list)
  
  if(length(genes2plot) > 0){
    for(g in intersect(ref_gene_list, spotclean_gene_list)){
      p_list[[length(p_list) + 1]] <- spatial_feature_plot_my(data_seurat_filt, g, title=paste(g, "before SpotClean"))
      p_list[[length(p_list) + 1]] <- spatial_feature_plot_my(decont_seur, g, title=paste(g, "after SpotClean"))
    }
  }
  
  # Write figures to file
  write_plot2file_my(patchwork::wrap_plots(p_list, nrow=length(intersect(ref_gene_list, spotclean_gene_list))+1, ncol=2), 
                     filename_out4, num_row=length(intersect(ref_gene_list, spotclean_gene_list))+1, num_col=2)
  
  # Only include genes/spots from filtered 10X in decontaminated
  data_seurat_init <- decont_seur[row.names(data_seurat_filt), colnames(data_seurat_filt)]
  
  # Add the original raw counts just in case we want to see data before SpotClean
  data_seurat_init@assays$Spatial_init <- data_seurat_filt@assays$Spatial
  
  
  return(data_seurat_init)
}


gen_spotdetect_df_my <- function(outs_dir, positions_csv, max_error = 50) {
  ## Generate the error stats for spot placement
  
  # Read SpaceRanger calc-ed positions
  if(file.exists(paste0(outs_dir, "/spatial/tissue_positions.csv"))) {
    tissue_filename <- paste0(outs_dir, "/spatial/tissue_positions.csv")
  } else if(file.exists(paste0(outs_dir, "/spatial/tissue_positions_list.csv"))) {
    tissue_filename <- paste0(outs_dir, "/spatial/tissue_positions_list.csv")
  } else {
    stop("Can't find tissue positions file.")
  }
  
  tissue_pos <- read.csv(tissue_filename, header = FALSE)
  if(tissue_pos[1,1] == "barcode") {
    tissue_pos <- read.csv(tissue_filename, header = TRUE)
  }
  names(tissue_pos) <- c("barcode", "tissue", "col", "row", "Y", "X")
  tissue_pos$type <- "provided"
  row.names(tissue_pos) <- tissue_pos$barcode
  
  # Read spotdetect.py calc-ed positions
  calc_pos <- read.csv(positions_csv, header = FALSE)
  names(calc_pos) <- c("X", "Y", "radius")
  
  # Convert pixel to um using scale factor
  scale_factors <- suppressWarnings(read.table(paste0(outs_dir, 
                                                      "/spatial/scalefactors_json.json"), 
                                               sep = ","))
  
  spot_diam <- scale_factors[,grepl("spot_diameter_fullres", scale_factors[1,])]
  spot_diam <- as.numeric(gsub("}", "", stringr::str_split(spot_diam, ":")[[1]][2]))
  mppx <- 55/spot_diam
  
  scale_factor <- scale_factors[,grepl("tissue_hires_scalef", scale_factors[1,])]
  scale_factor <- as.numeric(gsub("}", "", stringr::str_split(scale_factor, ":")[[1]][2]))
  
  calc_pos$X_scale <- calc_pos$X/scale_factor
  calc_pos$Y_scale <- calc_pos$Y/scale_factor
  calc_pos$type <- "calculated"
  
  # Merge table for calculations by matching closest spots
  joint_df <- calc_pos[,c("X_scale", "Y_scale", "type")]
  names(joint_df) <- c("X", "Y", "type")
  
  joint_df <- rbind(joint_df, tissue_pos[,c("X", "Y", "type")])
  dist_mat <- proxy::dist(joint_df[joint_df$type == "calculated",1:2],
                          y=joint_df[joint_df$type == "provided",1:2])
  crossmat <- `dim<-`(c(dist_mat), dim(dist_mat))
  row.names(crossmat) <- row.names(dist_mat)
  colnames(crossmat) <- colnames(dist_mat)
  
  dist_df <- reshape2::melt(as.matrix(crossmat))
  dist_df <- dist_df %>% dplyr::group_by(Var1) %>% dplyr::top_n(-1, value)
  
  all_df <- merge(dist_df, calc_pos, by.x = "Var1", by.y = 0)
  names(all_df) <- c("obs_spot", "barcode", "distance", "X", 
                     "Y", "obs_radius", "X_obs", "Y_obs", "type")
  all_df <- all_df[,c("obs_spot", "barcode", "distance", "obs_radius", 
                      "X_obs", "Y_obs")]
  all_df <- merge(all_df, tissue_pos[,c("barcode", "tissue", "X", "Y")], 
                  by = "barcode")
  
  # Calc error metrics
  all_df$distance_um <- all_df$distance * mppx
  
  all_df$X_dev <- all_df$X - all_df$X_obs
  all_df$Y_dev <- all_df$Y - all_df$Y_obs
  
  all_df$X_dev_um <- all_df$X_dev * mppx
  all_df$Y_dev_um <- all_df$Y_dev * mppx
  all_df$mppx <- mppx
  all_df$obs_diam <- 2 * all_df$obs_radius
  all_df$obs_diam_um <- all_df$obs_diam * mppx
  all_df$diam <- spot_diam
  all_df$diam_um <- spot_diam * mppx
  
  # Remove spots that are probably badly matched
  all_df <- all_df[all_df$tissue == 0,]
  all_df <- all_df[all_df$distance_um < max_error,]
  
  return(all_df)
}


generate_spotdetect_summary_fig_my <- function(error_df, sample_name, output_file) {
  ## Function to make spotdetect summary figure
  
  library(patchwork)
  library(ggplot2)
  
  pred <- ggplot(data.frame(X = c(error_df$X, error_df$X_obs),
                            Y = c(error_df$Y, error_df$Y_obs),
                            Type = c(rep("Expected", nrow(error_df)), 
                                     rep("Observed", nrow(error_df)))), 
                 aes(X, -Y, color = Type)) + 
    geom_point(alpha = 0.9, size = 0.3) + 
    theme_bw(base_size = 6) + 
    ggtitle(label = "", subtitle = "(diameters not representative)") + 
    labs(x = "", y = "") + scale_color_manual(values = c("black", "#F4C95D")) + 
    guides(colour = guide_legend(override.aes = list(size=3)))
  
  x_diff <- ggplot(error_df, aes(X_dev_um)) + geom_density() + 
    theme_bw(base_size = 6) + labs(x = "Difference (um)", y = "Density") + 
    geom_vline(xintercept = 0, color = 'red') + 
    ggtitle(label = "Deviation in X direction", 
            subtitle = paste0("Mean: ", 
                              round(mean(error_df$X_dev_um), 2), ", median: ", 
                              round(median(error_df$X_dev_um), 2), ", sd: ", 
                              round(sd(error_df$X_dev_um), 2)))
  
  y_diff <- ggplot(error_df, aes(Y_dev_um)) + geom_density() + 
    theme_bw(base_size = 6) + labs(x = "Difference (um)", y = "Density") + 
    geom_vline(xintercept = 0, color = 'red') + 
    ggtitle(label = "Deviation in Y direction", 
            subtitle = paste0("Mean: ", 
                              round(mean(error_df$Y_dev_um), 2), ", median: ", 
                              round(median(error_df$Y_dev_um), 2), ", sd: ", 
                              round(sd(error_df$Y_dev_um), 2)))
  
  tot_dist <- ggplot(error_df, aes(distance_um)) + geom_density() + 
    theme_bw(base_size = 6) + labs(x = "Distance (um)", y = "Density") + 
    ggtitle(label = "Distance from observed to expected", 
            subtitle = paste0("Mean: ", 
                              round(mean(error_df$distance_um), 2), ", median: ", 
                              round(median(error_df$distance_um), 2), ", max: ", 
                              round(max(error_df$distance_um), 2)))
  
  diam <- ggplot(error_df, aes(obs_diam_um)) + geom_density() + 
    theme_bw(base_size = 6) + labs(x = "Diameter (um)", y = "Density") + 
    ggtitle(label = "Observed diameter", 
            subtitle = paste0("Mean: ", 
                              round(mean(error_df$obs_diam_um), 2), ", median: ", 
                              round(median(error_df$obs_diam_um), 2), ", max: ", 
                              round(max(error_df$obs_diam_um), 2))) +
    geom_vline(xintercept = 55)
  
  patchwork <- ((x_diff | y_diff) / (tot_dist ) | pred)
  png(output_file, width = 9, height = 5, res = 300, units = "in")
  print(patchwork + patchwork::plot_annotation(
    title = sample_name,
    subtitle = 'Error between observed (image analysis) and expected (SpaceRanger calculated) spot positions'
  ))
  dev.off()
  
  return(NULL)
}
