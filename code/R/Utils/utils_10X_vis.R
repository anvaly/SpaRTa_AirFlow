# Author: Anna Lyubetskaya. Date: 20-04-16
# Helper functions for ST analysis


library(ggplot2)
source("code/utils/utils_ggplot.R")
source("code/R/Utils/utils_colormap.R")


plot_format_my <- function(p, title="", x="", y=""){
  ## For a random Seurat plot, update fonts, labels, and legend
  
  base_font <- 5
  
  p <- p + 
    ggplot2::labs(x=x, y=y, title=title) + 
    ggplot2::theme(plot.title = element_text(size=base_font+2),
                   legend.title = element_blank(), 
                   legend.text = element_text(size=base_font, hjust = 1, angle = 0),
                   legend.key.size = unit(0.25, "line"), 
                   legend.position = "right")
  
  return(p)
  
}


find_pt_size_factor_my <- function(data_seurat){
  ## Define the size of the circle for tissue visualizations
  
  if("user.pt.size.factor" %in% names(data_seurat@misc)){
    if(!is.na(data_seurat@misc[["user.pt.size.factor"]])){
      pt.size.factor <- data_seurat@misc[["user.pt.size.factor"]]
    }else{
      pt.size.factor <- 1.5
    }
  } else{
    pt.size.factor <- 1.5
  }
  
  return(pt.size.factor)
}


find_col_num_my <- function(data_seurat, var){
  ## Find the number of unique values in a specific variable
  
  return(length(unique(data_seurat@meta.data[, var])))
}


define_cols_for_var_my <- function(data_seurat, var, col_names=NULL, col_type="jet"){
  ## Define a color scheme for a specific variable and return a named list
  
  # Find variable values
  if(is.factor(data_seurat@meta.data[[var]])){
    levels_loc <- levels(data_seurat@meta.data[[var]])
  } else{
    levels_loc <- sort(unique(data_seurat@meta.data[[var]]))
  }
  
  # Number of unique values in a variable
  n <- length(levels_loc)
  
  # Establish color scheme
  cols <- define_cols_my(n=n, col_type=col_type)
  
  if(is.null(col_names)){
    names(cols) <- levels_loc
  } else{
    names(cols) <- col_names
  }
  
  return(cols)
}


set_spatial_min_max_my <- function(data_seurat, feature_list, min_val=0, max_val=1){
  ## Set a min and a max value to a Seurat object for a spatial plot
  
  # Identify coordinates for two corner spots
  barcode_coord_list <- data_seurat@images[[1]]@coordinates %>% 
    dplyr::arrange(row, col)
  
  # Select only those meta data features that are present in the Seurat object
  sig_list <- intersect(colnames(data_seurat@meta.data), feature_list)
  
  # Set minimum value
  data_seurat@meta.data[rownames(barcode_coord_list)[1], sig_list] <- min_val
  # Set maximum value
  data_seurat@meta.data[rownames(barcode_coord_list)[2], sig_list] <- max_val
  
  return(data_seurat)
}


violin_visualization_gene_list_my <- function(data_seurat, feature_list, output_file_name){
  ## Create a composite violin plot for a short feature list and a spatial Seurat object
  
  # Initialize plot list
  p_list <- list()
  
  # Create individual plots for each marker to visualize
  for(feature in feature_list){
    p_list[[feature]] <- Seurat::VlnPlot(data_seurat, 
                                         features = feature, 
                                         pt.size = 0.1) %>%
      plot_format(feature) +
      theme(legend.position = "none")
  }
  
  # Write the composite image to a file
  write_plot2file_my(patchwork::wrap_plots(p_list, nrow=1, ncol=length(feature_list)), output_file_name)
}


seurat_heatmap_my <- function(data_seurat, feature_list, plot_title=plot_title, filename=NULL){
  ## Create a heatmap using Seurat
  
  p <- Seurat::DoHeatmap(data_seurat, 
                         features = feature_list, 
                         #group.by = categ, 
                         label = T, 
                         assay = "SCT", 
                         #cells = cell_list,
                         draw.lines = T, 
                         group.bar.height = 0, 
                         raster = FALSE) + 
    theme(legend.position = "bottom") + 
    ggtitle(plot_title)
  
  # Write the composite image to a file
  if(!is.null(filename)){
    write_plot2file_my(p, filename)
  }
  
  return(p) 
}


dim_plot_my <- function(data_seurat, group.by, split.by = NULL, pt.size=0.5, reduction=NULL, cols=NULL, ncolumns = 1){
  ## Create a dimension plot
  
  # Establish color scheme for the reference
  if(is.null(cols)){
    cols <- define_cols_for_var_my(data_seurat, group.by)
  }
  
  p <- Seurat::DimPlot(object = data_seurat, 
                       group.by = group.by, 
                       pt.size = pt.size,
                       split.by = split.by,
                       reduction = reduction,
                       cols = cols, 
                       ncol = ncolumns,
                       raster = FALSE) %>%
    plot_format_my(title=paste(reduction, "of", group.by), x="Dim 1", y="Dim 2")
  
  return(p)
}


feature_plot_my <- function(data_seurat, var, pt.size=0.25, reduction=NULL, min.cutoff="q1", max.cutoff="q99"){
  ## Create a feature plot
  
  p <- Seurat::FeaturePlot(object = data_seurat, 
                           features = var, 
                           pt.size = pt.size,
                           reduction = reduction,
                           min.cutoff = min.cutoff,
                           max.cutoff = max.cutoff) %>%
    plot_format_my(title=paste(reduction, "of", var), x="Dim 1", y="Dim 2")
  
  return(p)
}


spatial_dim_plot_my <- function(data_seurat, cell_highlight=NULL, group.by=NULL, title="", cols=NULL, 
                                images=NULL, combine=TRUE, num_col=NULL){
  # Create a spatial plot highlighting specific spots
  
  ## THIS IS A TEMPORARY FIX TO THIS BUG: https://github.com/satijalab/seurat/issues/6179
  if(!is.null(group.by)){
    Seurat::Idents(data_seurat) <- data_seurat[[group.by]]
  }
  
  # Font size
  base_font <- 6
  
  if(!is.null(cell_highlight)){
    facet.highlight <- TRUE
    cols.highlight <- c("darkblue", "grey100")
  } else{
    facet.highlight <- FALSE
    cols.highlight <- NULL
    
    # Establish color scheme for the reference
    if(is.null(cols)){
      cols <- define_cols_for_var_my(data_seurat, group.by)
    }
  }
  
  # ncol doesn't seem to work
  p <- Seurat::SpatialDimPlot(data_seurat, 
                              # group.by = group.by, # THIS IS A TEMPORARY FIX TO THIS BUG: https://github.com/satijalab/seurat/issues/6179
                              cols = cols, 
                              pt.size.factor = find_pt_size_factor_my(data_seurat),
                              cells.highlight = cell_highlight, 
                              facet.highlight = facet.highlight, 
                              cols.highlight = cols.highlight, 
                              images = images, 
                              combine = combine,
                              stroke = 0,
                              ncol = num_col
                              )
  
  # Seurat doesn't wrap the plot if combine = FALSE
  if(combine == FALSE){
    p <- patchwork::wrap_plots(p)  
  }
  
  if(length(names(data_seurat@images)) == 1){  #  || length(images) == 1
    p <- p %>%
      plot_format_my(title=title, x="", y="")
  } else if(combine == TRUE){
    p <- p +
      patchwork::plot_layout(guides="collect") +
      patchwork::plot_annotation(title=title)
  }
  
  return(p)
}


spatial_feature_plot_my <- function(data_seurat, feature, min.cutoff="q1", max.cutoff="q99", title="", 
                                    crop=TRUE, slot="data", color=NULL, images=NULL){
  ## Create a single spatial feature plot
  
  p <- Seurat::SpatialFeaturePlot(object = data_seurat, 
                                  features = feature,
                                  alpha = 1, 
                                  ncol = 1, 
                                  pt.size.factor = find_pt_size_factor_my(data_seurat), 
                                  min.cutoff = min.cutoff, 
                                  max.cutoff = max.cutoff, 
                                  crop = crop, 
                                  slot = slot,
                                  images = images, 
                                  stroke = 0)
  
  if(!is.null(color)){
    p <- p +
      ggplot2::scale_fill_gradientn(colors = color)
  }
  
  if(length(names(data_seurat@images)) == 1 || length(images) == 1){
    p <- p %>%
      plot_format_my(title)
  } else{
    base_font <- 8
    
    p <- p +
      patchwork::plot_annotation(title=title)
  }
  
  return(p)
}


batch_spatial_feature_plot_my <- function(seurat_list, feature_list, output_file=NULL, title="",
                                          min.cutoff="q1", max.cutoff="q99", slot="data", plot_HE=FALSE, 
                                          width=NULL, height=NULL, plot_type=rep("f", length(feature_list)),
                                          image_list=NULL, guides=NULL, col_list=list()){
  
  ## Create a composite spatial plot for a short feature list and a spatial Seurat object
  
  # Initialize plot list
  p_list <- list()
  
  # Cycle through samples
  for(dataset in names(seurat_list)){
    
    # Select a dataset
    data_seurat <- seurat_list[[dataset]]

    if(is.null(image_list)){
      this_image <- names(data_seurat@images)
    }else{
      this_image <- image_list
    }
    
    for(im in intersect(this_image, names(data_seurat@images))){
      
      # Plot tissue underneath with no spot overlay
      if(plot_HE == TRUE){
        name <- paste(dataset, "\n", "H&E Image")
        p_list[[name]] <- Seurat::SpatialDimPlot(data_seurat, pt.size.factor = 0, images=im) + 
          ggtitle(name) + 
          theme(legend.position = "none", plot.title = element_text(size=8))
      }
      
      # Cycle through features to plot
      for(i in 1:length(feature_list)){
        feature <- feature_list[[i]]
        plot_t <- plot_type[[i]]
        name <- paste(dataset, "\n", im, "\n", feature)
        
        # Create a single plot of a feature and a sample
        if(feature %in% colnames(data_seurat@meta.data) || gsub("^sig.", "", feature) %in% rownames(data_seurat)){
          
          if(plot_t == "f"){
            p_list[[name]] <- spatial_feature_plot_my(data_seurat, feature, title=name, images=im,
                                                      min.cutoff=min.cutoff, max.cutoff=max.cutoff, slot=slot)
          } else if(plot_t == "d"){
            cols <- NULL
            if(feature %in% names(col_list)){
              cols <- col_list[[feature]]
            }
            
            # This function returns the wrong thing when combine=FALSE
            p_list[[name]] <- spatial_dim_plot_my(data_seurat, group.by=feature, title=name, images=im, combine=FALSE, cols=cols)
          }
          
        }
      }
      
    }
  }
  
  # I don't know why this works but it does...
  # Seurat has a more developed spatial feature plotting function thatn spatial dim plotting
  if(length(unique(plot_type)) == 1 && unique(plot_type) == "d"){
    guides <- "collect"
  }
  
  # Write the composite image to a file
  if(!is.null(output_file) & length(p_list) > 0){
    
    # Define image size
    num_row <- ceiling(length(p_list) / 6)
    if(length(p_list) < 6){
      num_col <- length(p_list)
    } else{
      num_col <- 6
    }
    
    # Print image to file
    write_plot2file_my(patchwork::wrap_plots(p_list, nrow=num_row, ncol=num_col, guides=guides) +
                         patchwork::plot_annotation(title=title),  # theme(plot.title=element_text(size=8))
                       output_file, num_row=num_row, num_col=num_col, width=width, height=height)
  }
  
  return(p_list)
}


Seurat_pca_umap_spatial_my <- function(data_seurat, var, output_name, col_names=NULL, cols=NULL, title=NULL, pt.size=0.01){
  ## For a given variable, plot PCA, UMAP, and spatial distributions
  
  # Establish color scheme for the reference
  if(is.null(cols)){
    cols <- define_cols_for_var_my(data_seurat, var, col_names)
  }
  
  # Plot the variable data on the original tissue slice as well PCA/UMAP representations
  if(length(cols) > 1){
    
    # Add factors or sort existing factors
    data_seurat@meta.data[var] <- factor(data_seurat@meta.data[[var]], levels=names(cols))
    
    # Plot a classic two-dimensional projection: UMAP
    p1 <- dim_plot_my(data_seurat, var, reduction="pca", cols=cols, pt.size=pt.size)
    
    # Plot a classic two-dimensional projection: PCA
    p2 <- dim_plot_my(data_seurat, var, reduction="umap", cols=cols, pt.size=pt.size)
    
    # Define image size
    im_num <- length(names(data_seurat@images))
    num_row <- ceiling(im_num / 10)
    if(im_num < 10){
      num_col <- im_num
    } else{
      num_col <- 10
    }
    
    # Dummy up color list
    col_list <- list(var = cols)
    names(col_list) <- var
    
    # Plot clusters on the original tissue slice
    # p3 <- spatial_dim_plot_my(data_seurat, group.by=var, cols=cols, title=title, num_col=num_col, combine=FALSE)
    p3 <- batch_spatial_feature_plot_my(list("Cohort" = data_seurat), c(var), title=title, plot_type=c("d"),
                                        col_list=col_list)
    
    if(length(names(data_seurat@images)) == 1){
      
      # Write combo plot to file
      filename <- paste0(output_name, "_clust")
      # write_plot2file_my(patchwork::wrap_plots(list(p1, p2, p3), nrow=1), filename, num_row=1, num_col=6)
      write_plot2file_my(patchwork::wrap_plots(list(p2, patchwork::wrap_plots(p3)), nrow=1), filename, num_row=1, num_col=3)
      
    } else if (length(names(data_seurat@images)) > 1){
      
      # Write PCA and UMAP plot to file
      filename <- paste0(output_name, "_clust")
      # write_plot2file_my(patchwork::wrap_plots(list(p1, p2), nrow=1), filename, num_row=2, num_col=12)
      write_plot2file_my(p2, filename, num_row=2, num_col=6)
      
      # Write tissue plot to file
      filename <- paste0(output_name, "_clust_tissue")
      write_plot2file_my(patchwork::wrap_plots(p3, nrow=num_row, ncol=num_col, guides="collect"),
                         filename, num_row=num_row, num_col=num_col)
      
    } else{
      
      print("No image detected!")
      
    }
    
  } else{
    cols <- NULL
  }
  
  return(cols)
}


harmony_plot_spatial_clusters_my <- function(data_seurat, var, output_name, col_names=NULL, cols=NULL, title=NULL, image = NULL){
  ## Plot Harmony results
  
  # Establish color scheme for the reference
  if(is.null(cols)){
    cols <- define_cols_for_var_my(data_seurat, var, col_names)
  }
  
  # Plot the variable data on the original tissue slice as well PCA/UMAP representations
  if(length(cols) > 1){
    
    # Add factors or sort existing factors
    data_seurat@meta.data[var] <- factor(data_seurat@meta.data[[var]], levels=names(cols))
    
    # Plot a classic two-dimensional projection: UMAP
    p1 <- dim_plot_my(data_seurat, var, reduction="harmony", cols=cols, pt.size=0.05)+ggplot2::theme(aspect.ratio = 1)
    
    # Plot a classic two-dimensional projection: PCA
    p2 <- dim_plot_my(data_seurat, var, reduction="umap", cols=cols, pt.size=0.05)+ggplot2::theme(aspect.ratio = 1)
    
    # Plot clusters on the original tissue slice
    p3 <- spatial_dim_plot_my(data_seurat, group.by=var, cols=cols, title=title, images = image)
    
    p4 <- Seurat::SpatialDimPlot(data_seurat, pt.size = 0, images = image) + 
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(x = "", y = "")
    
    # Write combo plot to file
    filename <- paste0(output_name, "_clust")
    write_plot2file_my(patchwork::wrap_plots(list(p1, p2, p3, p4), nrow=1), filename, num_row=1, num_col=4)
    
    
  } else{
    cols <- NULL
  }
  
  return(cols)
}


harmony_cohort_plot <- function(data_seurat, var, split = NULL, output_name, col_names=NULL, cols=NULL, num_columns = 1, image = NULL){
  ## Plot Harmony results
  
  # Establish color scheme for the reference
  if(is.null(cols)){
    cols <- define_cols_for_var_my(data_seurat, var, col_names)
  }
  
  # Add factors or sort existing factors
  data_seurat@meta.data[var] <- factor(data_seurat@meta.data[[var]], levels=names(cols))
  num_samples <- length(unique(data_seurat@meta.data[[split]]))
  num_rows <- ceiling(num_samples/num_columns)
  
  # Plot a classic two-dimensional projection: Harmony-PCa
  p1 <- dim_plot_my(data_seurat, var, reduction="harmony", cols=cols, pt.size=0.05, ncolumns = num_columns, split.by = split)+ggplot2::theme(aspect.ratio = 1)
  filename <-paste0(output_name, "_pca")
  write_plot2file_my(p1, filename = filename, num_col = num_columns, num_row = num_rows)
  
  # Plot a classic two-dimensional projection: UMAP
  p2 <- dim_plot_my(data_seurat, var, reduction="umap", cols=cols, pt.size=0.05, ncolumns = num_columns, split.by = split)+ggplot2::theme(aspect.ratio = 1)
  filename <-paste0(output_name, "_umap")
  write_plot2file_my(p2, filename = filename, num_col = num_columns, num_row = num_rows)
  
  # Plot clusters on the original tissue slice
  p3 <- spatial_dim_plot_my(data_seurat, group.by=var, cols=cols, title=title, combine = FALSE)
  p3 <- patchwork::wrap_plots(p3, ncol = num_columns, nrow = num_rows, guides = "collect")
  filename <-paste0(output_name, "_spatial")
  write_plot2file_my(p3, filename = filename, num_col = num_columns, num_row = num_rows)
}
