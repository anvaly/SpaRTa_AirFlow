# Author: Anna Lyubetskaya. Date: 19-09-20
# These are ggplot function wrappers
# Plot types: violin, bar, hist, and scatter

# http://www.sthda.com/english/wiki/be-awesome-in-ggplot2-a-practical-guide-to-be-highly-effective-r-software-and-data-visualization

# todo: Set a different default pallette: https://stackoverflow.com/questions/10504724/change-the-default-colour-palette-in-ggplot
# todo: Make scatter plot more flexible in regards to shape and size
# todo: Reorder the stacked bar plot by total of Y axis
# todo: Pass "x" and "y" objects as variables to function plot_adjust_font
# todo: Remove the black border from the bar plot
# todo: Make violin plot summary stat visualization look better
# todo: Add slope and pvalue to the scatter plot
# todo: Coordinate font sizes


library(ggplot2)
options(ggplot2.qualitative.colour="Set1")


define_cols_my <- function(n=NULL, col_type="jet", cols=NULL){
  ## Define color palette for the plot
  ## Color types = c("random", "spectral", "rich", "rainbow", "viridis", "rgb", "jet")
  
  if(is.null(cols)){
    if(col_type == "random"){
      cols <- randomcoloR::distinctColorPalette(n, runTsne = FALSE)
    } else if(col_type == "spectral"){
      cols <- rev(RColorBrewer::brewer.pal(n, "Spectral"))
    } else if(col_type == "rich"){
      cols <- gplots::rich.colors(n, palette="temperature")
    } else if(col_type == "rainbow"){
      cols <- rainbow(n)
    } else if(col_type == "viridis"){
      cols <- viridis::plasma(n)
    } else if(col_type == "rgb"){
      cols <- colorRamps::rgb.tables(n)
    } else if(col_type == "jet"){
      col_func <- colorRampPalette(c("navy", "blue", "cyan", "darkgreen", "green", 
                                     "yellow", "orange", "brown", "red", "magenta", "purple", "gray"))
      cols <- col_func(n)
    }
  }
  
  # To visualize:
  # barplot(1:n, col=cols)
  
  return(cols)
}


create_violin_plot_my <- function(df, x_label, y_label, fill_label, facet_var=NULL, filename="violin", labels=NULL){
  ## Violin plot
  ## http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
  
  if(nrow(df) > 1){
    # Create axes
    p <- ggplot(df, aes(x=!!sym(x_label), y=!!sym(y_label), fill=!!sym(fill_label)))
    
    # Add a violin plot
    p <- p + geom_violin()
    
    # Add mean and standard deviation
    p <- p + stat_summary(fun.data="mean_sdl", geom="pointrange", color="black", size=0.25,
                          position = position_dodge(width = 0.9))

    # Tighten axes
    p <- p + 
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))
    
    # Format the plot by adding labels and facets and adjusting the format; wite the result to file
    p <- format_and_write_my(p, df, x_label, fill_label, facet_var, filename, labels)
    
    return(p)
  }
}


create_box_plot_my <- function(df, x_label, y_label, fill_label=NULL, facet_var=NULL, filename="bar", labels=NULL, 
                               with_dots=FALSE, with_line=NULL, reorder_x=FALSE, outlier_shape=19, cols=NULL){
  ## Box plot
  ## http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization
  ## Standard options for outlier.shape = NA or 19
  
  if(nrow(df) > 1){
    # Create axes
    
    # Order X labels by Y values in descending order
    if(reorder_x == TRUE){
      p <- ggplot(df, aes(x = reorder(!!sym(x_label), !!sym(y_label), FUN=mean), 
                          y = !!sym(y_label), fill=!!sym(fill_label)))
    } else{
      p <- ggplot(df, aes(x = !!sym(x_label), y = !!sym(y_label), fill=!!sym(fill_label)))
    }
    
    # Add a box plot
    p <- p + geom_boxplot(width=0.9, outlier.shape=outlier_shape)
    
    # Add jittered dot plot
    if(with_dots == TRUE){
      p <- p + geom_point(size=0.25, shape=16, position=position_jitterdodge(jitter.width=0.2, dodge.width=1), color="black")
    }
    # Add line between dots
    if(!is.null(with_line)){
      p <- p + geom_line(aes(group=!!sym(with_line)), alpha=0.3, linetype = "longdash")
    }
    # Change category colors
    if(!is.null(cols)){
      p <- p +
        scale_fill_manual(values=cols)
    }
    
    # Tighten axes
    p <- p + 
      scale_x_discrete(expand = c(0.05, 0.05)) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))

    # Format the plot by adding labels and facets and adjusting the format; wite the result to file
    p <- format_and_write_my(p, df, x_label, fill_label, facet_var, filename, labels)
    
    return(p)
  }
}


create_bar_plot_my <- function(df, x_label, y_label, fill_label=NULL, facet_var=NULL, position="stack", filename="bar", 
                               labels=NULL, error_label=NULL, reorder_x=FALSE, cols=NULL){
  ## Bar plot
  ## https://ggplot2.tidyverse.org/reference/geom_bar.html
  
  if(nrow(df) > 1){
    # Create axes
    # Pre-order X labels by Y values upstream of this function
    if(reorder_x == TRUE){
      p <- ggplot(df, aes(x=reorder(!!sym(x_label), !!sym(y_label)), weight=!!sym(y_label), fill=!!sym(fill_label)))
    } else{
      p <- ggplot(df, aes(x=!!sym(x_label), weight=!!sym(y_label), fill=!!sym(fill_label)))
    }

    # Add a bar plot
    p <- p + geom_bar(position=position, width=0.9)
    
    # Add error bars if data provided
    if(!is.null(error_label)){
      p <- p + geom_errorbar(aes(ymin=!!sym(y_label)-!!sym(error_label), 
                                 ymax=!!sym(y_label)+!!sym(error_label)), 
                             width=.1, position=position_dodge(0.9), color="black") 
    }
    
    # Change category colors
    if(!is.null(cols)){
      p <- p +
        scale_fill_manual(values=cols)
    }
    
    # Tighten axes
    p <- p + 
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.01)))

    # Format the plot by adding labels and facets and adjusting the format; wite the result to file
    p <- format_and_write_my(p, df, x_label, fill_label, facet_var, filename, labels)
    
    return(p)
  }
}


create_hist_plot_my <- function(df, x_label, fill_label=NULL, facet_var=NULL, intercept=0, binwidth=1, filename="hist", labels=NULL, 
                                add_density=FALSE, log_scale=FALSE){
  ## Histogram
  ## http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization
  
  if(nrow(df) > 1){
    # Create axes
    if(!is.null(fill_label)){
      p <- ggplot(df, aes(x = !!sym(x_label), fill = !!sym(fill_label)))
    } else{
      p <- ggplot(df, aes(x = !!sym(x_label)))
    }
    
    # Add a bar plot
    # About density: https://stackoverflow.com/questions/6967664/ggplot2-histogram-with-normal-curve
    p <- p + geom_histogram(binwidth=binwidth) + 
      geom_vline(xintercept=intercept, color="black", linetype="dashed", size=0.2)

    if(add_density == TRUE){
      p <- p + geom_density(aes(y=binwidth*..count..), alpha=0.25, fill="lightblue1", color="lightblue1")
    }
    
    # Tighten axes
    p <- p + 
      scale_x_continuous(expand = expansion(mult = c(0, 0))) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
    
    if(log_scale == TRUE){
      p <- p + scale_y_log10()
    }
    
    # Format the plot by adding labels and facets and adjusting the format; write the result to file
    p <- format_and_write_my(p, df, x_label, fill_label, facet_var, filename, labels)
    
    return(p)
  }
}


create_scatter_plot_my <- function(df, x_label, y_label, fill_label, shape=21, size=2, facet_var=NULL, dot_labels=NULL, 
                                   filename="dot", labels=NULL, do_fit=NULL, stroke=0.25, cols=NULL){
  ## Scatter plot
  ## http://www.sthda.com/english/wiki/ggplot2-scatter-plots-quick-start-guide-r-software-and-data-visualization
  
  if(nrow(df) > 1){
    if(is.null(dot_labels)){
      dot_labels <- fill_label
    }
    
    if(is.character(size)){
      # Create axes
      p <- ggplot(df, aes(x = !!sym(x_label), y = !!sym(y_label), label = !!sym(dot_labels), size=!!sym(size)))
      # Add a scatter plot: size = 1; shape = 21 by default
      p <- p + geom_point(aes(fill=!!sym(fill_label)), shape=shape, stroke=stroke)  # shape=!!sym(shape)), size=!!sym(size)
    } else{
      # Create axes
      p <- ggplot(df, aes(x = !!sym(x_label), y = !!sym(y_label), label = !!sym(dot_labels)))
      # Add a scatter plot: size = 1; shape = 21 by default
      p <- p + geom_point(aes(fill=!!sym(fill_label)), shape=shape, size=size, stroke=stroke)
    }
    
    
    if(dot_labels != fill_label){
      p <- p + ggrepel::geom_text_repel(size=2, max.overlaps=100, fontface = "italic")
    }
    
    if(!is.null(do_fit)){
      p <- p + geom_smooth(method="lm", formula=as.formula("y ~ x"), col="black")
      
      if(do_fit == "log"){
        p <- p + geom_smooth(method="lm", formula=as.formula("y ~ log(x)"), col="blue")
      }
    }
    
    # Change category colors
    if(!is.null(cols)){
      p <- p +
        scale_fill_manual(values=cols) +
        scale_colour_manual(values=cols)
    }
    
    # Tighten axes
    p <- p + 
      scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))
    
    # Format the plot by adding labels and facets and adjusting the format; wite the result to file
    p <- format_and_write_my(p, df, x_label, fill_label, facet_var, filename, labels)
    
    return(p)
  }
}


add_labels_my <- function(in_plot, labels){
  ## Add user defined labels
  ## labels = c(x-label, y-label, plot-title)
  
  base_font <- 9
  
  if(!is.null(labels)){
    in_plot <- in_plot + 
      labs(x=labels[1], y=labels[2], title=labels[3]) + 
      theme(plot.title = element_text(size=base_font),
            axis.title = element_text(size=base_font))
  }
  
  return(in_plot)
}


add_facets2plot_my <- function(in_plot, facet_var, num_row){
  ## Facet a plot by a variable
  ## facet_var - a list with parameters: facet_var[1] = facet_wrap variable; facet_var[1] = ["fixed" or "free"] for scales_va
  
  # Adjust the font size of facets by the number of facets

  base_font <- 10

  in_plot <- in_plot + 
    facet_wrap(as.formula(paste("~", facet_var[1])), nrow=num_row, scales=facet_var[2])  + 
    theme(strip.text = element_text(size = base_font))
  
  return(in_plot)
}


adjust_axes_font_my <- function(in_plot, verbose=FALSE){
  ## Adjust axis text based on the number of labels or remove labels completely
  
  base_font <- 10

  # Find the number of X axis labels
  num_labels <- length(ggplot_build(in_plot)$layout$panel_params[[1]]$x$breaks)

  # Adjust X axis labels
  if(num_labels <= 10){
    font_new <- base_font
    in_plot <- in_plot + theme(axis.text.x = element_text(angle=50, hjust=1, size=font_new))
  } else if(num_labels < 150){
    font_new <- base_font-round(num_labels/20)
    in_plot <- in_plot + theme(axis.text.x = element_text(angle=50, hjust=1, size=font_new))
  } else{
    in_plot <- in_plot + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
  }
  
  # Find the number of Y axis labels
  num_labels <- length(ggplot_build(in_plot)$layout$panel_params[[1]]$y$breaks)

  # Adjust Y axis labels
  if(num_labels <= 6){
    font_new <- base_font
    in_plot <- in_plot + theme(axis.text.y = element_text(angle=0, hjust=1, size=font_new))
  } else if(num_labels < 150){
    font_new <- base_font-round(num_labels/10)
    in_plot <- in_plot + theme(axis.text.y = element_text(angle=0, hjust=1, size=font_new))
  } else{
    in_plot <- in_plot + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
  }
  
  if(verbose == TRUE){
    cat("Number of X labels", num_labels, "\n")
    cat("Number of Y labels", num_labels, "\n")
  }
  
  return(in_plot)
}


adjust_legend_my <- function(in_plot, df, x_label=NULL, fill_label=NULL, verbose=FALSE){
  ## Make a plot look prettier
  
  base_font <- 8

  if(is.null(fill_label) || x_label == fill_label){
    legend_size <- 1
  } else{
    legend_size <- nrow(unique(df[fill_label]))
  }
  
  # Adjust legend font size and shape
  # http://www.sthda.com/english/wiki/ggplot2-legend-easy-steps-to-change-the-position-and-the-appearance-of-a-graph-legend-in-r-software
  if(legend_size <= 1){
    in_plot <- in_plot + 
      theme(legend.position = "none")
  }
  else{
    in_plot <- in_plot + 
      theme(legend.text = element_text(size=base_font), legend.title = element_text(size=base_font), 
            legend.key.size = unit(0.5, "line")) + 
      guides(fill=guide_legend(ncol=1))
  }
  
  if(verbose == TRUE){
    cat("Legend size =", legend_size, "\n")
  }
  
  return(in_plot)
}


write_plot2file_my <- function(in_plot, filename, num_row=1, num_col=1, width=NULL, height=NULL){
  ## Write plot to file
  
  # Derive file name based on the extension
  device <- "png"
  filename_full <- paste(filename, ".", device, sep="")

  if(is.null(width) & is.null(height)){
    # Calculate the size of the figure based on the number of facets
    width <- 5 * ceiling(num_col / 1.5)
    height <- 4 * ceiling(num_row / 1.5)
  }

  # Set max width / height of an image to avoid errors
  if(width > 50){
    width <- 50
  }
  if(height > 50){
    height <- 50
  }
  
  # Write figure to file
  png(file=filename_full, width=width, height=height, units="in", res=300)
  
  dev.control(displaylist="enable")
  show(in_plot)
  dev.off()

  #+fig.width=width, fig.height=height
  show(in_plot)
  
  cat("Figure created:", filename_full, "\n")
}


format_and_write_my <- function(in_plot, df=NULL, x_label=NULL, fill_label=NULL, facet_var=NULL, filename=NULL, labels=NULL, verbose=FALSE){
  ## Format the plot by adding labels and facets and adjusting the format; wite the result to file

  # Make background white
  in_plot <- in_plot + theme_bw()
  
  # Calculate number of facets
  num_row <- 1
  num_col <- 1
  if(!is.null(facet_var)){
    facet_num <- nrow(unique(df[facet_var[1]]))

    # Change the number of columns based on the number of facets
    if(facet_num <= 5){
      num_row <- facet_num
    } else{
      num_row <- ceiling(facet_num / 3)
      num_col <- ceiling(facet_num / num_row)
    }
  }

  if(verbose == TRUE){
    cat("Number of facets =", facet_num, "\n")
  }

  # Add facets
  if(num_row != 1 || num_col != 1){
    in_plot <- add_facets2plot_my(in_plot, facet_var, num_row)
  }
  
  # Adjust axis text based on the number of labels
  in_plot <- adjust_axes_font_my(in_plot)
  
  # Remove legend if not necessary
  in_plot <- adjust_legend_my(in_plot, df, x_label, fill_label)
  
  # Add labels
  if(!is.null(labels)){
    in_plot <- add_labels_my(in_plot, labels)
  }

  # Write to file
  if(!is.null(filename)){
    write_plot2file_my(in_plot, filename, num_row, num_col)
  }
  
  return(in_plot)
}
