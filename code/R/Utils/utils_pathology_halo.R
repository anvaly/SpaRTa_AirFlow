# Author: Andrew Fisher. Date: 22-08-23


library(doParallel)
source("code/R/Utils/utils_10X_vis.R")


if(packageVersion("spatstat") < 1.65){
  update.packages("spatstat")
}


xml_inspect_regions_my <- function(xml_file){
  ## Read an XML file and list all region names
  
  print(xml_file)
  
  # open annot file
  fp <- xml2::read_xml(xml_file)
  
  # read annotation object
  annots <- xml2::xml_find_all(fp, './/Annotation')
  
  # Return annotation names
  return(xml2::xml_attr(annots,"Name"))
}


spot_annotation_overlap_my <- function(spot_x, spot_y, spot_radius, annotation_win, do_mask=FALSE){
  ## Computing overlap area between spots and annotations
  
  ## Add parallelization to handle large vectors
  cl <- parallel::makeCluster(parallel::detectCores() / 2)
  doParallel::registerDoParallel(cl)
  
  # Create a mask of the annotation
  if(do_mask == TRUE){
    annotation_win <- spatstat::as.mask(annotation_win, dimyx=c(1000, 1000))
  }
  
  spot_annotation_list <- foreach::foreach(i = 1:length(spot_x), .combine="c") %dopar% {
    # Calculate the overlap between each spot approximated as a disc with an polygon of a pathology-defined compartment
    spatstat.geom::area(spatstat.geom::intersect.owin(
      spatstat.geom::disc(radius = spot_radius, centre = c(spot_x[i], spot_y[i]), npoly = 256), 
      annotation_win))
  };
  
  parallel::stopCluster(cl)
  
  gc()
  
  return(spot_annotation_list)
}


spot_pathology_distance_my <- function(spot_x, spot_y, spot_radius, anno_polys, anno_class, do_fill=FALSE, make_plots=FALSE, plot_path=NULL, plot_name=NULL, seurat_img=NULL){
  ## Compute the distance from spots to annotations
  
  require(pbmcapply)
  
  ########
  # dev
  if(FALSE){
    spot_x <- X
    spot_y <- Y
    anno_polys <- polys
    cl <- anno_class
    do_fill=(cl %in% do_fill_classes)
    make_plots = TRUE
    plot_path = file.path(output_init_figs,"DigiPath")
    plot_name = paste0("Annotation_Distance - ",sample_name," - ",cl)
    seurat_img = GetImage(data_seurat, m = "raster")
  }
  ########
  
  # grab the relevant annotations
  annotation_win <- polys[[cl]]
  
  # fill holes of annotation if necessary
  if(do_fill == TRUE){
    polys_temp <- annotation_win
    
    # Modify the bdry list (list of polygon shapes as xy coordinates) by logically filtering out the holes
    polys_temp$bdry <- polys_temp$bdry[!as.logical(lapply(polys_temp$bdry, spatstat.utils::is.hole.xypolygon))]
    
    # Split the owin object into separate shapes
    shape_list <- purrr::map(polys_temp$bdry, ~ spatstat::owin(poly = .x))
    # Find the shape union
    polys_temp <- spatstat::union.owin(spatstat::as.solist(shape_list))
    
    annotation_win <- polys_temp
  }
  
  # Distance functions
  dist_fun_away <- spatstat::distfun.owin(annotation_win, invert = FALSE)
  dist_fun_into <- spatstat::distfun.owin(annotation_win, invert = TRUE)
  
  # Report a summary for the computed distance function
  #summary(dist_fun)
  
  # Visual representation of the pathology feature to which we are measuring distance
  if(make_plots){
    # ensure necessary parameters passed for plotting
    if(is.null(plot_path) | is.null(plot_name) | is.null(seurat_img)){
      warning("Unable to generate annotation distance plots; no path, seurat image, and/or filename was given to the function.")
    } else {
      # ensure output directory exists
      if (!dir.exists(plot_path)){
        dir.create(plot_path, recursive = TRUE)
      }
      png(file.path(plot_path, paste0(plot_name,".png")), width = 1000, height = 800)
      par(mfrow=c(2,2), mar = c(0, 0, 2, 1), oma = c(0,0,2,0))
      plot(c(0,1),c(0,1), type = "n", xlab = "", ylab = "", main = "", axes = FALSE, asp = 1)
      mtext("H&E", side = 3)
      rasterImage(seurat_img, 0, 1, 1, 0)
      plot(polys[["Tissue"]], main="")
      plot(annotation_win, add=TRUE, col='red')
      mtext(paste0(cl," - do_fill=",do_fill), side = 3)
      if(cl=="Tissue"){
        plot(dist_fun_away, eps=100, main="")
      }else{
        plot(dist_fun_away, W = polys[["Tissue"]], eps=100, main="")
      }
      mtext("Distance away from region", side = 3)
      plot(dist_fun_into, W = polys[["Tissue"]], eps=100, main="")
      mtext("Distance into region", side = 3)
      mtext(plot_name, side = 3, line=0, outer = TRUE, cex=2)
      dev.off()
    }
  }
  
  # Calculate distance between every spot center and the pathology feature
  # Parallelize to save time, keep a couple cores open to avoid crashes
  #############################
  # define function to parallelize
  dist_comp_f <- function(x_in, y_in, dist_in, dist_out, spot_radius){
    # Identify X and Y ranges for this spot
    x_list <- (x_in - spot_radius) : (x_in + spot_radius)
    y_list <- (y_in - spot_radius) : (y_in + spot_radius)
    
    # Create a grid of all coordinate values as a square around the spot center
    square_coord <- expand.grid(x_list, y_list)
    
    # Filter the square down to circle
    circle_coord <- square_coord[which(sqrt((square_coord$Var1 - x_in)^2 + (square_coord$Var2 - y_in)^2) <= spot_radius),]
    
    # Calculate mean distances across the spot
    away_dist_val <- mean(dist_out(list(x=circle_coord$Var1, y=circle_coord$Var2)))
    into_dist_val <- -1*mean(dist_in(list(x=circle_coord$Var1, y=circle_coord$Var2)))
    
    # Choose which distance to use; unless 0, will use the away
    if (away_dist_val==0){
      into_dist_val
    } else {
      away_dist_val
    }
  }
  
  distance_list <- pbmcmapply(dist_comp_f, spot_x, spot_y,
                              MoreArgs = list(dist_in=dist_fun_into, dist_out=dist_fun_away, spot_radius=spot_radius),
                              mc.preschedule = TRUE,
                              mc.cores = max(1, parallel::detectCores()-3))
  
  return(distance_list)
}


plot_pathology_ingestion_check_my <- function(polys, spot_x, spot_y, spot_radius, filename, main_class="Tissue", title="", cols=NULL, class_list){
  ## Plot classification windows and 10X spots from the Seurat object
  ## https://stackoverflow.com/questions/8929663/r-legend-placement-in-a-plot
  
  # Write plot to a file
  png(filename, width=6, height=4, units="in", res = 300)
  
  # Put the legend on the right side
  par(mar = c(0, 0, 1, 0))
  layout(matrix(c(1, 2), ncol=2), widths=c(4, 2), heights=c(4, 4))
  
  # Leave only classes present in the list of polygons
  class_list <- intersect(class_list, names(polys))
  
  # Establish color scheme if not provided by the user
  if(is.null(cols)){
    cols <- define_cols_my(n=length(class_list), col_type="jet")
    names(cols) <- class_list
  }
  
  # Plot the tissue outline
  if(!main_class %in% names(polys)){
    main_class <- names(polys)[1]
  }
  spatstat.geom::plot.owin(polys[[main_class]], ylim=rev(polys[[main_class]]$yrange), border=cols[main_class], main=title)
  
  # Plot the rest of the pathology classes
  for(i in 1:length(class_list)){
    cl <- class_list[i]
    if(cl != main_class){
      spatstat.geom::plot.owin(polys[[cl]], ylim=rev(polys[[main_class]]$yrange), border=cols[[cl]], col=cols[[cl]], add=TRUE)
    }
  }
  
  # Add Visium spots
  for(i in 1:length(X)){
    # Draw a disc around spot centers, find radius in the 10X Seurat object
    spot_disc_ppp <- spatstat.geom::disc(radius = spot_radius, centre = c(spot_x[i], spot_y[i]), npoly = 256)
    
    # Add the disc to the plot
    spatstat.geom::plot.owin(spot_disc_ppp, ylim=rev(spot_disc_ppp$yrange), border="black", lwd=0.25, add=TRUE)
  }
  
  # Put legend in the second "subplot" to the right of the image
  plot.new()
  
  # Add legend to the plot
  legend(x = "topleft", legend = class_list, col = cols[class_list], lty = 1, lwd = 2, cex = 0.5)
  
  # Save plot from device to an object
  # p <- recordPlot()
  
  dev.off()
  
  # return(p)
}


halo_annot_to_poly_my <- function(xml_file, coef=NULL, mpp=1, as_owin=TRUE, region_list=NULL,
                                  attr_dict=NULL) {
  ## This script transforms a HALO XML output file to a polygon
  
  ## Inputs:
  ##   xml_file: path to Halo .annotations file
  ##   mpp: Scaling factor to convert vertices from pixels to microns
  ##   as_owin: If TRUE, output will be a spatstat owin object. Otherwise, list(s) of vertices will be returned.
  ##   region_list (optional): list of names, if given only annotations with names matching those in the list will be included
  
  # open annot file
  fp <- xml2::read_xml(xml_file)
  
  # read annotation object
  annots <- xml2::xml_find_all(fp, './/Annotation')
  
  # filter nodeset based on input region_list
  if (!is.null(region_list)){
    # filter down nodeset based on region names
    annots <- annots[is.element(xml2::xml_attr(annots,"Name"), region_list)]
  }
  
  annots_attrs <- purrr::map(annots, xml2::xml_attrs)
  
  # read region objects
  regions <- purrr::map(annots, ~ xml2::xml_find_all(.x, './/Region'))
  region_attrs <- purrr::map(regions, xml2::xml_attrs)
  
  # Update attribute dictionary if user wants to collapse any of the classes
  if(!is.null(attr_dict)){
  for(i in 1:length(annots_attrs)){
    if(annots_attrs[[i]][["Name"]] %in% names(attr_dict)){
      prev_name <- annots_attrs[[i]][["Name"]]
      annots_attrs[[i]][["Name"]] <- attr_dict[[prev_name]]
    }
  }
  }

  # if an additional scaling factor provided, update mpp
  if(!is.null(coef)){
    mpp <- mpp * coef
  }
  
  # read vertex (V) objects for polygons
  vertices <- regions %>%
    purrr::map_depth(2, xml2::xml_find_all, './/V') %>%
    purrr::map_depth(3, xml2::xml_attrs, c('X', 'Y')) %>%
    purrr::map_depth(2, purrr::transpose) %>%
    purrr::map_depth(3, ~ as.numeric(.x) * mpp) %>%
    purrr::map_depth(2, rlang::set_names, c('x', 'y'))
  
  # Check for ellipses and rectangles in set of regions, convert to polys as needed
  if (any(grepl("Ellipse|Rectangle",region_attrs))){
    for (i in 1:length(region_attrs)){
      if (grepl("Ellipse|Rectangle",region_attrs[i])){
        for (j in 1:length(region_attrs[[i]])){
          if (region_attrs[[i]][[j]][["Type"]]=="Ellipse"){
            # Change the vertices in related list from bbox to poly approx counter-clockwise
            ellipse_xy <- shape::getellipse(rx = abs(diff(vertices[[i]][[j]]$x)),
                                            ry = abs(diff(vertices[[i]][[j]]$y)),
                                            mid = c(mean(vertices[[i]][[j]]$x),vertices[[i]][[j]]$y),
                                            dr = 0.01)
            vertices[[i]][[j]]$x <- ellipse_xy[,1]
            vertices[[i]][[j]]$y <- ellipse_xy[,2]
          } else if (region_attrs[[i]][[j]][["Type"]]=="Rectangle"){
            # Expand vertices list from 2pts to 4pts counter-clockwise
            rect_x <- vertices[[i]][[j]]$x
            rect_y <- vertices[[i]][[j]]$y
            vertices[[i]][[j]]$x <- c(max(rect_x),min(rect_x),min(rect_x),max(rect_x))
            vertices[[i]][[j]]$y <- c(min(rect_y),min(rect_y),max(rect_y),max(rect_y))
          }
        }
      }
    }
  }
  
  # Check for any pins in annotation regions and remove section (cannot make a window from points!)
  bool_filt <- !grepl("Pin|Ruler",region_attrs)
  vertices <- vertices[bool_filt]
  annots_attrs <- annots_attrs[bool_filt]
  region_attrs <- region_attrs[bool_filt]
  
  # convert vertices into polygons
  # NOTE: currently assuming all regions are polygon type
  polys <- purrr::map2(vertices, region_attrs, purrr::map2, function(xy, attrib) {
    if(any(grepl("IsNegative", region_attrs))){
      is_hole <- as.numeric(attrib[["IsNegative"]] == 1)
    }
    else if(any(grepl("NegativeROA", region_attrs))){
      is_hole <- as.numeric(attrib[["NegativeROA"]] == 1)
    }
    else{
      warning("Did not see attribute annotation for holes!")
    }
    # is_hole <- (as.numeric(attrib[["IsNegative"]]) == 1 | as.numeric(attrib[['NegativeROA']]) == 1)
    is_hole_ordered <- spatstat.utils::is.hole.xypolygon(xy)
    
    # make sure vertices are ordered properly
    if(xor(is_hole, is_hole_ordered))
      return(purrr::map(xy, rev))
    else
      return(xy)
  })
  
  # convert from list to spatstat owin
  if(as_owin){
    polys <- purrr::map(polys, ~ spatstat.geom::owin(poly = .x))
  }
  
  # set same layer names as annotation file
  names(polys) <- purrr::map(annots_attrs, ~ .x[['Name']])
  
  return(polys)
}


halo_to_ppp_my <- function(in_ImgObj, in_Anno=NULL, in_mpp=NULL, coreg_prefix="", in_anno_str="Tissue",
                           drop_marks=NULL, add_marks=NULL) {
  ## This script transforms a HALO nuclei segmentation output file to a spatstat ppp object
  
  ## Inputs:
  ##   in_ImgObj: path to Halo image object export csv file
  ##   in_Anno: path to Halo/Aperio exported annotation file to fetch ppp window boundary from, 
  ##            if NULL is given, then the ppp window will default to a bounding rectangle of all the points
  ##   in_mpp: Scaling factor to convert vertices from pixels to microns (if left NULL, keeps units in pixels)
  ##   coreg_prefix: if coregistered/custom coordinates are to be used instead of default 
  ##                 "XMin","XMax", etc., then give the prefix for the column name here. For example,
  ##                 "Stain1.Registered." will use the coregistered coordinates of the Stain1 section (HALO export option)
  ##   in_anno_str: name of the annotation region to be read in as the ppp window boundary (default is "Tissue")
  ##   drop_marks: vector of column names to drop from marks (to save storage space)
  ##   add_marks: named vector of mark columns to add to ppp
  
  
  if (is.null(in_mpp)) {
    in_mpp <- 1 # Default value
    ppp_units <- "px"
  } else {
    ppp_units <- "um"
  }
  
  read.func = data.table::fread
  if (stringr::str_ends(in_ImgObj, 'gz')) {
    read.func = function(f) data.table::fread(cmd=sprintf("gzcat '%s'", f))
  }
  
  # Read in ImgObj Table with specified coordinate columns
  df_ImgObj = purrr::map(in_ImgObj, read.func) %>%
    as.data.frame() %>%
    dplyr::bind_rows() %>%
    tibble::as_tibble(.name_repair = 'universal') %>%
    dplyr::mutate(X_center = 
                    dplyr::select(., matches(paste0("^(",coreg_prefix,"XMin|",coreg_prefix,"XMax)$"))) %>% 
                    # use pmap on this subset to get a vector of min from each row.
                    # dataframe is a list so pmap works on each element of the list that is to say each row
                    purrr::pmap_dbl(sum)
    ) %>%
    dplyr::mutate(X_center = X_center/2*in_mpp) %>%
    dplyr::mutate(Y_center = 
                    dplyr::select(., matches(paste0("^(",coreg_prefix,"YMin|",coreg_prefix,"YMax)$"))) %>% 
                    # use pmap on this subset to get a vector of min from each row.
                    # dataframe is a list so pmap works on each element of the list that is to say each row
                    purrr::pmap_dbl(sum)
    ) %>%
    dplyr::mutate(Y_center = Y_center/2*in_mpp) %>%
    dplyr::select(X_center, Y_center, everything())
  
  if (! is.null(drop_marks)) {
    for (n in drop_marks) df_ImgObj[[n]] = NULL
  }
  
  if (! is.null(add_marks)) {
    for (n in names(add_marks)) {
      df_ImgObj[[n]] = add_marks[n]
      if (is.character(df_ImgObj[[n]])) df_ImgObj[[n]] = as.factor(df_ImgObj[[n]])
    }
  }
  
  if (!is.null(in_Anno)){
    # Import annotation
    anno_owin <- halo_annot_to_poly(in_Anno, mpp = in_mpp, as_owin = TRUE, region_list = list(in_anno_str))
    my_win <- spatstat.geom::owin(anno_owin[[in_anno_str]])
  } else {
    # If no tissue annotations are used, then simply create a bounding rectangle for ppp window
    my_win <- spatstat.geom::owin(range(df_ImgObj$X_center), range(df_ImgObj$Y_center))
  }
  
  spat_ppp <- spatstat.geom::as.ppp(df_ImgObj, W = my_win, SIMPLIFY=FALSE)
  spatstat.geom::unitname(spat_ppp) <- ppp_units # assign units (px or um)
  
  return(spat_ppp)
}
