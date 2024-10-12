# Author: Andrew Fisher, Anna Lyubetskaya. Date: 22-09-28
# Function for fetching image patch from Concentriq via API call
#
# Assumes that Concentriq API credentials are stored in system

require(httr)
require(tiff)
require(EBImage)


calculate_structure_coord <- function(xy_list, expand_val=100){
  # Calculate the parameters of the structure represented as a rectangle and return a corresponding Concentriq call
  
  SW <- c(min(xy_list[[1]])-expand_val, min(xy_list[[2]])-expand_val)
  rect_width <- max(xy_list[[1]]) - min(xy_list[[1]]) + expand_val*2+1
  rect_height <- max(xy_list[[2]]) - min(xy_list[[2]]) + expand_val*2+1
  
  # Return Concentriq call for the structure
  frame_coord <- paste0(SW[1], ",", 
                        SW[2], ",", 
                        rect_width, ",", 
                        rect_height, "/", 
                        max(rect_width, rect_height), ",/0/default.tiff")
  
  return(frame_coord)
}


calculate_patch_coord <- function(xy_list, swatch_radius_px){
  # Calculate the parameters of the spot represented as a rectangle and return a corresponding Concentriq call
  
  frame_coord <- paste0(in_coord_px[1]-swatch_radius_px, ",",
                        in_coord_px[2]-swatch_radius_px, ",",
                        swatch_radius_px*2+1, ",",
                        swatch_radius_px*2+1, "/",
                        swatch_radius_px*2+1, ",/0/default.tiff")
  
  return(frame_coord)
}


fetch_image_patch <- function(in_concentriq_id, xy_list, in_radius_px, expand_coef = 3, expand_val = 100, 
                              in_rotation = 0, vis_type = "spot"){
  #####DEV#####
  if(FALSE){
    in_concentriq_id <- 993242 # unique Concentriq image ID
    #xy_list <- c(x=9200, y=4800) # c(x,y) of spot center, doesn't need to be named
    xy_list <- data.frame(x=c(26000, 26385, 26780, 25605, 26000, 26385, 26780, 25605, 26000, 26385),
                     y=c(10385, 10600, 10832, 10600, 10832, 11063, 11289, 11063, 11289, 11518))
    # or list of coordinates to form a structure
    in_radius_px <- 100 # radius in pixels of the patch to capture
    in_rotation <- 0 # unused for now, will be needed if coords require rotation
    expand_coef = 3 # multiplied to spot radius for single spot patch size (radius)
    expand_val <- 100 # number of pixels to add to the frame (for structures)
    vis_type <- "structure" # or "spot"
    # "spot" user provides a center spot coordinate of the image and a radius
    # "structure" user provides a data frame of center spot coordinates to be united into a single field of view
  }
  #############
  
  # update the swatch radius
  swatch_radius_px <- in_radius_px * expand_coef
  
  httr::set_config(httr::config(http_version = 0))
  
  # read in Concentriq credentials
  api_user <- Sys.getenv("CONCENTRIQ_API_RO_EMAIL")
  api_pass <- Sys.getenv("CONCENTRIQ_API_RO_PASSWORD")
  
  # construct the initial API call
  base <- 'https://concentriq-prod.rdcloud.bms.com/api/'
  endpoint <- 'images/'
  
  # First query Concentriq to get image id
  call1 <- paste0(base,endpoint,in_concentriq_id)
  get_img <- httr::GET(call1, authenticate(api_user, api_pass, type='basic'), config = httr::config(ssl_verifypeer = FALSE))
  api_return1 <- httr::content(get_img)
  
  # check for error in the request/return
  if(httr::http_error(get_img)) {
    stop(
      sprintf("Concentriq GET image information failed [%s]\n%s\n<%s>", 
              httr::status_code(get_img), 
              api_return1$error$name, 
              api_return1$error$message), 
      call. = FALSE)
  }
  
  # using the iip image server call, stream out the patch data
  iip_img_base <- gsub("info.json","",api_return1$data$imageData$imageSources[[1]]$imageServerUrl)
  
  if(vis_type == "structure"){
    frame_coord <- calculate_structure_coord(xy_list, expand_val=expand_val)
  } else if(vis_type == "spot"){
    frame_coord <- calculate_patch_coord(xy_list, swatch_radius_px)
  }
  
  call2 <- paste0(iip_img_base,frame_coord)
  
  # fetch the patch
  get_patch <- httr::GET(call2, authenticate(api_user, api_pass, type='basic'), config = httr::config(ssl_verifypeer = FALSE))
  
  # check for error in the request/return
  if(httr::http_error(get_img)) {
    stop(
      sprintf("Concentriq GET image patch failed [%s]", 
              httr::status_code(get_img)), 
      call. = FALSE)
  }
  
  # extract the image from the API return
  my_img_mat <- tiff::readTIFF(httr::content(get_patch))
  
  # send image matrix data into an EBImage object
  my_img <- EBImage::Image(my_img_mat, colormode = "Color")
  
  # orientation of image requires transpose
  my_img <- EBImage::transpose(my_img)
  
  # add a spot circle
  if(vis_type == "spot"){
    
    for(a in c(-1, 0, 1)){
      my_img <- EBImage::drawCircle(my_img, swatch_radius_px, swatch_radius_px, radius=in_radius_px + a, col="green")
    }
    
  } else if(vis_type == "structure"){
    
    SW <- c(min(xy_list[[1]])-expand_val, min(xy_list[[2]])-expand_val)
    for(i in 1:nrow(xy_list)){
      for(a in c(-1, 0, 1)){
        # (0;0) is in upper left corner? I need Andrew! ='(
        my_img <- EBImage::drawCircle(my_img, xy_list[[i, 1]]-SW[1], xy_list[[i, 2]]-SW[2],
                                      radius=in_radius_px + a, col="green")
      }
    }
    
  }
  
  return(my_img)
}


combine_and_vis_patches_my <- function(patch_list, title="", filename=NULL, brightness=0, contrast=1,
                                       do_resize=TRUE, nrow=2){
  
  # Identify minimal common swatch dimensions
  swatch_dims <- matrixStats::rowMins(sapply(patch_list, function(x) dim(x)))
  
  # Resize all images down to the minimum one
  if(do_resize == TRUE){
    patch_list <- lapply(patch_list, function(x) EBImage::resize(x, w=swatch_dims[[1]], h=swatch_dims[[2]]))
  }
  
  # Combine patches
  patch_image_list <- EBImage::tile((EBImage::combine(patch_list) + brightness) * contrast, 
                                    nx=ceiling(length(patch_list) / nrow), lwd=nrow, fg.col = "gray")
  
  # Display image or save image to file  
  if(is.null(filename)){
    EBImage::display(patch_image_list, method="raster", all=TRUE, nx=length(patch_list))
    # text(x = 20, y = 20, label = title)
  } else{
    
    EBImage::writeImage(patch_image_list, filename, type="png", quality = 100)
    
    img <- png::readPNG(filename)
    
    # Get image size
    h <- dim(img)[1]
    w <- dim(img)[2]
    
    # Open new file for output
    png(gsub(".png", ".png", filename), width=w, height=h)
    
    par(mar=c(0,0,2,0), xpd=NA, mgp=c(0,0,0), oma=c(0,0,0,0), ann=F)
    plot.new()
    plot.window(0:1, 0:1)
    
    # Fill plot with image
    usr <- par("usr")    
    rasterImage(img, usr[1], usr[3], usr[2], usr[4])
    
    # Add title
    title(title, line=-1, outer=TRUE)
    
    # Close image
    dev.off()
    
  }
  
  return(patch_image_list)
}
