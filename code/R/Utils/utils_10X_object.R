# Author: Anna Lyubetskaya. Date: 20-11-28


seurat_misc_to_meta <- function(seurat_list, misc_field_list=NULL){
  ## Move misc parameters to meta data slot
  
  # misc_field_list <- c("user.Sample_Name", "user.Group_Name", "user.Tissue", "user.Organism", 
  #                      "user.Protocol", "user.Date", "user.Area", "user.Project_ID", 
  #                      "user.pt.size.factor", "user.Primary")
  
  # Find all misc slots in the first object
  if(is.null(misc_field_list)){
    misc_field_list <- names(seurat_list[[1]]@misc)
  }
  
  # List of misc fields to keep
  if(!is.null(misc_field_list)){
    
    # Add sample names to barcode names and add batch variable to the meta data
    for(s in names(seurat_list)){
      # seurat_list[[s]] <- Seurat::RenameCells(seurat_list[[s]], add.cell.id=s)
      # seurat_list[[s]]@meta.data[["Batch"]] <- s
      
      # Add misc values to meta data to keep them through integration
      for(m in misc_field_list){
        if(m %in% names(seurat_list[[s]]@misc)){
          seurat_list[[s]]@meta.data[[m]] <- seurat_list[[s]]@misc[[m]]
        }
      }
    }
  }
  
  return(seurat_list)
}


add_misc_to_seurat_object_my <- function(data_seurat, field="pt.size.factor", dict_value=NULL, default_value=1.7){
  ## Add a misceleneous parameter to a Seurat object
  
  if(!is.null(dict_value)){
    data_seurat@misc[[field]] <- dict_value
  } else{
    print("Default value added")
    data_seurat@misc[[field]] <- default_value
  }
  
  return(data_seurat)
}
