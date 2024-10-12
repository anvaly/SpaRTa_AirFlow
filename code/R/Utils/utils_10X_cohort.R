# Author: Anna Lyubetskaya. Date: 21-01-18
# Functions focused on Seurat data integration
# Note: These functions require large instances minimum!


seurat_merge_my <- function(seurat_list){
  ## Merge a list of Seurat objects
  
  data_combo_seurat <- merge(seurat_list[[1]], y=seurat_list[2:length(seurat_list)], project="Spatial")
  
  return(data_combo_seurat)
}


seurat_merge_loop_my <- function(seurat_list){
  ## Merge a list of Seurat objects
  
  # Remove SCT slot
  for(s in names(seurat_list)){
    Seurat::DefaultAssay(seurat_list[[s]]) <- "Spatial"
    if(keep_SCT == FALSE){
      seurat_list[[s]][["SCT"]] <- NULL
    }
    seurat_list[[s]]@graphs <- list()
    seurat_list[[s]]@neighbors <- list()
    seurat_list[[s]]@reductions <- list()
  }
  
  # Merge Seurat objects
  data_combo_seurat <- seurat_list[[1]]
  # Cleanup to save resources
  seurat_list <- seurat_list[2:length(seurat_list)]
  
  while(length(seurat_list) >= 1){
    print(length(seurat_list))
    
    # Add another Seurat object
    data_combo_seurat <- merge(data_combo_seurat, seurat_list[[1]])
    
    # Cleanup to save resources
    if(length(seurat_list) > 1){
      seurat_list <- seurat_list[2:length(seurat_list)]
    } else{
      seurat_list <- list()
    }
    
    gc()
  }
  
  return(data_combo_seurat)
}


seurat_sct_integrate_my <- function(seurat_list, feature_num=2000, integration_method="cca", reference=NULL){
  ## Use SCT integration approach for merging data
  ## Carefully consider experimental design before integrating
  ## https://satijalab.org/seurat/v3.2/integration.html
  ## Integration type: CCA or RPCA
  
  gc()
  
  # Adjust future.globals.maxSize for method requirements
  options(future.globals.maxSize = 50000*1024^2)
  
  # The size of the smallest dataset to integrate
  size_min <- min(sapply(seurat_list, function(x) ncol(x)))
  
  # Adjust k.weight for integration down if the smallest dataset is below 100
  k.weight <- 100
  if(size_min < 100){
    k.weight <- size_min
  }
  
  # Find variable features
  features <- Seurat::SelectIntegrationFeatures(object.list = seurat_list, nfeatures = feature_num, 
                                                fvf.nfeatures = 10000)
  
  # Calculate Pearson residuals
  seurat_list <- Seurat::PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
  
  # If reference is provided, translate if from a list of sample names to a list of indices from the Seurat list
  if(!is.null(reference)){
    sample_list <- unname(sapply(seurat_list, function(x) unique(x@meta.data$user.Sample_Name)))
    reference <- intersect(reference, sample_list)
    if(!is.null(reference)){
      reference <- sort(c(unname(sapply(reference, function(x) which(grepl(x, sample_list))))))
    } else{
      print("The reference samples couldn't be found in the Seurat object!")
    }
  }
  
  if(integration_method == "cca_sct"){
    
    # Identify SCT anchors
    # Defaults: dims = 1:30, k.anchor = 5, k.filter = 200, k.score = 30 n.trees = 50
    seurat_anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
                                                     anchor.features = features, reference=reference, dims = 1:30,
                                                     k.anchor = 20, k.filter = 200, k.score = 30, n.trees = 50)
    
    # Integrate data
    data_combo_seurat <- Seurat::IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT", k.weight = k.weight,
                                               dims = 1:30)
    
  } else if(integration_method == "rpca_sct"){
    # Apply PCA to individual samples using a common feature set
    seurat_list <- lapply(X = seurat_list, FUN = Seurat::RunPCA, features = features)
    
    # Identify SCT anchors
    seurat_anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
                                                     anchor.features = features, dims = 1:30, 
                                                     reduction = "rpca", k.anchor = 20, reference=reference)
    
    # Integrate data
    data_combo_seurat <- Seurat::IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT", dims=1:30)
  }
  
  return(data_combo_seurat)
}


precast_integrate_my <- function(seurat_list, feature_num=2000, K=15){
  ## Following PRECAST vignette:
  ## https://feiyoung.github.io/PRECAST/articles/PRECAST.BreastCancer.html
  
  
  library(PRECAST)
  
  
  # Format Seurat objects for PRECAST
  for(s in names(seurat_list)){
    # Switch default assay to integrated
    Seurat::DefaultAssay(seurat_list[[s]]) <- "Spatial"
    
    # Put row and col in meta.data for PRECAST
    seurat_list[[s]]@meta.data <- cbind(seurat_list[[s]]@meta.data,
                                        seurat_list[[s]]@images[[s]]@coordinates[match(rownames(seurat_list[[s]]@meta.data),
                                                                                       rownames(seurat_list[[s]]@images[[s]]@coordinates)),])
  }
  
  # Detect number of cores
  core_num <- parallel::detectCores() / 8
  
  # Create PRECASTObject
  PRECASTObj <- PRECAST::CreatePRECASTObject(seurat_list, gene.number = feature_num, selectGenesMethod="HVGs")
  
  # Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting
  PRECASTObj <- PRECAST::AddAdjList(PRECASTObj, platform = "Visium")
  
  # Add a model setting in advance for a PRECASTObj object
  PRECASTObj <- PRECAST::AddParSetting(PRECASTObj, Sigma_equal = FALSE, verbose = TRUE, int.model = NULL,
                                       coreNum = core_num)
  
  # For function PRECAST, users can specify the number of clusters K
  # or set K to be an integer vector by using modified BIC(MBIC) to determine K
  # Given K
  PRECASTObj <- PRECAST::PRECAST(PRECASTObj, K = K)
  
  # Backup the fitting results in resList
  resList <- PRECASTObj@resList
  PRECASTObj <- PRECAST::selectModel(PRECASTObj)
  
  # Integrate data
  seuInt <- PRECAST::IntegrateSpaData(PRECASTObj, species = "human")
  
  # Scale data
  seuInt <- Seurat::ScaleData(seuInt)
  
  # Add tSNE and UMAP representations
  seuInt <- PRECAST::AddTSNE(seuInt, n_comp=2)
  seuInt <- PRECAST::AddUMAP(seuInt, n_comp=2)
  
  # Make sure that reductions are named as expected
  seuInt@reductions$pca <- seuInt@reductions$PRECAST
  seuInt@reductions$tsne <- seuInt@reductions$tSNE
  seuInt@reductions$umap <- seuInt@reductions$UMAP
  
  return(seuInt)
}


seurat_scvi_integration <- function(){
  ## THIS IS AN UNSUCCESSFUL ATTEMPT AT IMPLEMENTING SCVI ON SEURAT ST OBJECT
  
  ## SCVI: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/scvi_in_R.html
  ## Install and initialize python environment for scvi
  
  # system("pip install scvi-tools")
  # devtools::install_github("cellgeni/sceasy")
  
  library(reticulate)
  library(sceasy)
  
  sc <- import("scanpy", convert = FALSE)
  scvi <- import("scvi", convert = FALSE)
  np <- import("numpy", convert = FALSE)
  
  ## Transform data to Python from a merged object
  
  data_combo_seurat <- seurat_merge_my(seurat_list)
  
  # SCVI needs raw counts and can't deal with non-integers
  # data_combo_seurat@assays$Spatial@counts <- as.integer(data_combo_seurat@assays$Spatial@counts)
  
  # Perform SCTransform normalization
  # data_combo_seurat <- Seurat::SCTransform(data_combo_seurat, assay="Spatial", 
  #                                          vars.to.regress = vars_to_regress, variable.features.n=10000,
  #                                          return.only.var.genes = FALSE, verbose = FALSE, vst.flavor = "v2")
  
  # Perform log-norm on the raw slot because tutorial does it
  data_combo_seurat <- Seurat::NormalizeData(data_combo_seurat, assay="Spatial", normalization.method = "LogNormalize", 
                                             scale.factor = 10000)
  
  # Find variable features on raw slot because tutorial does it
  data_combo_seurat <- Seurat::FindVariableFeatures(data_combo_seurat, selection.method = "vst", nfeatures=feature_num,
                                                    assay="Spatial")
  
  # Identify top genes from log-norm most variable genes in the norm slot
  top_genes <- head(Seurat::VariableFeatures(data_combo_seurat, assay="Spatial"), feature_num)
  
  # Prune and simplify the Seurat object
  data_combo_scvi_seurat <- data_combo_seurat[top_genes]
  Seurat::DefaultAssay(data_combo_scvi_seurat) <- "Spatial"
  data_combo_scvi_seurat <- Seurat::DietSeurat(data_combo_scvi_seurat, assays=c("Spatial"))
  data_combo_scvi_seurat@images <- list()
  data_combo_scvi_seurat@meta.data <- data_combo_scvi_seurat@meta.data[c("user.Sample_Name")]
  
  # Convert Seurat to AnnData
  adata <- sceasy::convertFormat(data_combo_scvi_seurat, from="seurat", to="anndata", main_layer="counts", assay="Spatial", 
                                 drop_single_values=FALSE)
  # adata$X = adata$X$astype("int")
  print(adata)
  
  # Run SCVI
  
  # run setup_anndata
  scvi$model$SCVI$setup_anndata(adata)  # batch_key="user.Sample_Name"
  
  # create the model
  model = scvi$model$SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
  
  # train the model
  model$train(max_epochs = as.integer(400))
  
  # get the latent represenation
  latent = model$get_latent_representation()
  
  # put it back in our original Seurat object
  latent <- as.matrix(latent)
  rownames(latent) = colnames(data_combo_seurat)
  data_combo_seurat[["scvi"]] <- Seurat::CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(data_combo_seurat))
  
}
