# Author: Anna Lyubetskaya. Date: 22-02-25

# This script runs comparisons between one target gene and the rest
# The comparison is performed in the following ways:
# - Correlation between target gene and each gene expression profile across spots
# - Linear regression fit: outcome ~ GeneX expression + confounders
# - ElNet feature selection: outcome ~ GeneA expression + ... + GeneZ expression + confounders
# Differences between tissue sections can be accounted by including section identity as a model covariate or by normalizing the distance metric to within the section


## ENVIRONMENT ----


library(Seurat)
library(foreach)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`


source("code/utils/utils_ggplot.R")
source("code/utils/utils_signatures.R")
source("code/utils/utils_lin_regress.R")

source("code/R/Utils/utils_10X_matrix.R")
source("code/R/Utils/utils_10X_signatures.R")


## PARAMETERS ----


# Path to processed Seurat data
sample_name <- "PDAC108_path14_harmonyepi_rpca_sct"
sample_exclude <- NULL

# Distance variable to use for correlations
# If target is a signature, don't forget "sig." prefix
# "PDAC.collisson.classical", "PDAC.moffitt.basal", "PDAC.moffitt.activatedstroma", "BMS.CL.Hypoxia"
# "sig.PDAC.collisson.classical:sig.PDAC.moffitt.activatedstroma", "sig.PDAC.moffitt.basal:sig.PDAC.moffitt.activatedstroma"
target_var <- "sig.PDAC.moffitt.activatedstroma:sig.BMS.CL.Hypoxia"


# Normalize taget gene expression by sample
normalize_dist_by_sample <- FALSE
# Put either target_var or each individual gene in the outcome part of the linear model
outcome_definition <- "gene"  # "gene" or "target_var"
# Add interaction term to the model
interaction_term <- TRUE
# Peform lasso parameter selection
do_reguralize <- TRUE
regularize_type <- "elnet"  # elnet or ridge or lasso

# Variables to include in the ElNet model as confounders
# confounder_vars_num <- paste0("Pathology.", c("IntestineAdj", "MuscleAdj", "NormalAdj",  "LymphNode", "TLSAggregate", "TLSImmature", "TLSMature", 
#                                               "Adipose", "Vessel", "Muscle", "Nerve", "NonEpi", "BenignEpi", "Blood", "ExoEndo", "LuminalNec", "Tumor"), ".percent")
confounder_vars_num <- paste0("sig.", c("PDAC.P19.Acinar","PDAC.P19.Bcell","PDAC.P19.Ductal_1",
                                        "PDAC.collisson.classical", "PDAC.moffitt.basal",  # "PDAC.P19.Ductal_2",
                                        "PDAC.P19.Endocrine","PDAC.P19.Endothelial", "PDAC.moffitt.activatedstroma",
                                        # "PDAC.P19.Fibroblast",
                                        "PDAC.P19.Macrophage","PDAC.P19.Stellate","PDAC.P19.Tcell", 
                                        "BMS.CL.Hypoxia"))
confounder_vars_cat <- c("user.Sample_Name")  # e.g., NULL or c("user.Sample_Name", "user.Tissue")

# Filters to select genes
sct_threshold <- 1
spot_threshold <- 1000

# List of genes to highlight in the scatter plot
gene_highlight_list <- "NONE"  # paste(c("^COL\\d+", "^ITGA\\d+"), collapse="|")  # "NONE" as default

# Add the following signatures to the Seurat object
sig_select <- c("PDAC.P19.Acinar","PDAC.P19.Bcell","PDAC.P19.Ductal_1","PDAC.P19.Ductal_2",
                "PDAC.P19.Endocrine","PDAC.P19.Endothelial","PDAC.P19.Fibroblast",
                "PDAC.P19.Macrophage","PDAC.P19.Stellate","PDAC.P19.Tcell",
                "PDAC.collisson.classical", "PDAC.moffitt.basal",
                "PDAC.moffitt.activatedstroma", "BMS.CL.Hypoxia")  # NULL or c(...)

# Feature # in a spot threshold
feature_threshold <- 2000

# Subset to specific clusters for modeling if needed
cluster_select <- c(3, 10, 0, 4, 6, 9)

# Do a full ElNet
do_elnet <- FALSE


## PATHS ----


# Import additional meta data for the object - neighborhood properties
nb_meta_path <- NULL

# Path to a signature file
sig_path <- "data/import/Signatures/signatures_pdac_20230720.txt"

# Location of pre-processed data
input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, sample_name, "_", target_var, "_", length(c(confounder_vars_num, confounder_vars_cat)), "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)

# Log file
log_file <- paste0(output_path, "log.txt")
write(paste0("Cohort name: ", sample_name), log_file)
write(paste0("Selected clusters: ", paste0(cluster_select, collapse=", ")), log_file, append=TRUE)
write(paste0("UMI threshold = ", feature_threshold), log_file,append=TRUE)
write(paste0("Outcome = ", target_var), log_file, append=TRUE)
write(paste0("Regularization = ", regularize_type), log_file, append=TRUE)


## INGEST ADDITIONAL META DATA ----


# Add neighborhood scores
if(!is.null(nb_meta_path)){
  # Only pick the closest neighborhood
  nb <- 1.5
  
  # Create a wide tibble of neighborhood scores  
  meta_extra_df <- readr::read_delim(nb_meta_path, delim="\t") %>%
    dplyr::filter(NbDistance == nb) %>%
    dplyr::mutate(Feature = paste0(nb, "_", Feature)) %>%
    dplyr::select(-NbDistance) %>%
    df_long2wide_my(rows="Coordinate", cols="Feature", value="ScoreMean")
}


## INGEST DATA ----


# Seurat data
data_seurat <- readRDS(input_path)

# Remove spots if necessary
if(!is.null(feature_threshold)){
  barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data$nFeature_Spatial >= feature_threshold),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# Remove clusters if necessary
if(!is.null(cluster_select)){
  barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data[[data_seurat@misc$user.Clustering]] %in% cluster_select),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# Remove samples if necessary
if(!is.null(sample_exclude)){
  barcode_list <- rownames(data_seurat@meta.data[which(!data_seurat@meta.data[["user.Sample_Name"]] %in% sample_exclude),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}

# Abundant gene list
gene_list <- unique(seurat_select_abundant_genes_my(data_seurat, sct_threshold=0.5, spot_threshold=5, 
                                                    assay="SCT", slot="data"))

# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, gene_list)

# Invert signatures to get an annotated gene list
sig_gene_df <- invert_list_my(signature_list)


# Load signatures and filter them down to only well represented genes
signature_list <- read_filter_signatures_my(sig_path, gene_list, sig_names=sig_select)

# Add signature scores to a Seurat object
if(!is.null(sig_select)){
  data_seurat <- add_signature_scores_my(data_seurat, signature_list, prefix="sig.")
}


## WRANGLE DATA ----


# Abundant gene list
gene_list <- unique(seurat_select_abundant_genes_my(data_seurat, sct_threshold=sct_threshold, spot_threshold=spot_threshold, 
                                                    assay="SCT", slot="data"))

# Exclude MT/Ribo genes for ease of interpretation
gene_list <- setdiff(gene_list, gene_list[grep("^MT-|^RP[SL]|^GM\\d+", gene_list)])

# Make sure that the target variable is in the list
gene_list <- unique(c(gene_list, target_var))

# Add the number of selected genes to the log file
write(paste0("Spot number = ", ncol(data_seurat)), log_file, append=TRUE)
write(paste0("Gene thresholds (SCT threshold / spot number): ", sct_threshold, " / ", spot_threshold), log_file, append=TRUE)
write(paste0("Number of genes selected = ", length(gene_list)), log_file, append=TRUE)


# Subset data to only relevant genes
data_seurat <- subset(data_seurat, features=gene_list)

gc()

# Extract a wide matrix of SCT normalized expression values
# Subset data to only relevant genes
data_wide_df <- tibble::as_tibble(t(as.matrix(Seurat::GetAssayData(data_seurat, assay="SCT", slot="data"))), 
                                  rownames="Coordinate")

# Check that metadata and expression data have the same order out of paranoia
table(rownames(data_seurat@meta.data) == data_wide_df$Coordinate)

# Add the gene of interest as a meta data variable
if(!target_var %in% colnames(data_seurat@meta.data) && target_var %in% colnames(data_wide_df)){
  data_seurat@meta.data[[target_var]] <- data_wide_df[[target_var]]
  
  # Remove the target gene from the expression data
  data_wide_df <- data_wide_df %>%
    dplyr::select(-!!rlang::sym(target_var))
  
  gene_list <- setdiff(gene_list, c(target_var))
}


# Extract the relevant meta data
meta_data <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate")

# Add neighborhood meta data to the expression data
if(!is.null(nb_meta_path)){
  meta_data <- meta_data %>%
    dplyr::inner_join(meta_extra_df, by="Coordinate")
}

# Select meta columns of interest
meta_data <- meta_data %>% 
  dplyr::select(dplyr::all_of(intersect(colnames(meta_data), c("Coordinate", target_var, confounder_vars_num, confounder_vars_cat))))

# Re-normalize pathology distance by sample (% of max)
if(normalize_dist_by_sample == TRUE){
  dist_norm_df <- data_seurat@meta.data[c(target_var, "user.Sample_Name")] %>%
    dplyr::group_by(user.Sample_Name) %>%
    dplyr::mutate(MaxDist = max(!!rlang::sym(target_var)),
                  DISTANCE = !!rlang::sym(target_var) / MaxDist * 100)
  
  meta_data[target_var] <- dist_norm_df$DISTANCE
}

# Add meta data to the expression data
data_wide_df <- data_wide_df %>%
  dplyr::inner_join(meta_data, by="Coordinate")

# Cleanup before running analysis
rm(data_seurat)
gc()

# Calculate feature levels in the selected spots
data_gene_df <- data_wide_df %>%
  df_wide2long_my(key="Feature", val="Score") %>%
  dplyr::mutate(Score = as.numeric(Score)) %>%
  tidyr::drop_na() %>%
  dplyr::group_by(Feature) %>%
  dplyr::summarise(SpotMean = round(mean(Score), 2),
                   SpotCount = sum(ifelse(Score >= 1, 1, 0)))


## SANITY CHECKS BETWEEN VARIABLES ----


if(target_var %in% colnames(data_wide_df)){
  # Scatter of the the outcome variable vs continuous variables
  for(i in confounder_vars_num){
    filename <- paste0(output_path, "scatter_", i)
    create_scatter_plot_my(data_wide_df, x_label=target_var, y_label=i, 
                           fill_label=confounder_vars_cat[1], facet_var=c(confounder_vars_cat[1], "fixed"),
                           size=1, filename=filename, labels=NULL, do_fit=FALSE, stroke=0)
  }
}


## FIND TARGET CORRELATES ----


if(target_var %in% colnames(data_wide_df)){
  # Calculate correlation between the distance and each of the gene expression vectors
  cor_vec <- sapply(gene_list, function(x) 
    if(is.numeric(data_wide_df[[x]])){
      cor(data_wide_df[[x]], meta_data[[target_var]])
    } else{
      0
    })
  names(cor_vec) <- gene_list
  
  # Create a correlation tibble
  cor_df <- tibble::tibble(Symbol = gene_list,
                           R2 = cor_vec) %>%
    dplyr::mutate(R2 = round(R2, 3)) %>%
    dplyr::arrange(R2) %>%
    dplyr::left_join(sig_gene_df, by="Symbol")
  
  # Write correlations to file
  filename <- paste0(output_path, "correlation_result.txt")
  readr::write_delim(cor_df, filename, delim="\t")
  
  # Create a correlation histogram
  filename <- paste0(output_path, "correlation_hist")
  create_hist_plot_my(cor_df, x_label="R2", fill_label="InSignature", intercept=c(-0.25, 0, 0.25), 
                      binwidth=0.01, filename=filename, 
                      labels=c("R2", "Gene number", paste0("Correlation with", target_var)))
}


## INDIVIDUAL LMs FOR EACH GENE ----


## Perform necessary checks on a wide tibble before regression

# Factorize factor fields
confounder_vars_cat <- intersect(colnames(data_wide_df), confounder_vars_cat)
if(!is.null(confounder_vars_cat)){
  data_wide_df[confounder_vars_cat] <- lapply(data_wide_df[confounder_vars_cat], as.factor) 
}

# Explicitly define row names
if(length(rownames(data_wide_df)) != nrow(data_wide_df)){
  data_wide_df <- data_wide_df %>%
    tibble::column_to_rownames("Coordinate")
}
data_wide_df <- data_wide_df %>%
  dplyr::select(-Coordinate)

# Finalize the list of confounders
confounder_list <- intersect(c(confounder_vars_cat, confounder_vars_num), colnames(data_wide_df))


# Save gene LMs to a file as a list of models
filename <- paste0(output_path, "gene_lm_result.rds")
gene_model_list <- list()

# Go through all genes and perform individual fits: Target ~ Gene expression + confounders
if(!file.exists(filename)){
  
  # Define the general lm formula
  formula <- "outcome ~ ."
  
  # Add the target interactive term if necessary
  if(interaction_term == TRUE){
    formula <- paste(formula, "+", target_var)
  }
  
  # Write the formulate to the 
  write(paste0("Formula: ", formula), log_file,append=TRUE)
  write(paste0("Confounders, num: ", paste0(confounder_vars_num, collapse=", ")), log_file, append=TRUE)
  write(paste0("Confounders, cat: ", paste0(confounder_vars_cat, collapse=", ")), log_file, append=TRUE)
  
  
  # Clean the environment
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
  gc()
  
  
  # Parallelize the process for speed
  cl <- parallel::makeCluster(parallel::detectCores() / 4)
  doParallel::registerDoParallel(cl)
  
  gene_model_list <- foreach(gene = gene_list, .combine='c') %dopar% {
    if(gene %in% colnames(data_wide_df)){
      
      # Define the outcome and explicitly define row names
      data_wide_loc_df <- data_wide_df %>%
        dplyr::select(dplyr::all_of(intersect(colnames(data_wide_df), c(target_var, gene, confounder_vars_cat, confounder_vars_num)))) %>%
        tidyr::drop_na()
      
      # Select the outcome variable
      if(outcome_definition == "gene"){
        data_wide_loc_df <- data_wide_loc_df %>%
          dplyr::rename("outcome" = !!rlang::sym(gene))
      } else{
        data_wide_loc_df <- data_wide_loc_df %>%
          dplyr::rename("outcome" = !!rlang::sym(target_var))
      }
      
      
      # Define the general lm formula
      # I have to update formula for every gene because I edit it later
      formula <- "outcome ~ ."
      
      # Add the target interactive term if necessary
      if(interaction_term == TRUE){
        formula <- paste(formula, "+", target_var)
      }
      
      
      # Select parameters through regularization
      if(do_reguralize == TRUE){
        
        # Fit elnet or ridge model
        g <- run_all_models_my(data_wide_loc_df, c(regularize_type), formula=formula)
        
        # Extract zero model coefficients
        coef_remove <- analyze_model_my(g[[regularize_type]], analysis_type="elnet", lambda_type="lambda.1se") %>%
          dplyr::filter(abs(coefficients) == 0) %>%
          dplyr::pull(variable)
        
        # Update the lm input table
        data_wide_loc_df <- data_wide_loc_df %>%
          dplyr::select(dplyr::all_of(setdiff(colnames(data_wide_loc_df), gsub("`", "", coef_remove))))
        
        # Deal with the interaction term if one of its components got removed
        if(interaction_term == TRUE && length(intersect(strsplit(target_var, ":")[[1]], coef_remove)) > 0 ){
          formula <- gsub(paste(" \\+", target_var), "", formula)
        }
      }
      
      # Fit the linear model
      model_res <- stats::lm(formula, data_wide_loc_df)
      out <- list(gene = summary(model_res)$coefficients)
      names(out) <- gene
      
      out
      
    }
  }
  
  parallel::stopCluster(cl)
  
  # Save a list of models to file
  saveRDS(gene_model_list, filename)
} else{
  gene_model_list <- readRDS(filename)
}

# Remove unwanted symbols from model rownames or the downstream code won't work
for(g in names(gene_model_list)){
  rownames(gene_model_list[[g]]) <- gsub("`", "", rownames(gene_model_list[[g]]))
}


## ANALYZE THE RESULT OF INDIVIDUAL LMs FOR EACH GENE ----


## Analyze confounders


if(!is.null(confounder_vars_num)){
  
  coef_list <- list()
  
  # Go through user-defined quantitative confounders
  for(f in confounder_vars_num){
    
    # Analyze the confounder coefficients of the linear model
    coef_gene_df <- t(sapply(names(gene_model_list), function(g)
      if(f %in% rownames(gene_model_list[[g]])){
        gene_model_list[[g]][f, c(1,4)]
      } else{
        c(0, 1)
      }))
    
    # Make sure the column names are captured
    colnames(coef_gene_df) <- c("Estimate", "Pr(>|t|)")
    
    # Wrangle gene coefficients
    coef_gene_df <- coef_gene_df %>%
      tibble::as_tibble(rownames="Symbol") %>%
      dplyr::mutate(Score = round(-log10(`Pr(>|t|)` / length(gene_list)), 1),
                    Score = ifelse(is.infinite(Score), 320, Score))
    
    # Establish the estimate threshold as mean +/- 3SD
    estimate_threshold <- mean(abs(coef_gene_df$Estimate)) + sd(abs(coef_gene_df$Estimate)) * 3
    
    # Define significance
    coef_gene_df <- coef_gene_df %>%
      dplyr::mutate(IsSignificant = `Pr(>|t|)` / length(gene_list) <= 0.05  & abs(Estimate) >= estimate_threshold,
                    Estimate = round(Estimate, 3)) %>%
      dplyr::arrange(desc(IsSignificant), desc(Score)) %>%
      dplyr::select(-`Pr(>|t|)`)
    
    # Update column names
    colnames(coef_gene_df) <- c(colnames(coef_gene_df)[1], paste(colnames(coef_gene_df)[2:4], f, sep="_"))
    
    coef_list[[f]] <- coef_gene_df
  }
  
  # Create a wide table of confounder results
  coef_conf_df <- Reduce(dplyr::inner_join, coef_list)
  
  filename <- paste0(output_path, "gene_lm_result_confounders.txt")
  readr::write_delim(coef_conf_df, filename, delim="\t")
  
}


## Analyze the target


# Analyze the target result of the linear model
if(outcome_definition == "gene"){
  coef_gene_df <- t(sapply(names(gene_model_list), function(g)
    if(target_var %in% rownames(gene_model_list[[g]])){
      gene_model_list[[g]][target_var, c(1,4)]
    } else{
      c(0, 1)
    }))
} else{
  coef_gene_df <- t(sapply(names(gene_model_list), function(g)
    if(g %in% rownames(gene_model_list[[g]])){
      gene_model_list[[g]][g, c(1,4)]
    } else{
      c(0, 1)
    }))
}

# Make sure the column names are captured
colnames(coef_gene_df) <- c("Estimate", "Pr(>|t|)")

# Wrangle gene coefficients
coef_gene_df <- coef_gene_df %>%
  tibble::as_tibble(rownames="Symbol") %>%
  dplyr::mutate(Score = -log10(`Pr(>|t|)` / length(gene_list)),
                Score = ifelse(is.infinite(Score), 320, Score))

# Establish the estimate threshold as mean +/- 3SD
estimate_threshold <- mean(coef_gene_df$Estimate) + sd(coef_gene_df$Estimate) * 3

# Annotate and filter model results
if(estimate_threshold > 0){
  
  # Call significant coefficients
  coef_gene_df <- coef_gene_df %>%
    dplyr::mutate(IsSignificant = `Pr(>|t|)` / length(gene_list) <= 0.05  & abs(Estimate) >= estimate_threshold) %>%
    dplyr::arrange(desc(IsSignificant), desc(Score))
  
  # Establish a threshold for labeling genes on a plot
  score_threshold <- mean(coef_gene_df$Score) + sd(coef_gene_df$Score) * 5
  
  # Add labels for plotting and signature annotations
  coef_gene_df <- coef_gene_df %>%
    dplyr::mutate(Labels = ifelse(grepl(gene_highlight_list, Symbol), Symbol, ""),
                  Labels2Plot = ifelse(abs(Estimate) >= estimate_threshold & Score >= score_threshold, 
                                       Symbol, "")) %>%
    dplyr::left_join(sig_gene_df, by="Symbol")
  
  # Visualize the result of linear modeling of gene expression ~ distance to center
  filename <- paste0(output_path, "gene_lm_result_scatter")
  create_scatter_plot_my(coef_gene_df, x_label="Estimate", y_label="Score", 
                         fill_label="IsSignificant", shape=21, size=1, dot_labels="Labels2Plot", 
                         filename=filename, labels=NULL, do_fit=NULL, stroke=0)
  
  # Annotate and format lm results for output
  coef_gene_ann_df <- coef_gene_df %>%
    dplyr::select(-`Pr(>|t|)`, -InSignature, -Labels2Plot, -Num_Sigs) %>%
    dplyr::mutate(Estimate = round(Estimate, 2),
                  Score = round(Score, 1)) %>%
    dplyr::arrange(desc(IsSignificant), desc(Score), desc(Estimate)) %>%
    dplyr::left_join(data_gene_df, by=c("Symbol" = "Feature"))
  
  
  # Add the information about quantitative confounders of the model
  if(!is.null(confounder_vars_num)){
    coef_gene_ann_df <- coef_gene_ann_df %>%
      dplyr::left_join(coef_conf_df, by="Symbol") %>%
      dplyr::mutate_if(is.numeric, round, 3)
  }
  
  # Write to file
  filename <- paste0(output_path, "gene_lm_result_scatter.txt")
  readr::write_delim(coef_gene_ann_df, filename, delim="\t")
  
  # Select genes of interest for targeted follow up
  coef_gene_ann_filt_df <- coef_gene_ann_df %>%
    dplyr::filter(IsSignificant == TRUE) %>%
    dplyr::arrange(desc(Estimate))
  
  filename <- paste0(output_path, "gene_lm_result_ann_scatter.txt")
  readr::write_delim(coef_gene_ann_filt_df, filename, delim="\t")
  
}


## ELNET OF TARGET ON EXPRESSION ----


if(do_elnet == TRUE){
  
  # Clean the environment
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
  gc()
  
  
  # Define the outcome and explicitly define row names
  data_wide_loc_df <- data_wide_df %>% 
    dplyr::rename("outcome" = !!rlang::sym(target_var))
  
  # Run ElNet
  filename <- paste0(output_path, "dist_elnet_result.rds")
  if(!file.exists(filename)){
    model_res <- run_all_models_my(data_wide_loc_df, c("elnet"), output_loc=NULL, formula=NULL)
    saveRDS(model_res, filename)
  } else{
    model_res <- readRDS(filename)
  }
  
  model_res$elnet
  
  # Analyze feature importance
  # https://topepo.github.io/caret/variable-importance.html
  elnet_imp <- caret::varImp(model_res$elnet)
  
  # Analyze the result of ElNet
  coef_df <- analyze_model_my(model_res$elnet, analysis_type="elnet", lambda_type="lambda.1se") %>%
    dplyr::left_join(tibble::as_tibble(elnet_imp$importance, rownames="variable"), by="variable") %>%
    dplyr::left_join(sig_gene_df, by=c("variable" = "Symbol")) %>%
    dplyr::mutate(Overall = round(Overall, 2))
  
  # Number of top features to visualize
  n <- 30
  
  # Visualize the result of ElNet
  filename <- paste0(output_path, "dist_elnet_result")
  coef_df <- visualize_model_my(coef_df, filename, feature_list=gene_list, n=n)
  
  # Write correlations to file
  filename <- paste0(output_path, "dist_elnet_result_coefficients_filt.txt")
  readr::write_delim(coef_df %>%
                       dplyr::arrange(-Overall) %>%
                       dplyr::filter(Overall >= 5), filename, delim="\t")
  
  # Plot importance
  filename <- paste0(output_path, "dist_elnet_result_imp_bar")
  p <- create_bar_plot_my(coef_df %>%
                            dplyr::filter(Overall >= 25), x_label="variable", y_label="Overall", 
                          fill_label="IsLineage", filename=filename, reorder_x=TRUE, 
                          labels=c("Feature", "Importance", "ElNet, coefficient importance"))
  
  # Default feature imporance visualization
  filename <- paste0(output_path, "dist_elnet_result_imp_line.png")
  png(filename, width = 10, height = 5, units = "in", res = 300)
  plot(elnet_imp, top = n)
  dev.off()
  
}
