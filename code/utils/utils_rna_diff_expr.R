# Author: Anna Lyubetskaya. Date: 19-10-10
# Perform the differential expression analysis using DESeq2
# Libraries used: DeSeq2, MASS, sfsmisc

# Published on 10/08/2019: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# RNA-seq analysis for detecting quantitative trait-associated genes: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4829873/
# https://data.princeton.edu/R/linearModels


source("code/utils/utils_ggplot.R")
source("code/utils/utils_in_out.R")


edgeR_normalize_my <- function(eset){
  ## Normalize an eset
  
  # Create a DGEList Object for edgeR
  dge <- edgeR::DGEList(counts=Biobase::exprs(eset), genes = rownames(eset))
  # Calculate normalization factors
  dge <- edgeR::calcNormFactors(dge, method="TMM")
  # Get counts per million.  
  Biobase::exprs(eset) <- edgeR::cpm(dge, log = T)
  
  return(eset)
}


dea_from_wide_tibble_my <- function(expr_df, meta_df, params, name_out="~/dea", do_batch=FALSE){
  ## Filter rows and columns of the data matrix and perform PCA
  ## Rows = features, columns = samples
  
  ## params <- list(sample = "DepMap_ID", 
  ##                sample_filter = c(),
  ##                feature = "Symbol", 
  ##                feature_filter = c(),
  ##                batch = "lineage",
  ##                condition = "IsMutated",
  ##                contrasts = list(c("IsMutated", "TRUE", "FALSE")),
  ##                num_cores = 1)
  
  # Filter raw counts in the wide format for select genes and samples for the DEA
  # Explicitly name rows from a select column: necessary for DEA
  data_select <- expr_df %>%
    dplyr::filter(!!sym(params$feature) %in% params$feature_filter) %>%
    dplyr::select(dplyr::all_of(c(params$feature, params$sample_filter))) %>%
    tibble::column_to_rownames(var = params$feature)
  
  # Sort data by column names
  data_select <- data_select[, order(colnames(data_select))]
  
  # Extract only DEA-relevant sample information
  meta_select <- meta_df %>%
    # dplyr::group_by_at(c(params$sample, params$condition, params$batch)) %>%
    dplyr::filter(!!sym(params$sample) %in% params$sample_filter) %>% 
    dplyr::select(dplyr::all_of(c(params$sample, params$condition, params$batch))) %>%
    dplyr::arrange(!!sym(params$sample)) %>%
    tibble::column_to_rownames(var = params$sample)
  
  # Count number of samples per group
  group_summary <- meta_select %>% 
    dplyr::group_by_at(vars(c(params$condition, params$batch))) %>% 
    dplyr::tally()
  
  # Write number of samples per group to file
  readr::write_delim(group_summary, path=paste0(name_out, "_groups.csv"), delim = "\t")
  
  
  # Perform DEA using DESeq2
  filename_rda <- paste0(name_out, ".rda")
  if(file.exists(filename_rda)){
    load(file=filename_rda)  # Loads DESeq2 Data Set called dds
  } else{
    dds <- dea_model_my(data_select, meta_select, params[c("condition", "batch", "num_cores")], filename_rda, do_batch=do_batch)
  }
  
  print(DESeq2::resultsNames(dds))
  
  # Extract contrast-specifc information; analyze the result of every comparison of interest
  data_contrasts <- list()
  for(contrast in params$contrasts){
    contrast_string <- paste(contrast, collapse="_")
    filename_contrast <- paste0(name_out, "_", contrast_string)
    
    cat("Contrast =", contrast_string, "\n")
    
    if(file.exists(paste0(filename_contrast, ".tsv"))){
      data_contrasts[[contrast_string]] <- read_file2df_my(paste0(filename_contrast, ".tsv"), delim="\t", in_cols=NULL)
    } else{
      data_contrasts[[contrast_string]] <- dea_compare_contrasts_my(dds, contrast, filename_contrast, 
                                                                    significance_threshold=params$pv_threshold, 
                                                                    foldchange_threshold=params$fc_threshold)
    }
  }
  
  return(list("data" = data_select, 
              "meta" = meta_select, 
              "dds" = dds,
              "contrasts" = data_contrasts))
}


dea_model_my <- function(data_df, sample_info, params, name_out, do_batch=FALSE){
  ## Perform DEA using DESeq2
  
  ## params <- list(batch = "lineage", 
  ##                condition = "IsMutated",
  ##                contrasts = list(c("IsMutated", "TRUE", "FALSE")),
  ##                num_core = 1)
  ## contrasts are passed in the form "condition_B_v_A"
  
  BiocParallel::register(BiocParallel::MulticoreParam(params$num_core))
  
  cat("Performing DEA: condition =", params$condition, "batch =", params$batch, "\n")
  
  # Values in the input matrix have to be raw counts (integers) and columns of the countData should be same (and in the same order) as rows in the colData
  data_matrix <- data.matrix(data_df)
  mode(data_matrix) <- "integer"
  
  # Make contrast variables into factors
  for(p in params$condition){
    sample_info[[p]] <- as.factor(sample_info[[p]])
  }
  
  if(is.null(params$batch) || do_batch == FALSE){
    formula <- as.formula(paste("~", paste(params$condition, collapse=" + ")))
  } else{
    formula <- as.formula(paste("~", params$batch, "+", paste(params$condition, collapse=" + ")))
  }
  
  # Create a DESeq2 object
  # design: to benefit from the default settings, put the variable of interest at the end of the formula and make sure the control level is the first level
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = data_matrix,
                                        colData = sample_info,
                                        design = formula)
  
  # Perform DEA analysis
  dds <- DESeq2::DESeq(dds, parallel=TRUE)
  # Write the entire DESeq2 data structure to an RDA object
  save(dds, file=name_out)
  
  return(dds)
}


dea_compare_contrasts_my <- function(dds, contrast, name_out, significance_threshold=0.05, foldchange_threshold=1){
  ## Using a DESeq2 data structure, extract and compare specific pairs of variables
  ## Write a human readable output and create plots

  # Extract a contrast
  data_contrast <- DESeq2::results(dds, contrast=contrast, tidy=TRUE) %>%
    tibble::as_tibble() %>%
    dplyr::rename("Symbol" = "row") %>%
    dplyr::arrange(pvalue)
  
  # Write individual contrast to a file
  # readr::write_delim(data_contrast, path=paste0(name_out, ".tsv"), delim="\t")
  
  data_contrast <- data_contrast %>%
    dplyr::mutate(neglog10padj = round(-log10(padj), 2), 
                  log2baseMean = round(log2(baseMean), 2),
                  IsSignificant = padj <= significance_threshold & abs(log2FoldChange) >= foldchange_threshold,
                  log2FoldChange = round(log2FoldChange, 2),
                  Direction = ifelse(log2FoldChange > 0, "UP", "DN"),
                  Contrast = paste0(contrast[1], ":", contrast[2], "/", contrast[3]),
                  stat = round(stat, 2)) %>%
    dplyr::select(Symbol, IsSignificant, Direction, log2FoldChange, neglog10padj, log2baseMean, stat, Contrast) %>%
    dplyr::arrange(desc(IsSignificant), desc(neglog10padj))
  
  data_contrast[["IsSignificant"]] <- factor(data_contrast[["IsSignificant"]], levels=c(TRUE, FALSE))
  data_contrast[["Direction"]] <- factor(data_contrast[["Direction"]], levels=c("UP", "DN"))
  
  ## This takes a lot of time (due to lfcShrink) and I am not sure it adds a lot of value
  # Shrink log fold changes associated with a given condition
  #res_shrink <- DESeq2::lfcShrink(dds, contrast=contrast, type="normal")  # Types: apeglm, ashr, normal
  
  # Plot DEA result: FC v adjusted p-value and mean counts v FC
  # Plotting shrunken log2 fold changes removes noise from low count genes
  #data2plot <- tibble::as_tibble(res_shrink) %>% 
  #  dplyr::mutate("-log10padj" = -log10(padj), 
  #                log2baseMean = log2(baseMean),
  #                IsSignificant = padj <= significance_threshold)
  
  count_up <- data_contrast %>%
    dplyr::filter(IsSignificant == TRUE & log2FoldChange >= foldchange_threshold) %>%
    dplyr::pull(Symbol)
  
  count_dn <- data_contrast %>%
    dplyr::filter(IsSignificant == TRUE & log2FoldChange <= -foldchange_threshold) %>%
    dplyr::pull(Symbol)
  
  plot_title <- paste0("DEA of ", contrast[1], ": ", contrast[2], " / ", contrast[3], "\n",
                       "Up = ", length(count_up), ". Down = ", length(count_dn))
  
  # Scatter plots:
  # log2 FC v -log10 adjusted p-value
  # log2 base mean v fold change
  for(axis_vals in list(c("log2FoldChange", "neglog10padj"), c("log2baseMean", "log2FoldChange"))){
    filename <- paste(name_out, axis_vals[[1]], axis_vals[[2]], sep="_")
    
    create_scatter_plot_my(data_contrast, x_label=axis_vals[[1]], y_label=axis_vals[[2]], fill_label="IsSignificant", 
                           filename=filename, labels=c(axis_vals[[1]], axis_vals[[2]], plot_title),
                           do_fit=NULL, size=0.5)
  }
  
  return(data_contrast)
}
