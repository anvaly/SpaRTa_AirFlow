
#####


run_min = FALSE

source("code/utils/utils_lin_regress.R")


# Wide tibble of all measurements
seurat_values_wide_df <- tibble::as_tibble(t(as.matrix(Seurat::GetAssayData(data_seurat, assay="SCT"))), rownames="id") %>%
  dplyr::rename(Coordinate = id)

# Find genes that have more than one value across samples
genes_changing <- seurat_values_wide_df %>%
  df_wide2long_my(key="Symbol", val="Count") %>%
  dplyr::group_by(Symbol) %>%
  dplyr::summarise(Count_unique = dplyr::n_distinct(Count),
                   Count_median = median(Count)) %>%
  dplyr::filter(Count_unique > 2 &
                  !grepl("RP11-", Symbol)) %>%
  dplyr::arrange(desc(Count_median), desc(Count_unique)) %>%
  dplyr::pull(Symbol)


########################################


output_folders <- create_output_subfolders_my(output_folder, c("group_dea_elnet_kenzie"))

# Find the ratio between two signatures of interest
sig_ratio_df <- sig_select_df %>%
  dplyr::mutate(Ratio_spot = log2(!!sym(sig1) / !!sym(sig2)),
                Ratio_area = log2(Sig_area1 / Sig_area2),
                LibSize_log2 = log2(nCount_Spatial)) %>%
  tidyr::drop_na() %>%
  #dplyr::filter(nCount_Spatial <= mean(nCount_Spatial) + sd(nCount_Spatial),
  #              Category_sig1 == "Hi") %>%
  dplyr::select(Coordinate, LibSize_log2, Ratio_spot, Ratio_area)


elnet_fit_per_gene <- function(gene){
  # Cleanup data
  data_model_df <- sig_ratio_df %>%
    dplyr::inner_join(seurat_values_wide_df[c("Coordinate", gene)], by="Coordinate") %>%
    dplyr::rename("outcome" = !!rlang::sym(gene)) %>%
    tidyr::drop_na()
  
  # Format input data for modeling
  data_model_df <- data_model_df %>%
    tibble::column_to_rownames("Coordinate")
  
  # Try to fit various models
  model <- run_models_my(data_model_df, c("elnet"), filename, formula=formula)
  
  # Analyze performance of various linear regressions
  filename <- paste0(output_folders$group_dea, "elnet")
  coef_df <- analyze_model_my(model, c("elnet"))
  
  # Simplify column names
  colnames(coef_df) <- gsub("elnet:.+", gene, colnames(coef_df))
  
  return(coef_df)
}


formula <- "outcome ~ Ratio_spot*Ratio_area + LibSize_log2"

cl <- parallel::makeCluster(parallel::detectCores() - 1)

doParallel::registerDoParallel(cl)

coef_list <- foreach::foreach(gene = genes_changing, .combine="c") %dopar% {
  elnet_fit_per_gene(gene)
};

parallel::stopCluster(cl)

result <- coef_list %>% 
  purrr::reduce(dplyr::full_join, by="variable")

result_df <- result %>%
  df_wide2long_my(key="Symbol", val="Score")

result_stat_df <- result %>%
  df_wide2long_my(key="Symbol", val="Score") %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(mean = mean(Score),
                   sd = sd(Score),
                   min = min(Score),
                   max = max(Score))



########################################


output_folders <- create_output_subfolders_my(output_folder, c("group_dea_elnet2"))

# Find the ratio between two signatures of interest
sig_ratio_df <- sig_select_df %>%
  dplyr::mutate(Ratio = round(log2(!!sym(sig1) / !!sym(sig2)), 2)) %>%
  tidyr::drop_na() %>%
  dplyr::filter(nCount_Spatial <= mean(nCount_Spatial) + sd(nCount_Spatial),
                Category_sig1 == "Hi") %>%
  dplyr::select(Coordinate, nCount_Spatial, Ratio)

# Cleanup data
data_model_df <- sig_ratio_df %>%
  dplyr::inner_join(seurat_values_wide_df[c("Coordinate", genes_changing)], by="Coordinate") %>%
  dplyr::rename("outcome" = "Ratio") %>%
  tidyr::drop_na()

# Format input data for modeling
data_model_df <- data_model_df %>%
  tibble::column_to_rownames("Coordinate")

# Try to fit various models
model <- run_models_my(data_model_df, c("elnet"), filename)

# Analyze performance of various linear regressions
filename <- paste0(output_folders$group_dea, "elnet")
coef_df <- analyze_model_my(model, c("elnet"), filename)

# Simplify column names
colnames(coef_df) <- gsub("elnet:.+", "coef", colnames(coef_df))

# Select non-zero coefficients
coef_filt_df <- coef_df %>% 
  dplyr::filter(coef != 0)

# Find intersection between non-zero coefficients and our signatures
intersect_list <- sapply(signature_list, function(x) intersect(coef_filt_df$variable, x))
intersect_df <- tibble::tibble(Signature_name = names(intersect_list),
                               Intersect_len = unname(sapply(intersect_list, function(x) length(x))),
                               Intersect = unname(sapply(intersect_list, function(x) paste0(x, collapse=";"))))

# Combine signature length information with the signature genes found in modeling effort
sig_lengths_df <- sig_intersect_df %>% 
  dplyr::select(Signature_name, Sig_length) %>% 
  unique %>%
  dplyr::inner_join(intersect_df, by="Signature_name") %>%
  dplyr::mutate(Ratio = round(Intersect_len / Sig_length * 100)) %>%
  dplyr::arrange(desc(Ratio)) %>%
  dplyr::select(Signature_name, Sig_length, Intersect_len, Ratio, Intersect)

readr::write_delim(sig_lengths_df, path=paste0(output_folders$group_dea, "coefficients_by_signature.txt"), delim = "\t")

coef_filt_plot_df <- coef_filt_df %>% 
  dplyr::mutate(InSignature = variable %in% unique(unlist(unname(signature_list))))

# Compare coefficients proposed by various linear regression fits
p <- create_bar_plot_my(coef_filt_plot_df %>%
                          dplyr::top_n(20, abs(coef)), 
                        x_label="variable", y_label="coef", 
                        fill_label="InSignature", position="dodge", 
                        filename=paste0(output_folders$group_dea, "coefficients_filt"), 
                        labels=c("Gene Symbol", "ElNet coefficient", ""), 
                        reorder_x=TRUE)

readr::write_delim(coef_filt_plot_df, path=paste0(output_folders$group_dea, "coefficients_final.txt"), delim = "\t")


##############


output_folders <- create_output_subfolders_my(output_folder, c("group_dea_elnet1"))

# Find the ratio between two signatures of interest
sig_ratio_df <- sig_select_df %>%
  dplyr::mutate(Ratio = round(log2(!!sym(sig1) / !!sym(sig2)), 2)) %>%
  tidyr::drop_na() %>%
  dplyr::select(Coordinate, nCount_Spatial, Ratio) %>%
  dplyr::filter(nCount_Spatial <= mean(nCount_Spatial) + sd(nCount_Spatial))

# Cleanup data
data_model_df <- sig_ratio_df %>%
  dplyr::inner_join(seurat_values_wide_df[c("Coordinate", genes_changing)], by="Coordinate") %>%
  dplyr::rename("outcome" = "Ratio") %>%
  tidyr::drop_na()

# Format input data for modeling
data_model_df <- data_model_df %>%
  tibble::column_to_rownames("Coordinate")

# Try to fit various models
model <- run_models_my(data_model_df, c("elnet"), filename)

# Analyze performance of various linear regressions
filename <- paste0(output_folders$group_dea, "elnet")
coef_df <- analyze_model_my(model, c("elnet"), filename)

# Simplify column names
colnames(coef_df) <- gsub("elnet:.+", "coef", colnames(coef_df))

# Select non-zero coefficients
coef_filt_df <- coef_df %>% 
  dplyr::filter(coef != 0)

# Find intersection between non-zero coefficients and our signatures
intersect_list <- sapply(signature_list, function(x) intersect(coef_filt_df$variable, x))
intersect_df <- tibble::tibble(Signature_name = names(intersect_list),
                               Intersect_len = unname(sapply(intersect_list, function(x) length(x))),
                               Intersect = unname(sapply(intersect_list, function(x) paste0(x, collapse=";"))))

# Combine signature length information with the signature genes found in modeling effort
sig_lengths_df <- sig_intersect_df %>% 
  dplyr::select(Signature_name, Sig_length) %>% 
  unique %>%
  dplyr::inner_join(intersect_df, by="Signature_name") %>%
  dplyr::mutate(Ratio = round(Intersect_len / Sig_length * 100)) %>%
  dplyr::arrange(desc(Ratio)) %>%
  dplyr::select(Signature_name, Sig_length, Intersect_len, Ratio, Intersect)

readr::write_delim(sig_lengths_df, path=paste0(output_folders$group_dea, "coefficients_by_signature.txt"), delim = "\t")

coef_filt_plot_df <- coef_filt_df %>% 
  dplyr::mutate(InSignature = variable %in% unique(unlist(unname(signature_list))))

# Compare coefficients proposed by various linear regression fits
p <- create_bar_plot_my(coef_filt_plot_df %>%
                          dplyr::top_n(20, abs(coef)), 
                        x_label="variable", y_label="coef", 
                        fill_label="InSignature", position="dodge", 
                        filename=paste0(output_folders$group_dea, "coefficients_filt"), 
                        labels=c("Gene Symbol", "ElNet coefficient", ""), 
                        reorder_x=TRUE)

readr::write_delim(coef_filt_plot_df, path=paste0(output_folders$group_dea, "coefficients_final.txt"), delim = "\t")