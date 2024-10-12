# Author: Anna Lyubetskaya. Date: 23-01-18

# Re-annotate results of an existing model


## ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_ggplot.R")


## PARAMETERS ----


# Add the rest of the evaluated variables from the model
add_confounders <- FALSE
# Keep the following X first columns of a file
keep_n <- 8

# Threshold on effect estimate
sd_threshold <- 2


## PATHS ----


# Work directory
work_path <- "XXXX"

# Find input files
work_list <- paste0(dir(work_path, full.names=TRUE), "/")
work_list <- work_list[which(!grepl(".txt", work_list))]

# Various annotation files
annotations_human_path <- c("XXXX",
                            "XXXX",
                            "XXXX",
                            "XXXX",
                            "XXXX",
                            "XXXX",
                            "XXXX",
                            "XXXX",
                            # "XXXX",
                            "XXXX"
                            )

annotations_mouse_path <- c("XXXX",
                            "XXXX",
                            "XXXX",
                            "XXXX",
                            "XXXX",
                            "XXXX",
                            "XXXX",
                            "XXXX",
                            "XXXX",
                            "XXXX",
                            "XXXX",
                            "XXXX"
                            )


## ITERATE ----


for(p in work_list){
  
  # Model results
  input_model <- paste0(p, "gene_lm_result_scatter.txt")
  input_confounders <- paste0(p, "gene_lm_result_confounders.txt")
  
  # Name of the model
  model <- gsub(".+/", "", gsub("/$", "", p))
  
  
  ## INGEST DATA ----
  
  
  # Read the model
  coef_gene_df <- unique(readr::read_delim(input_model, delim="\t")[, 1:keep_n])
  
  
  ## ANNOTATE PREVIOUS MODEL RESULTS ----
  
  
  # Pull additional external annotations to build up evidence
  for(f in annotations_human_path){
    print(f)
    
    coef_gene_df <- coef_gene_df %>%
      dplyr::left_join(readr::read_delim(f) %>%
                         dplyr::mutate(Symbol = toupper(Symbol)), by="Symbol")
  }
  
  # Pull additional external annotations to build up evidence
  for(f in annotations_mouse_path){
    coef_gene_df <- coef_gene_df %>%
      dplyr::left_join(readr::read_delim(f), by="Symbol_mouse")
  }
  
  # Add the information about quantitative confounders of the model
  if(add_confounders == TRUE){
    coef_conf_df <- readr::read_delim(input_confounders, delim="\t")
    
    coef_gene_ann_df <- coef_gene_df %>%
      dplyr::left_join(coef_conf_df, by="Symbol")
  }

  
  # Establish the estimate threshold as mean +/- 2SD
  estimate_threshold <- mean(coef_gene_df$Estimate) + sd(coef_gene_df$Estimate) * sd_threshold
  score_threshold <- round(mean(coef_gene_df$Score))
  
  # Round numbers and crop
  coef_gene_ann_df <- coef_gene_df %>%
    dplyr::mutate(IsSignificant = ifelse(abs(Estimate) >= estimate_threshold & Score >= score_threshold,
                                         TRUE, FALSE)) %>%
    dplyr::mutate_if(is.numeric, round, 4) %>%
    dplyr::arrange(desc(IsSignificant), desc(Estimate), desc(Score))
  
  # Write to file
  filename <- paste0(p, "upd_lm_res_ann.txt")
  readr::write_delim(coef_gene_ann_df, filename, delim="\t")
  
  
  # Visualize the result of linear modeling of gene expression ~ distance to center
  filename <- paste0(p, "upd_lm_res_scatter")
  create_scatter_plot_my(coef_gene_df, x_label="Estimate", y_label="Score", 
                         fill_label="IsSignificant", shape=21, size=1, dot_labels=NULL, 
                         filename=filename, labels=NULL, do_fit=NULL, stroke=0)
  

  # Select genes of interest for targeted follow up
  coef_gene_ann_filt_df <- coef_gene_ann_df %>%
    dplyr::filter(IsSignificant == TRUE) %>%
    dplyr::arrange(desc(Estimate), desc(Score)) %>%
    dplyr::relocate(InVivo_Screen_Included)
  
  
  # Write to file
  filename <- paste0(p, "upd_lm_res_filt.txt")
  readr::write_delim(coef_gene_ann_filt_df, filename, delim="\t")
  
  
  # Select columns of interest for targeted follow up
  coef_gene_ann_filt_df <- coef_gene_ann_filt_df %>%
    dplyr::select(-IsSignificant, -Labels, 
                  -GrowthPatternLinAdh, -GrowthPatternLinOrg, -GrowthPatternAllAdh, -GrowthPatternAllOrg,
                  -DepMapAllAdh_ExprMean, -DepMapAllAdh_SensCount, -DepMapAllAdh_DepCount, 
                  -DepMapAllOrg_ExprMean, -DepMapAllOrg_SensCount, -DepMapAllOrg_DepCount,
                  -DepMapAllAdh_ExprSD, -DepMapAllAdh_ExprCount, -DepMapAllAdh_MutCount,
                  -DepMapAllOrg_ExprSD, -DepMapAllOrg_ExprCount, -DepMapAllOrg_MutCount,
                  -DepMapLinAdh_ExprSD, -DepMapLinAdh_ExprCount, -DepMapLinAdh_MutCount,
                  -DepMapLinOrg_ExprSD, -DepMapLinOrg_ExprCount, -DepMapLinOrg_MutCount,
                  -Ensembl_gene_ID_human, -Ensembl_gene_ID_mouse, -InVitro_Screen_Included,
                  -Manguso22_KPC_ICB_v_WT_log2FC, -Manguso22_KPC_ICB_v_WT_nlog10pv,
                  -Manguso22_Panc02_ICB_v_WT_log2FC, -Manguso22_Panc02_ICB_v_WT_nlog10pv,
                  -Manguso22_KPC_Output_v_Input_log2FC, -Manguso22_KPC_Output_v_Input_nlog10pv,
                  -Manguso22_Panc02_Output_v_Input_log2FC, -Manguso22_Panc02_Output_v_Input_nlog10pv,
                  -Manguso22_KPC_WT_v_NGS_log2FC, -Manguso22_KPC_WT_v_NGS_nlog10pv,
                  -Manguso22_Panc02_WT_v_NGS_log2FC, -Manguso22_Panc02_WT_v_NGS_nlog10pv,
                  -Manguso22_KPC_WT_v_Input_log2FC, -Manguso22_KPC_WT_v_Input_nlog10pv,
                  -Manguso22_Panc02_WT_v_Input_log2FC, -Manguso22_Panc02_WT_v_Input_nlog10pv,
                  -CLexpr_PANC0203, -CLexpr_PANC0403, -CLDep_PANC0203, -CLDep_PANC0403,
                  -DepMapLinAdh_DepCount, -DepMapLinOrg_DepCount
    ) %>%
    dplyr::rename(`InScreen?` = InVivo_Screen_Included) %>%
    unique() %>%
    dplyr::group_by(Symbol) %>%
    dplyr::mutate(H2M_mapping = dplyr::n_distinct(Symbol_mouse),
                  Model = gsub("PDAC108_path14_harmonyepi_rpca_sct_|sig.|BMS.CL.|PDAC.|collisson.|moffitt.|", "", model),
                  Direction = ifelse(Estimate > 0, "UP", "DN")) %>%
    dplyr::select(Model, `InScreen?`, Symbol, Symbol_mouse, H2M_mapping, IsTF, 
                  Estimate, Score, SpotMean, SpotCount,	KP2CLparentalExprNorm, Sig_Name,
                  DepMapLinAdh_ExprMean, DepMapLinAdh_SensCount, DepMapLinOrg_ExprMean, DepMapLinOrg_SensCount,
                  Manguso22_KPC_WT_v_Output_log2FC, Manguso22_KPC_WT_v_Output_nlog10pv, Manguso22_Panc02_WT_v_Output_log2FC, Manguso22_Panc02_WT_v_Output_nlog10pv, Direction) %>%
    dplyr::arrange(desc(Direction), desc(`InScreen?`), desc(IsTF), desc(Estimate), desc(Score))
  
  # Write to file
  filename <- paste0(p, "upd_lm_res_filt_part.txt")
  readr::write_delim(coef_gene_ann_filt_df, filename, delim="\t")
  

}
