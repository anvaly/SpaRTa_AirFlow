# Author: Anna Lyubetskaya. Date: 23-01-22

# Re-annotate results of an existing model


## ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`


## PARAMETERS ----


# Threshold on effect estimate
sd_threshold <- 2


## PATHS ----


# Work directory
work_path <- "XXXX"
input_file <- "upd_lm_res_filt.txt"

work_list <- paste0(dir(work_path, full.names=TRUE), "/")
work_list <- work_list[which(!grepl(".txt", work_list))]


## ITERATE ----


invivo_gather_list <- list()
invitro_gather_list <- list()

for(p in work_list){
  
  # Name of the model
  model <- gsub("/", "", gsub(paste0(work_path, "|", input_file), "", p))
  
  # Model results
  input_model <- paste0(p, input_file)

  
  ## INGEST DATA ----
  
  
  # Read the model
  coef_gene_df <- unique(readr::read_delim(input_model, delim="\t")) %>%
    dplyr::filter(Estimate > 0)

  
  ## GATE ----
  
  
  gate_list <- list()
  
  ###########################
  
  # Gate 1
  gate_list[[1]] <- coef_gene_df %>%
    dplyr::filter((Manguso22_KPC_WT_v_Output_log2FC < 0 & Manguso22_KPC_WT_v_Output_nlog10pv > 1.3) | (Manguso22_Panc02_WT_v_Output_log2FC < 0 & Manguso22_Panc02_WT_v_Output_nlog10pv > 1.3)) %>%
    dplyr::mutate(Manguso22SigDn = ifelse(Manguso22_KPC_WT_v_Output_log2FC < Manguso22_Panc02_WT_v_Output_log2FC, Manguso22_KPC_WT_v_Output_log2FC, Manguso22_Panc02_WT_v_Output_log2FC)) %>%
    dplyr::select(Symbol, Symbol_mouse, Manguso22SigDn) %>%
    unique()

  # Gate 2
  gate_list[[2]] <- coef_gene_df %>%
    dplyr::filter((Manguso22_KPC_WT_v_Output_log2FC > 0 & Manguso22_KPC_WT_v_Output_nlog10pv > 1.3) | (Manguso22_Panc02_WT_v_Output_log2FC > 0 & Manguso22_Panc02_WT_v_Output_nlog10pv > 1.3)) %>%
    dplyr::mutate(Manguso22SigUp = ifelse(Manguso22_KPC_WT_v_Output_log2FC > Manguso22_Panc02_WT_v_Output_log2FC, Manguso22_KPC_WT_v_Output_log2FC, Manguso22_Panc02_WT_v_Output_log2FC)) %>%
    dplyr::select(Symbol, Symbol_mouse, Manguso22SigUp) %>%
    unique()

  # Gate 3
  gate_list[[3]] <- coef_gene_df %>%
    dplyr::filter((Manguso22_KPC_WT_v_Output_log2FC < -0.35) | (Manguso22_Panc02_WT_v_Output_log2FC < -0.35)) %>%
    dplyr::mutate(Manguso22MayBe = ifelse(Manguso22_KPC_WT_v_Output_log2FC < Manguso22_Panc02_WT_v_Output_log2FC, Manguso22_KPC_WT_v_Output_log2FC, Manguso22_Panc02_WT_v_Output_log2FC)) %>%
    dplyr::select(Symbol, Symbol_mouse, Manguso22MayBe)
  
  ###########################
  
  # Gate 4
  gate_list[[4]] <- coef_gene_df %>%
    dplyr::filter(CLexpr_PANC0203 >= 2 | CLexpr_PANC0403 >= 2) %>%
    dplyr::select(Symbol, Symbol_mouse, CLexpr_PANC0203, CLexpr_PANC0403) %>%
    unique()

  # Gate 5
  gate_list[[5]] <- coef_gene_df %>%
    dplyr::filter(KP2CLparentalExprNorm >= 2) %>%
    dplyr::mutate(KP2Expr = KP2CLparentalExprNorm) %>%
    dplyr::select(Symbol, Symbol_mouse, KP2Expr)
  
  ###########################
  
  # Gate 6
  gate_list[[6]] <- coef_gene_df %>%
    dplyr::filter(DepMapAllAdh_ExprCount >= 5) %>%
    dplyr::select(Symbol, Symbol_mouse, DepMapAllAdh_ExprCount)

  # # Gate 7
  # gate_list[[7]] <- coef_gene_df %>%
  #   dplyr::filter(KP2SyngBulk_Ductal_2 >= 0.1) %>%
  #   dplyr::mutate(KP2SyngCorr = KP2SyngBulk_Ductal_2) %>%
  #   dplyr::select(Symbol, Symbol_mouse, KP2SyngCorr)
  
  ###########################
  
  # Gate 8
  gate_list[[7]] <- coef_gene_df %>%
    dplyr::filter(DepMapLinAdh_SensCount >= 10 | DepMapLinAdh_DepCount >= 5) %>%
    dplyr::select(Symbol, Symbol_mouse, DepMapLinAdh_SensCount)
  
  # Gate 9
  gate_list[[8]] <- coef_gene_df %>%
    dplyr::filter(DepMapAllAdh_SensCount >= 100 | DepMapAllAdh_DepCount >= 50) %>%
    dplyr::select(Symbol, Symbol_mouse, DepMapAllAdh_SensCount)
  
  # Gate 10
  gate_list[[9]] <- coef_gene_df %>%
    dplyr::filter(DepMapAllOrg_SensCount >= 10 | DepMapAllOrg_DepCount >= 5) %>%
    dplyr::select(Symbol, Symbol_mouse, DepMapAllOrg_SensCount)
  
  # Gate 11
  gate_list[[10]] <- coef_gene_df %>%
    dplyr::filter(DepMapLinOrg_SensCount >= 5 | DepMapLinOrg_DepCount >= 2) %>%
    dplyr::select(Symbol, Symbol_mouse, DepMapLinOrg_SensCount)
  
  # Gate 12
  gate_list[[11]] <- coef_gene_df %>%
    dplyr::filter(CLDep_PANC0203 <= -0.5 | CLDep_PANC0403 <= -0.5) %>%
    dplyr::select(Symbol, Symbol_mouse, CLDep_PANC0203, CLDep_PANC0403)
  
  ###########################

  # Gate 13
  gate_list[[12]] <- coef_gene_df %>%
    dplyr::filter(IsTF == "Yes") %>%
    dplyr::select(Symbol, Symbol_mouse, IsTF) %>%
    unique()
  
  ###########################
  
  # Gate 14
  gate_list[[13]] <- coef_gene_df %>%
    dplyr::group_by(Symbol, Symbol_mouse) %>%
    dplyr::ungroup() %>%
    dplyr::select(Symbol, Symbol_mouse, InVitro_Screen_Included) %>%
    unique()

  # Gate 15
  gate_list[[14]] <- coef_gene_df %>%
    dplyr::group_by(Symbol, Symbol_mouse) %>%
    dplyr::ungroup() %>%
    dplyr::select(Symbol, Symbol_mouse, InVivo_Screen_Included) %>%
    unique()
  
  ###########################
  
  # Join all gates
  gate_df <- purrr::reduce(gate_list, dplyr::full_join) %>%
    dplyr::mutate(Model = model) %>%
    dplyr::mutate_if(is.numeric , tidyr::replace_na, replace = 0) %>%
    dplyr::mutate(FLAG = grepl("KRT|LAM|ITG|COL|MUC", Symbol)) %>%
    dplyr::arrange(FLAG, Manguso22SigDn, desc(Manguso22SigUp), Manguso22MayBe, IsTF, 
                   desc(CLexpr_PANC0203), desc(CLexpr_PANC0403), desc(KP2Expr))
  
  # Write to file
  filename <- paste0(p, "upd_lm_res_gated.txt")
  readr::write_delim(gate_df, filename, delim="\t")

  ###########################
  
  
  # Join all gates
  gate_invitro_df <- gate_df %>%
    dplyr::filter(FLAG == FALSE & (CLexpr_PANC0203 >= 2 | CLexpr_PANC0403 >= 2)) %>%
    dplyr::mutate(Priority1 = ifelse(Manguso22SigDn < 0 | Manguso22SigUp > 0 | Manguso22MayBe < 0, "P1", ""),
                  Priority2 = ifelse(!is.na(IsTF), "P2", ""),
                  # Priority3 = ifelse(KP2SyngCorr >= 0.5, "P3", ""),
                  Priority3 = ifelse(Priority1 == "" & Priority2 == "", CLexpr_PANC0203, ""),
                  Priority_InVitro = paste0(Priority1, Priority2, Priority3)) %>%
    dplyr::arrange(desc(Priority1), desc(Priority2), desc(Priority3), 
                   desc(max(CLexpr_PANC0203, CLexpr_PANC0403))) %>%
    dplyr::select(Symbol, Priority_InVitro, InVitro_Screen_Included, Model) %>%
    unique()
  
  # Write to file
  filename <- paste0(p, "upd_lm_res_gated_invitro.txt")
  readr::write_delim(gate_invitro_df, filename, delim="\t")

  
  # Join all gates
  gate_invivo_df <- gate_df %>%
    dplyr::filter(FLAG == FALSE & KP2Expr > 0 & !is.na(Symbol_mouse)) %>%
    dplyr::mutate(Priority1 = ifelse(Manguso22SigDn < 0 | Manguso22SigUp > 0 | Manguso22MayBe < 0, "P1", ""),
                  Priority2 = ifelse(!is.na(IsTF), "P2", ""),
                  # Priority3 = ifelse(KP2SyngCorr >= 0.5, "P3=0.5", ""),
                  # Priority4 = ifelse(KP2SyngCorr >= 0.1 & KP2SyngCorr < 0.5, "P3=0.25", ""),
                  Priority3 = ifelse(Priority1 == "" & Priority2 == "", KP2Expr, ""),
                  Priority_InVivo = paste0(Priority1, Priority2, Priority3)) %>%
    dplyr::arrange(desc(Priority1), desc(Priority2), desc(Priority3), desc(KP2Expr)) %>%
    dplyr::select(Symbol_mouse, Priority_InVivo, InVivo_Screen_Included, Model) %>%
    unique()
  
  # Write to file
  filename <- paste0(p, "upd_lm_res_gated_invivo.txt")
  readr::write_delim(gate_invivo_df, filename, delim="\t")

  ###########################
  
  
  # Write to file
  filename <- paste0(p, "upd_lm_res_gated.txt")
  readr::write_delim(gate_df %>%
                       dplyr::left_join(gate_invitro_df, by="Symbol") %>%
                       dplyr::left_join(gate_invivo_df, by="Symbol_mouse") %>%
                       dplyr::relocate(Priority_InVitro, Priority_InVivo), filename, delim="\t")
  
  ###########################
  
  invivo_gather_list[[p]] <- gate_invivo_df
  invitro_gather_list[[p]] <- gate_invitro_df
  
}


# Accumulate data across all models
# Filter for genes we actually already profiled
invivo_df <- dplyr::bind_rows(invivo_gather_list) %>%
  # dplyr::filter(InVivo_Screen_Included == "Y") %>%
  dplyr::group_by(Symbol_mouse, Priority_InVivo, InVivo_Screen_Included) %>%
  dplyr::summarise(Models = paste(Model, collapse="; "))

invitro_df <- dplyr::bind_rows(invitro_gather_list) %>%
  # dplyr::filter(InVitro_Screen_Included == "Y") %>%
  dplyr::group_by(Symbol, Priority_InVitro, InVitro_Screen_Included) %>%
  dplyr::summarise(Models = paste(Model, collapse="; "))

# Write to file
filename <- paste0(work_path, "lm_res_gated_invivo_", sd_threshold, ".txt")
readr::write_delim(invivo_df, filename, delim="\t")

filename <- paste0(work_path, "lm_res_gated_invitro_", sd_threshold, ".txt")
readr::write_delim(invitro_df, filename, delim="\t")
