# Author: Anna Lyubetskaya. Date: 20-01-23
# This code is based on Gene.Rmd from ModularClinicalAnalysis by Silpa Suthram
# This code uses trailESET_functions.R, edgeR_Normalize.R, and addGeneSetScores.R from ModularClinicalAnalysis.

# Clean tryCatch example: https://stackoverflow.com/questions/12193779/how-to-write-trycatch-in-r
# On error handling and exceptions in R: http://adv-r.had.co.nz/beyond-exception-handling.html


source_path <- "C:/Users/lyubetsa/Documents/Git/Modular_Clinical_Reports/code/R/"

source("code/utils/utils_in_out.R")

source(paste0(source_path, "trialESET_functions.R"))
#source(paste0(source_path, "edgeR_Normalize.R"))
#source(paste0(source_path, "addGeneSetScores.R"))
#source(paste0(source_path, "crossTrial_OR_Table.R"))
#source(paste0(source_path, "ResponseAssoc_BoxPlot.R"))


environment_rXpress_my <- function(){
  ## Check necessary package versions
  ## If install.packages doesn't work, try grabbing them from the BRAN website: http://bran.pri.bms.com/
  
  # You need SiteMinderBMS 1.6.6 or newer in order to load datasets from XPRESS
  if(packageVersion("SiteMinderBMS") < '1.6.6'){
    install.packages("SiteMinderBMS")
  }
  
  # You need rXpress 1.4.21 or newer so that we can use the function rXpress::getXpressEset
  if(packageVersion("rXpress") < '1.4.21'){
    install.packages("rXpress")
  }
}


# Setup environment if working with these functions
environment_rXpress_my()


meta_ifelse <- function(matrix, value){
  ## A utility function to conditionally look up a value
  
  if(value %in% colnames(matrix)){
    return(matrix[[value]])
  } else{
    return(NA)
  }
}


getTrialEset_tryCatch_my <- function(trial_id, xpressIDMapTable, clinVarMapTable) {
  ## Wrap getTrialEset function in a tryCatch wrapper
  
  counts_data <- tryCatch(
    {
      message(paste0(trial_id, " OBTAINED", "\n"))
      getTrialEset(trialID = trial_id, xpressMapTable = xpressIDMapTable, clinMapTable = clinVarMapTable)
    },
    error=function(e) {
      message(paste0(trial_id, " FAILED\n", e, "\n"))
      return(NULL)
    },
    warning=function(w) {
      message(paste0(trial_id, " WARNING\n", w, "\n"))
      return(NULL)
    },
    finally={
      # message("next")
    }
  )
  
  return(counts_data)
}


subset_trial_eset_my <- function(normEset, timepoint_value, arm_value){
  ## If given a time point and treatment arm, subset a full expression matrix to only those patients
  
  # Trial timepoints
  timepoint_list <- Biobase::pData(normEset)$Timepoint
  # Filter for samples from a particular timepoint
  normEset <- normEset[, which(timepoint_list == timepoint_value)]
  
  # Trial treatments
  treatment_list <- Biobase::pData(normEset)$TreatmentStandard
  
  cat("\nTrial=", trial_id, "; timepoints=", unique(timepoint_list), "; treatments=", unique(treatment_list), "\n")
  # Check if samples for a treatment arm of interest exist
  if(arm_value %in% unique(treatment_list)){
    # Filter for a treatment arm
    normEset <- normEset[,which(treatment_list == arm_value)]
  }
  
  return(normEset)
}


load_trial_eset_from_list_my <- function(trial_list, xpressIDMapTable, clinVarMapTable, arm_value="Nivo", timepoint_value="Pre-Treatment", 
                                         geneset_list=NULL, do_qc=TRUE, do_norm=TRUE){
  ## Create an list of eSet objects for a list of trials and pre-process the data.
  ## Choose only specific timepoint and treatment arm.
  
  cat("Trials loading: ", trial_list, "\n")
  
  esetList = list()
  
  for(trial_id in trial_list){
    # Read trial data using getTrialEset wrapped in a tryCatch
    eset <- getTrialEset_tryCatch_my(trial_id, xpressIDMapTable, clinVarMapTable)
    # print(eset)
    
    # In case of failure to load or non-eSet return (e.g. EdgeSeq option returns a list of two eSets)
    if(!is.null(eset)){
      # QC trial data
      if(do_qc == TRUE){
        eset <- qcTrialEset(eset, EstLibrarySizeOver = 5000000)
      }
      
      # Normalize counts, output CPMs
      if(do_norm == TRUE){
        eset <- edgeR_Normalize(eset) 
      }
      
      # Subset the expression set
      if(!is.null(timepoint_value)){
        eset <- subset_trial_eset_my(eset, timepoint_value, arm_value)
      }
      
      # Add geneset scores
      eset <- addGeneSetScores(eset, geneset_list)
      
      # Add final eSet to a list
      esetList[[trial_id]] <- eset
    }
  }
  
  cat("Trials loaded: ", names(esetList), "\n")
  
  return(esetList)
}


odds_ratio_table_box_my <- function(eset_list, gene, confounders=c("SPECTYPE"), path_out="box_association.png"){
  ## Calculate odds ratio table and plot boxplots of association of response with biomarker expression
  
  # Calculate odds ratio table
  or_table <- crossTrial_OR_Table(eset_list, marker = gene, confounders_phen = confounders)
  print(or_table)
  
  # Plot boxplots association response with expression of a biomarker
  p <- ResponseAssoc_BoxPlot(eset_list, marker = gene, OR_Table = or_table)
  png(path_out, width = 1000, height = 1000)
  plot(p)
  dev.off()
}


read_mutation_data_from_kenzie_my <- function(trial_list=NULL, qc=TRUE, essential_fields=TRUE){
  ## Load all mutation data Kenzie processed in one long tibble
  ## Columns in input files: USUBJID, Gene, AAChange, Frequency, SIFT, PolyPhen, IMPACT, CLIN_SIG
  ## qc flag checks that neither AAChange nor Frequency fields are NA which leaves 78% of mutations
  
  # A default list of trials based on v1.6
  if(is.null(trial_list)){
    trial_list <- c("CA209-017", "CA209-025", "CA209-026", "CA209-032", "CA209-038", "CA209-040", "CA209-057",
                    "CA209-064", "CA209-066", "CA209-067", "CA209-141", "CA209-142", "CA209-238", "CA209-275")
  }
  
  # Location of Kenzie's analysis of clinical mutation data
  input_mutations <- "XXXX"
  
  # Upload all mutation information from trials of interest in a tibble
  if(is.null(trial_list)){
    mutations_df <- read_dir2file_my(input_mutations, in_regex=".txt$")
  } else{
    mutations_df <- read_dir2file_my(input_mutations, in_regex=paste(trial_list, collapse="|"))
  }
  
  if(qc == TRUE){
    mutations_df <- mutations_df %>%
      dplyr::filter(!is.na(AAChange) & !is.na(Frequency))
  }
  
  if(essential_fields == TRUE){
    mutations_df <- mutations_df %>%
      dplyr::select(USUBJID, Gene, AAChange)
  }
  
  return(mutations_df %>%
           unique) 
}


collect_meta_from_stash_my <- function(clinVarMapTable=NULL, trial_list=NULL, patient_id_col="USUBJID"){
  ## Collect meta data from stash
  ## Patient ID column: patient_id_col = "USUBJID"
  
  # Clinical data global path
  input_clinical <- "XXXX"
  
  # If clinical data table is not loaded, load it now
  if(is.null(clinVarMapTable)){
    clinVarMapTable <- read.csv("XXXX", stringsAsFactors = F)
  }
  
  # If a list of trial IDs is not provided, grab from the clinical data table
  if(is.null(trial_list)){
    trial_list <- clinVarMapTable %>%
      dplyr::pull("Trial_ID") %>%
      unique
  }
  
  # For each trial load important clinical data
  # Standardize the clinical data using the meta data file as a guide
  clinical_data <- list()
  for(trial in trial_list){
    # Find important standard variables for this trial
    variables_df <- clinVarMapTable %>%
      dplyr::filter(Trial_ID == trial)
    
    # Extract indication for this trial
    indication <- variables_df %>%
      dplyr::filter(variable == "TumorType_Organ") %>%
      dplyr::pull(value)
    
    # Find the location of clinical data for this trial
    clinical_file <- variables_df %>%
      dplyr::filter(variable == "clinFilePath") %>%
      dplyr::pull(value)
    
    # Load all clinical data for this trial
    # Some clinical data is tab-delimited and some comma-delimitted
    if(file.exists(clinical_file)){
      clinical_data[[trial]] <- read_file2df_my(clinical_file, delim="\t")
      if(length(colnames(clinical_data[[trial]])) == 1){
        clinical_data[[trial]] <- read_file2df_my(clinical_file, delim=",")
      }
      
      # Find columns suggested by the meta data guide among loaded data
      variables_select_df <- variables_df %>%
        dplyr::filter(value %in% colnames(clinical_data[[trial]]))
      
      # Named list of important clinical meta data variables
      variables_list <- setNames(variables_select_df$value, variables_select_df$variable)
      
      # Cleanup the clinical data by renaming columns as suggested by the meta data guide
      clinical_data[[trial]] <- clinical_data[[trial]] %>%
        dplyr::select(c(patient_id_col, unname(variables_list))) %>%
        dplyr::rename(!!! variables_list)
      
      # If indication is not mentioned explicitly, add a value
      if(!"TumorType_SampleLevel" %in% clinical_data[[trial]]){
        clinical_data[[trial]][["TumorType_SampleLevel"]] <- indication
      }
      
      # Remove patients that failed screen or had no indication
      clinical_data[[trial]] <- clinical_data[[trial]] %>%
        dplyr::filter(!is.na(TumorType_SampleLevel) & TreatmentArm != "Screen Failure")
    }
  }
  
  # All clinical data in one tibble
  clinical_df <- plyr::rbind.fill(clinical_data)
  
  cat("Disease as reported in TumorType_SampleLevel:", 
      paste(clinical_df$TumorType_SampleLevel %>% unique(), sep=", "), "\n")
  
  return(clinical_df)
}


collect_maf_files_from_stash_my <- function(){
  ## Find all MAF files listed in cBioPortal folders on stash
  
  # Location of our mutation files.
  global_path1 <- "XXXX"
  global_path2 <- "XXXX"
  
  # Collect all files of interest.
  file_list <- c(system(paste0("ls ", global_path1), intern=TRUE),
                 system(paste0("ls ", global_path2), intern=TRUE))
  
  
  return(file_list)
}


clinical_eset_to_tibble_my <- function(eset_list, gene_interest_list, geneset_list=NULL){
  ## From a list of eSet, gather expression measurements and relevant patient meta information into a one long tibble
  

  # Find Ensembl gene ID for the gene symbol of interest using the first of clinical trials loaded
  symbol2ensembl_df <- Biobase::fData(eset_list[[1]]) %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(Symbol %in% gene_interest_list) %>% 
    dplyr::select(vars, Symbol)
  
  # Gather all relevant expression and patient data into the following tibbles for the genes of interest
  expression_interest_df <- tibble::as_tibble(list("EnsemblID"=NA, "Patient_ID"=NA, "log2CPM_norm"=NA, "Trial_ID"=NA))
  
  # Sample meta data columns
  patients_rna_df <- tibble::as_tibble(list("USUBJID"=NA, "Patient_ID"=NA, "TumorType_Organ"=NA, "TreatmentStandard"=NA, 
                                            "Timepoint"=NA, "BOR"=NA, "OS"=NA, "PFS"=NA, "OS.CNSR"=NA, "PFS.CNSR"=NA,
                                            "TMB"=NA, "Platform"=NA))

  # Initiate columns for signatures
  if(!is.null(geneset_list)){
    for(sig in names(geneset_list)){
      patients_rna_df[[sig]] <- NA
    }
  }

  # Remove annotation ID from patient ID
  string_trim <- "GRCh\\d+ERCC-ensembl\\d+-genes-"
  
  for(trial_id in names(eset_list)){
    # Grab an expression matrix for a given trial
    data_matrix <- Biobase::exprs(eset_list[[trial_id]])
    
    # Grab the expression vector from a given trial for a gene of interest
    expr_vector <- data_matrix[rownames(data_matrix) %in% symbol2ensembl_df$vars | 
                                 rownames(data_matrix) %in% symbol2ensembl_df$Symbol, ]
    
    if(length(expr_vector) > 1){
      
      ## Collect expression data across all trials for genes of interest  
      
      # Turn expression vector into a tibble
      loc_df <- tibble::as_tibble(expr_vector, rownames="id") %>%
        dplyr::rename("EnsemblID" = "id") %>%
        df_wide2long_my(key="Patient_ID", val="log2CPM_norm")
      
      # Explicitly add trial name
      loc_df[["Trial_ID"]] <- trial_id
      
      # Collect all expression data in one tibble
      expression_interest_df <- rbind(expression_interest_df, loc_df)
      
      ## Collect patient IDs
      
      # Grab patient meta data for a given trial
      patient_matrix <- Biobase::pData(eset_list[[trial_id]])
      
      # Grab the patient vector from a given trial
      pat_df <- tibble::as_tibble(list("USUBJID" = patient_matrix$USUBJID, 
                                       "Patient_ID" = rownames(patient_matrix),
                                       "TumorType_Organ" = patient_matrix$TumorType_Organ,
                                       "TreatmentStandard" = patient_matrix$TreatmentStandard,
                                       "Timepoint" = patient_matrix$Timepoint,
                                       "BOR" = meta_ifelse(patient_matrix, "BOR"),
                                       "OS" = meta_ifelse(patient_matrix, "OS"), 
                                       "PFS" = meta_ifelse(patient_matrix, "PFS"),
                                       "OS.CNSR" = meta_ifelse(patient_matrix, "OS.CNSR"), 
                                       "PFS.CNSR" = meta_ifelse(patient_matrix, "PFS.CNSR"),
                                       "TMB" = meta_ifelse(patient_matrix, "TMB"),
                                       "Platform" = Biobase::fData(eset_list[[trial_id]])$Platform[1]))
      
      # Add signature information
      if(!is.null(geneset_list)){
        for(sig in names(geneset_list)){
          pat_df[[sig]] <- patient_matrix[[sig]]
        }
      }
      
      # Collect all patient IDs in one tibble
      patients_rna_df <- rbind(patients_rna_df, pat_df)
    }
  }
  
  return(expression_interest_df %>%
           dplyr::inner_join(patients_rna_df, by="Patient_ID") %>%
           dplyr::inner_join(symbol2ensembl_df, by=c("EnsemblID" = "vars")) %>%
           dplyr::mutate(Patient_ID = gsub(string_trim, "", Patient_ID)) %>%
           dplyr::ungroup() %>%
           unique)
}
