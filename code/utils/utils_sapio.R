# Author: Anna Lyubetskaya. Date: 20-02-03
# Functions to work with our LIMS system - Sapio
# Managed through BMS internal SapioR package https://biogit.pri.bms.com/tilfordc/P02131_Sapio_R_API


environment_bran_SapioR_my <- function(){
  ## Install SapioR from BRAN
  ## https://biogit.pri.bms.com/tilfordc/P02131_Sapio_R_API
  
  if(!"SapioR" %in% rownames(installed.packages())){
    source('http://bran.pri.bms.com/resources/configureRepo.R')
    install.packages('SapioR')
  }
}


extract_meta_data_sapio_my <- function(id_list, output_file=NULL){
  ## From a presumed list of sample IDs, discover their patients ID and gather specified meta fields
  
  environment_bran_SapioR_my()
  
  # Find Sapio ID types for a list of IDs
  # Important: This doesn't always work for some reason
  id_types <- SapioR::classifySapioIds(id_list)
  # Explicitly tag IDs that were not identified
  id_types[is.na(as.vector(id_types))] <- "NotFound"
  
  # Extract input samples Sapio IDs
  # sample_list <- names(id_types[as.vector(id_types) == "Sample"])
  sample_list <- names(id_types)
  
  # Select the following sample meta data fields
  sample_fields <- c("Sample.SampleId", "Sample.LongitudinalTimepoint", "Sample.StudyId")  # "Sample.SampleName"
  # Select the following patient meta data fields
  patient_fields <- c("Patient.PatientId", "CambridgePatient.LIMSRegistrationType", "CambridgePatient.Sex",
                      "CambridgePatient.CancerSite", "CambridgePatient.CancerType", "CambridgePatient.StageAtDX")
  
  if(length(sample_list) > 0){
    # Annotate input samples using their Sapio IDs
    sample_meta <- SapioR::sapioSampleMetadata(sample_list, trimNA=TRUE, trimInvariant=TRUE) %>%
      dplyr::select(sample_fields)
    
    # Find all studies corresponding to input samples
    study_list <- sample_meta$Sample.StudyId %>%
      paste(collapse=",") %>%
      strsplit(",") %>%
      unlist() %>%
      unique
    
    # Find all patients included in studies corresponding to input samples
    studies2patients <- SapioR::sapioPatientsForStudy(study_list)
    
    # Find patients corresponding to input samples
    patients2samples <- SapioR::sapioSamplesForPatients(studies2patients$Patient.PatientId) %>% 
      dplyr::filter(Sample.SampleId %in% sample_list)
    
    # List patients corresponding to input samples
    patient_list <- patients2samples$Patient.PatientId %>% 
      unique
    
    # Annotate patients corresponding to input samples
    patient_meta <- SapioR::sapioPatientMetadata(patient_list, trimNA=TRUE, trimInvariant=TRUE) %>%
      dplyr::select(patient_fields)
    
    # Join all information in a wide tibble
    meta_data_df <- tibble::as_tibble(list("IdType" = as.vector(id_types), "Sample.SampleId" = names(id_types))) %>% 
      dplyr::left_join(sample_meta, by="Sample.SampleId") %>%
      dplyr::left_join(tibble::as_tibble(patients2samples), by="Sample.SampleId") %>%
      dplyr::left_join(patient_meta, by="Patient.PatientId") %>%
      dplyr::arrange(IdType, CambridgePatient.LIMSRegistrationType, Patient.PatientId, Sample.LongitudinalTimepoint)
    
    # Simplify column names
    colnames(meta_data_df) <- gsub("^[^\\.]+\\.", "", colnames(meta_data_df))
    
    # Write resulting meta data table to file for reference
    if(!is.null(output_file)){
      readr::write_delim(meta_data_df, path=output_file, delim = "\t", append=FALSE, col_names = TRUE)
    }
  }
  
  return(meta_data_df)
}
