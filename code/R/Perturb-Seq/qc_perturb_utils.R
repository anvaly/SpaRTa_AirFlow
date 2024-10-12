# Author: Anna Lyubetskaya. Date: 20-11-10


ingest_dataset <- function(file_name, meta_df, col_types=NULL){
  ## Ingest a standard 10X output file across a cohort
  
  # A list of full paths to summary metrics file
  files_loc <- paste0(meta_df$FullPath, file_name)
  
  # Read in all QC files
  data_df <- read_dir2file_my(files_loc, delim=",", in_cols=col_types) %>%
    dplyr::mutate(File = gsub(file_name, "", File)) %>%
    dplyr::inner_join(meta_df, by=c("File" = "FullPath")) %>%
    dplyr::select(-File)
  
  return(data_df)
}


crispr_analysis_barplots <- function(data_df, output_folder, file_name, pvalue_threshold=0.05){
  ## Generate various barplots for perturbation_efficiencies_by_feature.csv or perturbation_efficiencies_by_target.csv
  
  f <- gsub(".csv", "", file_name) 
  
  data_loc1_df <- data_df %>%
    dplyr::select(dplyr::any_of(c("Sample_Name", "Perturbation", "Target Gene", "Log2 Fold Change", "p Value", "Cells with Perturbation"))) %>%
    dplyr::mutate(IsSignificant = `p Value` <= pvalue_threshold,
                  Name = paste(`Perturbation`, `Cells with Perturbation`))
  
  filename <- paste0(output_folder, f, "1")
  create_bar_plot_my(data_loc1_df, x_label="Perturbation", y_label="Log2 Fold Change", fill_label="IsSignificant", 
                     facet_var=c("Sample_Name", "fixed"), filename=filename, 
                     labels=c("Perturbation", "Log2 Fold Change", "Log2 Fold Change of all Perturbations"), reorder_x=TRUE)
  
  
  y_values <- c("Mean UMI Count Among Cells with Perturbation",
                "Mean UMI Count Among Cells with Non-Targeting Guides")
  
  data_loc2_df <- data_df %>%
    dplyr::filter(`p Value` <= pvalue_threshold) %>%
    dplyr::select(dplyr::all_of(c("Sample_Name", "Perturbation", y_values))) %>%
    df_wide2long_my(key="Parameter", val="Value", start_col=3) %>%
    dplyr::mutate(Parameter = gsub("Mean UMI Count Among Cells with ", "", Parameter))
  
  filename <- paste0(output_folder, f, "2")
  create_bar_plot_my(data_loc2_df, x_label="Perturbation", y_label="Value", fill_label="Parameter", position="dodge",
                     facet_var=c("Sample_Name", "free_y"), filename=filename, 
                     labels=c("Perturbation", "Mean UMI", "Mean UMI for Perturbed Cells and Non-Targeting Cells for Significant Perturbations"), reorder_x=FALSE)
  
  y_values <- c("Cells with Perturbation", "Cells with Non-Targeting Guides")
  
  data_loc3_df <- data_df %>%
    dplyr::filter(Perturbation %in% unique(data_loc2_df$Perturbation)) %>%
    dplyr::select(dplyr::all_of(c("Sample_Name", "Perturbation", y_values))) %>%
    df_wide2long_my(key="Parameter", val="Value", start_col=3) %>%
    dplyr::mutate(Parameter = gsub("Mean UMI Count Among Cells with ", "", Parameter))
  
  filename <- paste0(output_folder, f, "3")
  create_bar_plot_my(data_loc3_df, x_label="Perturbation", y_label="Value", fill_label="Parameter", position="dodge",
                     facet_var=c("Sample_Name", "free_y"), filename=filename, 
                     labels=c("Perturbation", "Cell Count", "Cell Count for Perturbed Cells and Non-Targeting Cells for Significant Perturbations"), reorder_x=FALSE)
}
