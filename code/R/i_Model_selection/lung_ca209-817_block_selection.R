# Author: Anna Lyubetskaya. Date: 23-04-17
# Sample selection for spatial profiling from CA209-817
# Confluence page: https://biodoc.pri.bms.com:8443/display/CT/CA209-817+%28CM817%29+-+NSCLC+Biomarker+Project

# Approach:
# --- Only select blocks with more than 10 sections remaining.
# --- Necrosis <= 10%; Tumor% >= 10%; Area size at about half of Visium area.
# --- Collect mutation status


## SETUP ENVIRONMENT ----


# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_in_out.R")
source("code/utils/utils_ggplot.R")
source("code/utils/utils_heatmap.R")


## PARAMETERS ----


# Inventory
inventory_sheet <- "FFPE"
inventory_col <- "Potential sections remaining"
inventory_val <- c("10", "15", "20", "25", "30", "35", "40", "45", "50", "11-20", "21-30", "30+")

# Block size
size_sheet <- "MOS958_CD3+CD8+GrB"
size_col <- "Tissue Area (mm2)"
size_val <- 6.5*6.5/2  # Half Visium slide

# Block composition
composition_sheet <- "MOS003_H&E"
composition_col1 <- "% Necrosis"
composition_val1 <- 10
composition_col2 <- "% Cancer by Area"
composition_val2 <- 10


## PATHS ----


## Input files
input_path_inventory <- "XXXX"
input_path_ihc <- "XXXX"
input_path_mutation_manifest <- "XXXX"
input_path_mutations <- c("XXXX",
                          "XXXX",
                          "XXXX")
input_path_cd8 <- c("XXXX")
input_clinical <- "XXXX"


# Output folder
output_path <- "XXXX"
dir.create(output_path, showWarnings=FALSE)

# Write the log
filename_log <- paste0(output_path, "log.txt")


## INGEST DATA ----


## INVENTORY


# Read in inventory
inventory_df <- readxl::read_excel(input_path_inventory, sheet=inventory_sheet)

# Print options for number of sections left
# table(inventory_df[[inventory_col]])
write(paste("Number of patients in inventory:", length(unique(inventory_df$Subject)), "\n"), filename_log)

# Select blocks with sufficient amount of sections left
inventory_df <- inventory_df %>%
  dplyr::filter(!!rlang::sym(inventory_col) %in% inventory_val)

write(paste("Number of patients in inventory after filter:", length(unique(inventory_df$Subject)), "\n"), filename_log, append=TRUE)


## BLOCK SIZE


# Read in inventory
size_df <- readxl::read_excel(input_path_ihc, sheet=size_sheet, skip=1)

# Print options for number of sections left
# table(round(as.numeric(size_df[[size_col]])))
write(paste("Number of blocks measured:", length(unique(size_df$`Subject ID`)), "\n"), filename_log, append=TRUE)

# Select blocks with sufficient amount of sections left
size_df <- size_df %>%
  dplyr::filter(as.numeric(!!rlang::sym(size_col)) >= size_val)

# Ensure columns are numeric
size_df[[size_col]] <- round(as.numeric(size_df[[size_col]]), 2)

write(paste("Number of blocks measured after filtering:", length(unique(size_df$`Subject ID`)), "\n"), filename_log, append=TRUE)


## BLOCK COMPOSITION


# Read in inventory
composition_df <- readxl::read_excel(input_path_ihc, sheet=composition_sheet, skip=1)

# Print options for number of sections left
# table(composition_df[[composition_col1]])
# table(composition_df[[composition_col2]])
write(paste("Number of blocks evaluated:", length(unique(composition_df$`Subject ID`)), "\n"), filename_log, append=TRUE)

# Select blocks with sufficient amount of sections left
composition_df <- composition_df %>%
  dplyr::filter(as.numeric(!!rlang::sym(composition_col1)) <= composition_val1 & as.numeric(!!rlang::sym(composition_col2)) >= composition_val2)

# Ensure columns are numeric
composition_df[[composition_col1]] <- round(as.numeric(composition_df[[composition_col1]]), 2)
composition_df[[composition_col2]] <- round(as.numeric(composition_df[[composition_col2]]), 2)

write(paste("Number of blocks evaluated after filtering:", length(unique(composition_df$`Subject ID`)), "\n"), filename_log, append=TRUE)


# MUTATIONS


# Read in mutation manifest
mutations_meta_df <- readxl::read_excel(input_path_mutation_manifest)

# Read in mutations
mutations_list <- list()

for(f in input_path_mutations){
  mutations_list[[f]] <- readr::read_delim(f, delim="\t") %>% 
    dplyr::rename(USUBJID = `SUBJECT ID`, VARIANT_TYPE = `VARIANT-TYPE`,
           STATUS = `SOMATIC STATUS/FUNCTIONAL IMPACT`, HGVSP = `SV-PROTEIN-CHANGE`, CNA_TYPE = `CNA-TYPE`, CNA_RATIO =`CNA-RATIO`) %>% 
    dplyr::select(USUBJID, VARIANT_TYPE, GENE, STATUS, HGVSP, CNA_TYPE, CNA_RATIO) %>%
    dplyr::mutate(USUBJID = as.character(USUBJID))
}

# Merge and filter mutation information
mutations_df <- dplyr::bind_rows(mutations_list) %>%
  dplyr::filter(STATUS %in% c("known", "likely")) %>%
  dplyr::mutate(Count = 1)

write(paste("Number of patients with mutation:", length(unique(mutations_df$USUBJID)), "\n"), filename_log, append=TRUE)
write(paste("Number of mutation:", nrow(mutations_df), "\n"), filename_log, append=TRUE)


## CD8


# Infiltration status
cd8_df <- readr::read_delim(input_path_cd8, delim="\t")


## CLINICAL

# Clinical data
clinical_df <- readr::read_delim(input_clinical, delim=",")


## INTEGRATE DATA ----


# Merge all data
data_df <- inventory_df %>%
  dplyr::select(dplyr::all_of(c("Site", "Subject", inventory_col))) %>%
  dplyr::inner_join(size_df %>%
                      dplyr::select(dplyr::all_of(c("Subject ID", size_col))), by=c("Subject" = "Subject ID")) %>%
  dplyr::inner_join(composition_df %>%
                      dplyr::left_join(cd8_df, by=c("Mosaic ID" = "ID")) %>%
                      dplyr::select(dplyr::all_of(c("Subject ID", composition_col1, composition_col2, "CD8 TOPOLOGY", "PERCENT CD8 TOTAL"))), by=c("Subject" = "Subject ID")) %>%
  dplyr::mutate(`Partner Subject ID` = paste0("CA209817-", gsub("^0+", "", Site), "-", gsub("^0+", "", Subject))) %>%
  dplyr::inner_join(mutations_meta_df %>%
                      dplyr::select(dplyr::all_of(c("Partner Subject ID", "Accession Number", "Gender", "Date Received",
                                                    "Pathology %TN", "TMB Score"))), by="Partner Subject ID") %>%
  dplyr::inner_join(mutations_df, by=c("Accession Number" = "USUBJID")) %>%
  dplyr::inner_join(clinical_df, by=c("Partner Subject ID" = "USUBJID")) %>%
  dplyr::select(-`Partner Subject ID`, -'Accession Number') %>%
  dplyr::mutate(MutationOfInterest = GENE %in% c("KRAS", "STK11", "KEAP1", "TP53", "CDKN2A", "NFE2L2")) %>%
  dplyr::rename(CancerAreaPerc = `% Cancer by Area`) %>%
  unique()

data_df[["TMB Score"]] <- round(as.numeric(data_df[["TMB Score"]]), 2)


sort(table(data_df$GENE))

write(paste("Number of patients with mutation after filtering:", length(unique(data_df$Subject)), "\n"), filename_log, append=TRUE)
write(paste("Number of mutation after filtering:", nrow(data_df), "\n"), filename_log, append=TRUE)


# Write to file
filename <- paste0(output_path, "table_patients_select_by_mutation.txt")
readr::write_delim(data_df, filename, delim="\t")


# Summarize mutation information by patient
data_pat_df <- data_df %>% 
  dplyr::arrange(GENE) %>% 
  dplyr::group_by(Site, Subject, `Potential sections remaining`, `Tissue Area (mm2)`, `% Necrosis`, CancerAreaPerc, Gender, 
                  `Date Received`, `Pathology %TN`, `TMB Score`, `CD8 TOPOLOGY`, `PERCENT CD8 TOTAL`,
                  AGE, OS, SMKNG, CRFHIST, DSCELLSP, DSENTST, ECOGPS) %>% 
  dplyr::summarise(Alterations = paste(unique(GENE), collapse="; "))


# Write to file
filename <- paste0(output_path, "table_patients_select.txt")
readr::write_delim(data_pat_df, filename, delim="\t")

# Count all mutation pairs
mutation_pair_list <- list()
for(i in 1:nrow(data_pat_df)){
  mutation_pat_list <- strsplit(data_pat_df[[i, "Alterations"]], "; ")[[1]]
  if(length(mutation_pat_list) >= 2){
    mutation_pair_list[[i]] <- tidyr::unite(data.frame(t(combn(mutation_pat_list, 2))), "Mut_pair") %>%
      dplyr::mutate(Subject = data_pat_df[[i, "Subject"]])
  }
}

mut_pairs_df <- dplyr::bind_rows(mutation_pair_list) %>%
  dplyr::group_by(Mut_pair) %>%
  dplyr::summarise(Patient_list = paste(unique(sort(Subject)), collapse="; "),
                   Count = dplyr::n_distinct(Subject)) %>%
  dplyr::arrange(desc(Count), Mut_pair)


## VISUALIZE DATA ----


# Boxplot of mutations by patient
p <- create_bar_plot_my(data_df, x_label="GENE", y_label="Count", fill_label="Subject", 
                        filename=NULL, labels=c("Mutations", "Subject number", "Mutations by subject"), reorder_x=TRUE)

filename <- paste0(output_path, "bar_mut_by_pat")
write_plot2file_my(p, filename, num_row=1, num_col=3.5)


# Boxplot of mutation pairs by patient
p <- create_bar_plot_my(mut_pairs_df %>%
                          dplyr::filter(Count > 1), x_label="Mut_pair", y_label="Count", fill_label="Patient_list", 
                        filename=NULL, labels=c("Mutation pairs", "Subject number", "Mutation pairs by subject"), reorder_x=TRUE)

filename <- paste0(output_path, "bar_mutpairs_by_pat")
write_plot2file_my(p, filename, num_row=1, num_col=3.5)


# Create a heatmap of mutations by patient
params <- list(cell_value = "Count",
               row_label = "GENE",
               col_label = "Subject",
               distance = "pearson",
               row_annotation = c("MutationOfInterest"),
               col_annotation = NULL,  # c("CancerAreaPerc"),
               range = c(0, 1),
               colors = c("white", "royalblue4"),
               column_order = sort(unique(data_df$Subject)),
               row_order = sort(unique(data_df$GENE)))

filename <- paste0(output_path, "hm_mut_by_pat.png")
create_heatmap_my(data_df, params, row_list=NULL, col_list=NULL, col_meta_df=data_df, row_meta_df=data_df, 
                  filename=filename, width=8, height=12)
