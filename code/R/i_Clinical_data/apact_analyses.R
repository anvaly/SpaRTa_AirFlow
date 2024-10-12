# Authors: Anna Lyubetskaya. Date: 23-08-08
# Load and wrangle APACT data, create survival plots


## ENVIRONMENT ----


# Allow piping throughout the package
`%>%` <- magrittr::`%>%`


## PARAMETERS ----


## PATHS ----


apact_rds <- "XXXX"
gene_file <- "data/import/References/ensembl_human_mouse_homologs_230111.txt"

output_folder <- "XXXX"


## INGEST DATA ----


# Read the combined RDS from Andy B
data_apact <- readRDS(apact_rds)
apact_tpm <- data_apact$APACT$log2TPM
apact_clin <- data_apact$APACT$meta

rownames(apact_tpm)


# Read gene annotation
gene_df <- readr::read_delim(gene_file, delim="\t") %>%
  dplyr::select(Ensembl_gene_ID_human, Symbol_human) %>%
  unique() %>%
  tidyr::drop_na()

# Find unique gene symbols
gene_count <- table(gene_df$Symbol_human)
gene_list <- names(gene_count)[gene_count == 1]

# Gene list without gene symbol duplicates (because APACT only has gene symbols in this source)
gene_df <- gene_df %>%
  dplyr::filter(Symbol_human %in% gene_list & Symbol_human %in% rownames(apact_tpm))


## PARSE CLINICAL COLUMNS ----


# Sort through clinical columns and remove obviously unnecessary
clin_cols <- colnames(apact_clin)
cols_remove <- grepl("GSVA_.*|RNAseqExpr_.*|multiqc_.*|scaledCopyNumber_.*|copyNumber_.*|HLI_.*|Subtypes_.*|Purist_.*|.*Probability.*", clin_cols)
clin_cols <- clin_cols[!cols_remove]
apact_clin <- apact_clin[clin_cols]

# Identify columns with a lot of empty cells
clin_cols_filled <- c()
for(col in clin_cols){
  if(length(which(apact_clin[[col]] == "")) <= nrow(apact_clin) / 100 &&
     length(which(is.na(apact_clin[[col]]) | is.null(apact_clin[[col]]))) <= nrow(apact_clin) / 2){
    clin_cols_filled <- c(clin_cols_filled, col)
  }
}

clin_cols <- clin_cols_filled
apact_clin <- apact_clin[clin_cols]

# Overall response
cols_os <- clin_cols[grepl("OS_.*", clin_cols)]
cols_os

# This study has disease free survival not progression free survival
cols_dfs <- clin_cols[grepl("DFS_.*", clin_cols)]
cols_dfs

# PDL1
cols_pdl1 <- clin_cols[grepl("PDL1_.*", clin_cols)]
cols_pdl1

for(col in cols_pdl1){
  print(col)
  print(table(apact_clin[[col]]))
}

# Trim down to the remaining columns
clin_cols <- setdiff(clin_cols, c(cols_os, cols_dfs, cols_pdl1,
                                  "to_exclude", "ARM", "SEX", "AGE"))
apact_clin <- apact_clin[clin_cols]


for(col in c("to_exclude", "ARM", "OS_AVAL", "OS_CNSR",
             "DFS_AVAL", "DFS_CNSR", "PDL1_tumor", "CHHIST",
             "CHTNMT", "CHTNMN", "CHTNMM")){
  print(col)
  print(table(data_apact$APACT$meta[[col]]))
}


## PREPARE CLINICAL DATA ----


apact_clin <- data_apact$APACT$meta
apact_clin <- apact_clin[!cols_remove]

# Format the clinical data to spec
apact_clin_df <- apact_clin %>%
  dplyr::select(-USUBJID) %>%
  dplyr::rename(USUBJID = base_ID,
                TreatmentArm = ARM,
                OS = OS_AVAL,
                OS.CNSR = OS_CNSR,
                DFS = DFS_AVAL,
                DFS.CNSR = DFS_CNSR,
                TumorPDL1percent = PDL1_tumor,
                TumorType_Long = CHHIST) %>%
  dplyr::filter(to_exclude == "no" & TreatmentArm != "ScreenFailure") %>%
  dplyr::mutate(STUDYID = "APACT",
                clinFilePath = apact_rds,
                TNM_Stage = paste0(CHTNMT, CHTNMN, CHTNMM),
                TumorType_Short = "PDAC",
                TumorType_Organ = "Pancreas") %>%
  dplyr::select(USUBJID, TreatmentArm, SEX, AGE, STUDYID,
                OS, OS.CNSR, DFS, DFS.CNSR, TumorPDL1percent,
                TumorType_Long, TumorType_Short, TumorType_Organ,
                TNM_Stage, STUDYID, clinFilePath, Accession_ID) %>%
  tidyr::drop_na()

filename <- paste0(output_folder, "/APACT_clinical_data.txt")
readr::write_delim(apact_clin_df, filename, delim="\t")


## CREATE ESET ----


# Create an annotated data frame object of gene annotations
featureData <- new("AnnotatedDataFrame", data = data.frame(row.names = gene_df$Symbol_human,
                                                           Ensembl_gene_ID = gene_df$Ensembl_gene_ID_human))

# Create an annotated data frame object of sample phenotypes
phenoData <- new("AnnotatedDataFrame", data = apact_clin_df %>% 
                   # tibble::column_to_rownames("USUBJID") %>%
                   as.data.frame())

# Create an eSet object
data_eset <- Biobase::ExpressionSet(assayData = as.matrix(apact_tpm)[gene_df$Symbol_human, apact_clin_df$USUBJID],
                                    phenoData = phenoData,
                                    featureData = featureData)


## SURVIVAL ANALYSIS ----


# Genes of interest to perform survival analysis on
genes_interest <- toupper(c("Cdk6", "Cmpk1", "Ddit4", "Hk2", "Hnf4a", "Pmepa1", "Rac1", "Smad3"))

# Symbol to rownames
apact_tpm_df <- tibble::as_tibble(apact_tpm[gene_df$Symbol_human, apact_clin_df$USUBJID], rownames="Symbol") %>%
  tibble::column_to_rownames("Symbol")

for(g in genes_interest){
  HR <- summary(survival::coxph(survival::Surv(apact_clin_df[["OS"]], event = apact_clin_df[["OS.CNSR"]]) ~ as.numeric(apact_tpm_df[g,]) + apact_clin_df$TreatmentArm ))$coefficients[2]
  pval <- summary(survival::coxph(survival::Surv(apact_clin_df[["OS"]], event = apact_clin_df[["OS.CNSR"]]) ~ as.numeric(apact_tpm_df[g,]) + apact_clin_df$TreatmentArm ))$coefficients[5]
  
  cat(g, HR, pval, "\n")
}


## FUNCTION ----


write_plot2file_my <- function(in_plot, filename, num_row=1, num_col=1){
  ## Write plot to file
  
  # Derive file name based on the extension
  device <- "png"
  filename_full <- paste(filename, ".", device, sep="")
  
  # Calculate the size of the figure based on the number of facets
  width <- 5 * ceiling(num_col / 1.5)
  height <- 4 * ceiling(num_row / 1.5)
  
  # Write figure to file
  png(file=filename_full, width=width, height=height, units="in", res=600)
  dev.control(displaylist="enable")
  show(in_plot)
  dev.off()
  
  #+fig.width=width, fig.height=height
  show(in_plot)
  
  cat("Figure created:", filename_full, "\n")
}


## SURVIVAL ANALYSIS CURVES ----


# https://rpkgs.datanovia.com/survminer/reference/ggsurvplot.html
# https://www.theanalysisfactor.com/the-six-types-of-survival-analysis-and-challenges-in-learning-them/
# https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/
# https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html
# https://cran.r-project.org/web/packages/survminer/vignettes/Playing_with_fonts_and_texts.html

indication <- "PDAC"

# Join clinical and mutation data
for(endpoint in c("OS", "DFS")){
  # Time to event variable
  time <- endpoint
  # Censure boolean variable: 1 = lost to follow-up, i.e. the patient is alive
  event <- paste0(endpoint, ".CNSR")
  
  # For each gene group of interest, plot a survival curve
  for(g in genes_interest){
    
    # Gene profile stats
    gene_mean <- mean(as.numeric(apact_tpm_df[g,]))
    gene_sd <- sd(as.numeric(apact_tpm_df[g,]))
    
    cat(gene_mean, gene_sd)
    
    # Gene expression profile
    gene_exp_df <- t(apact_tpm_df[g,]) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("USUBJID")
    
    # Categorize expression into low, medium, high
    gene_exp_df["ExpStatus"] <- ""
    gene_exp_df[gene_exp_df[[g]] >= gene_mean + gene_sd, "ExpStatus"] <- "High"
    gene_exp_df[gene_exp_df[[g]] > gene_mean - gene_sd & gene_mean + gene_sd > gene_exp_df[[g]], "ExpStatus"] <- "Med"
    gene_exp_df[gene_exp_df[[g]] <= gene_mean - gene_sd, "ExpStatus"] <- "Low"
    table(gene_exp_df["ExpStatus"])
    
    # Join clinical data with gene-specific mutation data, assign mutation status, and select an endpoint
    data_joined <- apact_clin_df  %>%
      dplyr::inner_join(gene_exp_df, by="USUBJID") %>%
      dplyr::mutate(event = !!rlang::sym(event) == 0, time = !!rlang::sym(endpoint), GeneGroup = g) %>%
      dplyr::select(event, time, GeneGroup, ExpStatus, TreatmentArm, USUBJID, TumorType_Short)
    
    # Count the number of patients in each group of interest: treatment arm + mutation status of the gene group of interest
    data_joined_stat <- data_joined %>%
      dplyr::group_by(GeneGroup, ExpStatus, TreatmentArm, TumorType_Short) %>%
      dplyr::summarize(n_patients = dplyr::n_distinct(USUBJID)) %>%
      dplyr::mutate(Endpoint = endpoint)
    
    data_joined_stat
    
    # Write group sizes to a file
    readr::write_delim(data_joined_stat, path=paste0(output_folder, "/patient_counts.txt"), 
                       delim = "\t", append=TRUE, col_names = TRUE)
    
    # If patient groups have sufficient number of patients, draw survival curves
    if(min(data_joined_stat$n_patients) > 0){
      # Perform a survival fit
      fit <- survival::survfit(survival::Surv(time, event) ~ ExpStatus + TreatmentArm, data = data_joined)
      
      # Plot a survival curve
      figure_title <- paste0(indication, " patients ", time, " stratified by ", g)
      survival_plot <- survminer::ggsurvplot(fit, data = data_joined,
                                             conf.int = TRUE, conf.int.alpha = 0.05, pval = TRUE,
                                             risk.table = TRUE, risk.table.height=0.3, risk.table.y.text=FALSE,
                                             palette="Paired", legend = "right", title = figure_title)
      
      # Save a survival curve to file
      filename <- paste0(output_folder, paste("/survival", indication, time, g, sep="_"))
      write_plot2file_my(print(survival_plot), filename, num_row=2, num_col=4)
    }
    
  }
}
