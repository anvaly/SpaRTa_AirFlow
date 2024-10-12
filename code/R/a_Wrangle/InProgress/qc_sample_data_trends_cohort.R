# Author: Anna Lyubetskaya. Date: 22-05-06

# Evaluate data trends within ST datasets at a cohort level

# Proposed approach to promiscuous genes:
# - Measure the % of spots that express each gene in each sample splitting them by location on the slide (under, adjacent, and outside the tissue section).
# --- Let's call this measure spot occupancy.
# --- Promiscuous genes will be among genes that are identified in a lot of spots; genes that are highly expressed will be also identified in a lot of spots because they diffuse more readily.
# - Calculate mean spot occupancy for each gene and each slide area by protocol and compare gene occupancy between the two protocols.
# --- If a gene is identified everywhere simply because it's highly expressed in PDAC, both protocols will reflect that.
# --- Probe-specific effects will be apparent in our FFPE-probes data but not FF-polyA data. Averaging across all samples allows us to get away from sample/area-specific effects.
# - Identify genes that are expressed in a lot of spots across multiple samples in FFPE protocol but not FF protocol.
# --- These are our most-likely unspecific probes with high diffusion rate.
# - Use unspecific probes as a covariate to our SCT normalization.
# --- These genes correspond to unspecific probes that managed to get everywhere inside and outside of the tissue section and thus potentially reflect diffusion rate rather than true biology.
# --- In order to get away from overnormalizing to one specific gene, we will have a cumulative metric reflecting all these unspecific probes together rather than only one of them, similar to % mito and ribo content.


## ENVIRONMENT ----


# Allow piping throughout the package
`%>%` <- magrittr::`%>%`

source("code/utils/utils_stats.R")
source("code/utils/utils_in_out.R")
source("code/utils/utils_tibble.R")
source("code/utils/utils_ggplot.R")

source("code/R/Utils/utils_10X_in_out.R")


## PARAMETERS ----


# Cohort name
cohort_name <- "PDAC.*"

# Percent of spots outside of tissue that need to express a gene to have it excluded
percent_positive_threshold <- 95


## PATHS ----


# Input file with sample meta data
meta_file <- "XXXX"

# Output folder
input_path_init <- "XXXX"

# Create output folder
output_path_init <- "XXXX"
dir.create(output_path_init, showWarnings = FALSE)

# File with probe offtarget information
probe_file <- "XXXX"


## INGEST DATA ----


# Read in offtarget probe data
probe_df <- read_probe_offtargets(probe_file)

# Read sample file
cohort_meta_df <- readr::read_delim(file=meta_file, delim="\t")  %>% 
  dplyr::filter(grepl(cohort_name, Sample_Name) & Quality_Status == TRUE & !grepl("RTLFBC", Sample_Name))

cohort_meta_ffpe <- cohort_meta_df %>% dplyr::filter(Protocol == "FFPE-probes")
cohort_meta_ffpe_v2 <- cohort_meta_df %>% dplyr::filter(Protocol == "FFPE-probes_v2")
cohort_meta_ff <- cohort_meta_df %>% dplyr::filter(grepl("FF-polyA", Protocol))

# Identify files with the pre-processed sample-level data trends
file_list_FFPE <- paste0(input_path_init, cohort_meta_ffpe$Sample_Name, "/table_mean_v_occupancy_", cohort_meta_ffpe$Sample_Name, ".txt")
file_list_FF <- paste0(input_path_init, cohort_meta_ff$Sample_Name, "/table_mean_v_occupancy_", cohort_meta_ff$Sample_Name, ".txt")
file_list_FFPE_v2 <- paste0(input_path_init, cohort_meta_ffpe_v2$Sample_Name, "/table_mean_v_occupancy_", cohort_meta_ffpe_v2$Sample_Name, ".txt")

# Read all data trend stands into a single tibble
data_df_ffpe <- read_dir2file_my(file_list_FFPE) %>%
  dplyr::mutate(Sample = gsub(input_path_init, "", File),
                Sample = gsub("/table.+", "", Sample)) %>%
  dplyr::select(-File) %>%
  dplyr::mutate(Protocol = "FFPE",
                HasOffTarget = ifelse(grepl(":", Symbol), TRUE, FALSE),
                Symbol = gsub(":.+", "", Symbol),
                Protocol_UnderTissue = paste0(Protocol, UnderTissue))

data_df_ffpe_v2 <- read_dir2file_my(file_list_FFPE_v2) %>%
  dplyr::mutate(Sample = gsub(input_path_init, "", File),
                Sample = gsub("/table.+", "", Sample)) %>%
  dplyr::select(-File) %>%
  dplyr::mutate(Protocol = "FFPE_v2",
                HasOffTarget = ifelse(grepl(":", Symbol), TRUE, FALSE),
                Symbol = gsub(":.+", "", Symbol),
                Protocol_UnderTissue = paste0(Protocol, UnderTissue))

data_df_ff <- read_dir2file_my(file_list_FF) %>%
  dplyr::mutate(Sample = gsub(input_path_init, "", File),
                Sample = gsub("/table.+", "", Sample)) %>%
  dplyr::select(-File) %>%
  dplyr::mutate(Protocol = "FF",
                HasOffTarget = ifelse(grepl(":", Symbol), TRUE, FALSE),
                Symbol = gsub(":.+", "", Symbol),
                Protocol_UnderTissue = paste0(Protocol, UnderTissue))

data_df <- dplyr::bind_rows(data_df_ffpe, data_df_ff, data_df_ffpe_v2)

# Summarize cohort trends by gene and protocol
gene_df <- data_df %>%
  dplyr::group_by(Symbol, Protocol, UnderTissue, Protocol_UnderTissue) %>%
  dplyr::summarise(PercentPositive_Med = round(median(PercentPositive), 1),
                   PercentPositive_SD = round(sd(PercentPositive), 1),
                   ExpressionMeanLogNormUMI_Med = round(median(ExpressionMeanLogNormUMI), 2),
                   ExpressionMeanLogNormUMI_SD = round(sd(ExpressionMeanLogNormUMI), 2),
                   ExpressionMeanLogUMI_Med = round(median(ExpressionMeanLogUMI), 2),
                   ExpressionMeanLogUMI_SD = round(sd(ExpressionMeanLogUMI), 2),
                   HasOffTarget = unique(HasOffTarget)) %>%
  dplyr::ungroup()


## VISUALIZE DATA ----

# Normal ffpe vs ff
# Loop through each of the three spot metrics
for(value in c("ExpressionMeanLogNormUMI_Med", "ExpressionMeanLogUMI_Med", "PercentPositive_Med")){
  
  # Create a wide table of gene-level data trends
  gene_wide_df <- gene_df %>%
    df_long2wide_my(rows="Symbol", cols="Protocol_UnderTissue", value=value) %>%
    dplyr::mutate(HasOffTarget = Symbol %in% probe_df$probe_id) %>%
    # dplyr::mutate_if(is.numeric, dplyr::funs(tidyr::replace_na(., 0))) %>%
    tidyr::drop_na()
  
  # Create a scatter plot comparing FF and FFPE protocols for each segment of the image: under, adjacent, outside of the tissue
  p_list <- list()
  for(p in c("OutsideTissue", "AdjacentTissue", "UnderTissue")){
    
    # Calculate linear regression and residuals, identify outliers
    gene_wide_df <- lm_residual_outliers(gene_wide_df, paste0("FF", p), paste0("FFPE", p), name=p, sd_num=4)
    
    # Factorize categories
    gene_wide_df[paste0("IsOutlier", p)] <- factor(gene_wide_df[[paste0("IsOutlier", p)]], 
                                              levels=c(paste0("FF", p, "Outlier"), paste0("FFPE", p, "Outlier"), "Fit"))
    
    # Comparison between FF and FFPE data of mean occupancy for each gene outside of tissue
    p_list[[p]] <- create_scatter_plot_my(gene_wide_df, x_label=paste0("FF", p), y_label=paste0("FFPE", p), 
                                 fill_label=paste0("IsOutlier", p), size=1, filename=NULL, stroke=0.1, do_fit=TRUE) + 
      geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=0.5)
    
    if(value == "PercentPositive_Med"){
      p_list[[p]] <- p_list[[p]] + 
        xlim(0,100) + ylim(0,100)
    }

  }

  # Write the comparison to file
  filename <- paste0(output_path_init, "table_FF_v_FFPE_", value, ".txt")
  readr::write_delim(gene_wide_df, filename, delim="\t")
  
  # Write figures to file
  filename <- paste0(output_path_init, "scatter_FF_v_FFPE_", value)
  write_plot2file_my(patchwork::wrap_plots(p_list, nrow=1, num_col=4), filename, num_row=1, num_col=4)

  
  # Define promiscuous probes
  gene_select_df <- gene_wide_df %>%
    dplyr::filter(IsOutlierAdjacentTissue == "FFPEAdjacentTissueOutlier" |
                    IsOutlierUnderTissue == "FFPEUnderTissueOutlier" |
                    IsOutlierOutsideTissue == "FFPEOustideTissueOutlier")
  
  # Write the promiscuous genes / probes to file
  filename <- paste0(output_path_init, "table_promiscuous_", value, ".txt")
  readr::write_delim(gene_select_df, filename, delim="\t")
  
}

# Normal ffpe vs ffpe_V2
# Loop through each of the three spot metrics
for(value in c("ExpressionMeanLogNormUMI_Med", "ExpressionMeanLogUMI_Med", "PercentPositive_Med")){
  
  # Create a wide table of gene-level data trends
  gene_wide_df <- gene_df %>%
    df_long2wide_my(rows="Symbol", cols="Protocol_UnderTissue", value=value) %>%
    dplyr::mutate(HasOffTarget = Symbol %in% probe_df$probe_id) %>%
    # dplyr::mutate_if(is.numeric, dplyr::funs(tidyr::replace_na(., 0))) %>%
    tidyr::drop_na()
  
  # Create a scatter plot comparing FF and FFPE protocols for each segment of the image: under, adjacent, outside of the tissue
  p_list <- list()
  for(p in c("OutsideTissue", "AdjacentTissue", "UnderTissue")){
    
    # Calculate linear regression and residuals, identify outliers
    gene_wide_df <- lm_residual_outliers(gene_wide_df, paste0("FFPE_v2", p), paste0("FFPE", p), name=p, sd_num=4)
    
    # Factorize categories
    gene_wide_df[paste0("IsOutlier", p)] <- factor(gene_wide_df[[paste0("IsOutlier", p)]], 
                                                   levels=c(paste0("FFPE_v2", p, "Outlier"), paste0("FFPE", p, "Outlier"), "Fit"))
    
    # Comparison between FF and FFPE data of mean occupancy for each gene outside of tissue
    p_list[[p]] <- create_scatter_plot_my(gene_wide_df, x_label=paste0("FFPE_v2", p), y_label=paste0("FFPE", p), 
                                          fill_label=paste0("IsOutlier", p), size=1, filename=NULL, stroke=0.1, do_fit=TRUE) + 
      geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=0.5)
    
    if(value == "PercentPositive_Med"){
      p_list[[p]] <- p_list[[p]] + 
        xlim(0,100) + ylim(0,100)
    }
    
  }
  
  # Write the comparison to file
  filename <- paste0(output_path_init, "table_FFPE_v2_v_FFPE_", value, ".txt")
  readr::write_delim(gene_wide_df, filename, delim="\t")
  
  # Write figures to file
  filename <- paste0(output_path_init, "scatter_FFPE_v2_v_FFPE_", value)
  write_plot2file_my(patchwork::wrap_plots(p_list, nrow=1, num_col=4), filename, num_row=1, num_col=4)
  
  
  # Define promiscuous probes
  gene_select_df <- gene_wide_df %>%
    dplyr::filter(IsOutlierAdjacentTissue == "FFPEAdjacentTissueOutlier" |
                    IsOutlierUnderTissue == "FFPEUnderTissueOutlier" |
                    IsOutlierOutsideTissue == "FFPEOustideTissueOutlier")
  
  # Write the promiscuous genes / probes to file
  filename <- paste0(output_path_init, "table_promiscuous_FFPE_v2_FFPE_", value, ".txt")
  readr::write_delim(gene_select_df, filename, delim="\t")
  
}

## ffpe v2 vs ff ##
# Loop through each of the three spot metrics
for(value in c("ExpressionMeanLogNormUMI_Med", "ExpressionMeanLogUMI_Med", "PercentPositive_Med")){
  
  # Create a wide table of gene-level data trends
  gene_wide_df <- gene_df %>%
    df_long2wide_my(rows="Symbol", cols="Protocol_UnderTissue", value=value) %>%
    dplyr::mutate(HasOffTarget = Symbol %in% probe_df$probe_id) %>%
    # dplyr::mutate_if(is.numeric, dplyr::funs(tidyr::replace_na(., 0))) %>%
    tidyr::drop_na()
  
  # Create a scatter plot comparing FF and FFPE protocols for each segment of the image: under, adjacent, outside of the tissue
  p_list <- list()
  for(p in c("OutsideTissue", "AdjacentTissue", "UnderTissue")){
    
    # Calculate linear regression and residuals, identify outliers
    gene_wide_df <- lm_residual_outliers(gene_wide_df, paste0("FF", p), paste0("FFPE_v2", p), name=p, sd_num=4)
    
    # Factorize categories
    gene_wide_df[paste0("IsOutlier", p)] <- factor(gene_wide_df[[paste0("IsOutlier", p)]], 
                                                   levels=c(paste0("FF", p, "Outlier"), paste0("FFPE_v2", p, "Outlier"), "Fit"))
    
    # Comparison between FF and FFPE data of mean occupancy for each gene outside of tissue
    p_list[[p]] <- create_scatter_plot_my(gene_wide_df, x_label=paste0("FF", p), y_label=paste0("FFPE_v2", p), 
                                          fill_label=paste0("IsOutlier", p), size=1, filename=NULL, stroke=0.1, do_fit=TRUE) + 
      geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=0.5)
    
    if(value == "PercentPositive_Med"){
      p_list[[p]] <- p_list[[p]] + 
        xlim(0,100) + ylim(0,100)
    }
    
  }
  
  # Write the comparison to file
  filename <- paste0(output_path_init, "table_FF_v_FFPE_v2_", value, ".txt")
  readr::write_delim(gene_wide_df, filename, delim="\t")
  
  # Write figures to file
  filename <- paste0(output_path_init, "scatter_FF_v_FFPE_v2_", value)
  write_plot2file_my(patchwork::wrap_plots(p_list, nrow=1, num_col=4), filename, num_row=1, num_col=4)
  
  
  # Define promiscuous probes
  gene_select_df <- gene_wide_df %>%
    dplyr::filter(IsOutlierAdjacentTissue == "FFPEv2AdjacentTissueOutlier" |
                    IsOutlierUnderTissue == "FFPEv2UnderTissueOutlier" |
                    IsOutlierOutsideTissue == "FFPEv2OustideTissueOutlier")
  
  # Write the promiscuous genes / probes to file
  filename <- paste0(output_path_init, "table_promiscuous_ffpe_v2_vs_ff", value, ".txt")
  readr::write_delim(gene_select_df, filename, delim="\t")
  
}

# # Loop through each of the three spot metrics
# for(value in c("ExpressionMeanLogNormUMI_Med", "ExpressionMeanLogUMI_Med", "PercentPositive_Med")){
# 
#   # Create a wide table of gene-level data trends
#   gene_wide_df <- gene_df %>%
#     df_long2wide_my(rows="Symbol", cols="Protocol_UnderTissue", value=value) %>%
#     dplyr::mutate(HasOffTarget = Symbol %in% probe_df$probe_id) %>%
#     dplyr::filter(HasOffTarget == TRUE) %>%
#     tidyr::drop_na()
#   # Create a wide table of gene-level data trends
#   gene_wide_df <- gene_df %>%
#     df_long2wide_my(rows="Symbol", cols="Protocol_UnderTissue", value=value) %>%
#     dplyr::mutate(HasOffTarget = Symbol %in% probe_df$probe_id) %>%
#     # dplyr::mutate_if(is.numeric, dplyr::funs(tidyr::replace_na(., 0))) %>%
#     tidyr::drop_na()
# 
#   # Create a scatter plot comparing FF and FFPE protocols for each segment of the image: under, adjacent, outside of the tissue
#   p_list <- list()
#   for(p in c("OutsideTissue", "AdjacentTissue", "UnderTissue")){
# 
#     # Calculate linear regression and residuals, identify outliers
#     gene_wide_df <- lm_residual_outliers(gene_wide_df, paste0("FF", p), paste0("FFPE", p), name=p, sd_num=4)
# 
#     # Factorize categories
#     gene_wide_df[paste0("IsOutlier", p)] <- factor(gene_wide_df[[paste0("IsOutlier", p)]],
#                                                    levels=c(paste0("FF", p, "Outlier"), paste0("FFPE", p, "Outlier"), "Fit"))
# 
#     # Comparison between FF and FFPE data of mean occupancy for each gene outside of tissue
#     p_list[[p]] <- create_scatter_plot_my(gene_wide_df, x_label=paste0("FF", p), y_label=paste0("FFPE", p),
#                                           fill_label=paste0("IsOutlier", p), size=1, filename=NULL, stroke=0.1, do_fit=TRUE) +
#       geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=0.5)
# 
#     if(value == "PercentPositive_Med"){
#       p_list[[p]] <- p_list[[p]] +
#         xlim(0,100) + ylim(0,100)
#     }
# 
#   }
# }