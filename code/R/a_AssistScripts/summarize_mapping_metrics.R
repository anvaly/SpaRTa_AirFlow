
library(tidyverse)

data_dir <- "XXXX"
filename_out <- paste0(data_dir, "PDAC_metrics_aggregated")
filename_out_ffpe <- paste0(filename_out, "_probes.csv")
filename_out_ff <- paste0(filename_out, "_poly_a.csv")

meta_file <- paste0(data_dir, "full_pdac_meta_data.txt")
meta_df <- read_delim(file=meta_file, delim="\t") %>% filter(grepl("probes", Protocol)) %>% filter(grepl("PDAC", Sample_Name))

for(ii in (5:nrow(meta_df))){
  data_dir <- meta_df[[ii, "FullPath"]]
  filename_in <- paste0(data_dir, "metrics_summary.csv")
  if(ii == 5){
    metrics_ffpe_df <- read_csv(file = filename_in) 
  }else{
    temp_df <- read_csv(file = filename_in)
    metrics_ffpe_df <- add_row(metrics_ffpe_df, temp_df)
  }
}

write_csv(metrics_ffpe_df, file=filename_out_ffpe)
cat("Saved File:", filename_out_ffpe, "\n")

meta_FF_df <- read_delim(file=meta_file, delim="\t") %>% filter(grepl("poly", Protocol))
for(ii in (1:nrow(meta_FF_df))){
  data_dir <- meta_FF_df[[ii, "FullPath"]]
  filename_in <- paste0(data_dir, "metrics_summary.csv")
  if(ii == 1){
    metrics_ff_df <- read_csv(file = filename_in)
  }else{
    temp_df <- read_csv(file = filename_in)
    metrics_ff_df <- add_row(metrics_ff_df, temp_df)
  }
}
write_csv(metrics_ff_df, file=filename_out_ff)
