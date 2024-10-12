library(tidyverse)
head_dir <- "XXXX"
probe_files <- dir(head_dir, "probe_set.csv", recursive = TRUE)

#initialize the dataframe
data <- read_csv(paste0(head_dir, probe_files[1]), comment = "#")

#loop through, and left join
for(ii in 2:length(probe_files)){
	file_name <- paste0(head_dir, probe_files[ii])
	temp_data <- read_csv(paste0(head_dir, probe_files[ii]), comment = "#")
	data <- full_join(data, temp_data, by = c("gene_id", "probe_seq", "probe_id"))
}

logic_mat <- data[,4:77]
rowSums(logic_mat)
