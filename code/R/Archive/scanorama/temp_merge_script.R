
`%>%` <- magrittr::`%>%`
library(Seurat)
source("./code/R/Utils/utils_10X_in_out.R")
source("./code/R/Utils/utils_10X_object.R")

input_path <- "XXXX"
seurat_path <- "XXXX"
all_seurat_files <- dir(seurat_path, full.names = TRUE)

meta_df <- readr::read_delim(paste0(input_path, "full_pdac_meta_data.txt"), delim = "\t")
selected <-meta_df %>% dplyr::filter(Quality_Status == TRUE & Protocol != "FFPE-probes_v2")

cohort_regex <- stringr::str_c(selected$Sample_Name, collapse = "|")

seurat_file_names <- all_seurat_files[grepl(cohort_regex, all_seurat_files)]

misc_field_list <- c("user.Sample_Name", "user.Group_Name", "user.Tissue", 
                     "user.Date", "user.Area", "user.Project_ID", "user.Primary", 
                     "user.Protocol", "user.Organism", "user.pt.size.factor", 
                     "user.Clustering", "user.Slide") #,"Pathology.Tissue",
                     # "Pathology.Exocrine",
                     # "Pathology.Benign Glands",
                     # "Pathology.Blood",
                     # "Pathology.Luminal Debris",
                     # "Pathology.Blood Vessels",
                     # "Pathology.Stroma",
                     # "Pathology.Epithelium",
                     # "Pathology.AnnotationTotal",
                     # "Pathology.Total",
                     # "Pathology.Exocrine.percent",
                     # "Pathology.Benign Glands.percent",
                     # "Pathology.Blood.percent",
                     # "Pathology.Luminal Debris.percent",
                     # "Pathology.Blood Vessels.percent",
                     # "Pathology.Stroma.percent",
                     # "Pathology.Epithelium.percent",
                     # "Pathology.None.percent",
                     # "Pathology.Group")


# Ingest a set of RDS Seurat objects288
start_time <- proc.time()
seurat_list <- read_rds_list_simple_my(seurat_file_names)
proc.time() - start_time

for(name in names(seurat_list)){
  # names(seurat_list[[name]]@images) <- name
  seurat_list[[name]]@meta.data <- seurat_list[[name]]@meta.data[which(!grepl("SCT_snn_res.", names(seurat_list[[name]]@meta.data)))]
  seurat_list[[name]]@active.assay <- "Spatial"
  seurat_list[[name]][["SCT"]] <- NULL
}
# Move misc parameters to meta data slot
seurat_list <- seurat_misc_to_meta(seurat_list, misc_field_list=misc_field_list)

#
data_seurat_merged <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)], project = "pdac_cohort_ffpe", merge.data = TRUE)
# save_filename = "XXXX"
# saveRDS(data_seurat_merged, file = save_filename, compress = FALSE)

data_seurat_merged@active.assay <- "Spatial"
# data_seurat_merged[["SCT"]] <- NULL

save_filename = "XXXX"

saveRDS(data_seurat_merged, file = save_filename, compress = FALSE)
