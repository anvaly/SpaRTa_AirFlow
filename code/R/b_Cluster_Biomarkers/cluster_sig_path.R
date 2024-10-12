# Author: Anna Lyubetskaya. Date: 23-10-13
# Compare pathology and signature composition for each de novo cluster


## ENVIRONMENT ----


library(Seurat)

# Allow piping throughout the package.
`%>%` <- magrittr::`%>%`

source("code/utils/utils_tibble.R")
source("code/utils/utils_ggplot.R")
source("code/R/Utils/utils_10X_vis.R")


## PARAMETERS ----


# Rename classes to better labels
label_rename <- TRUE


# Name of the analysis to use in folder/file names
# TLS_LN, Mac_in_Fib, Tcell_in_Fib, Bcell_in_Fib
run_name <- "Peng"

# Path to processed Seurat data
# PDAC108_path14_harmonyepi_rpca_sct, PDAC108_path14_5K_harmony
cohort_name <- "PDAC108_path14_5K_harmony"


# Fix order of pathology niches and signatures
path_order <- c("ExoEndo", "IntestineAdj", "BenignEpi", "Tumor", "LuminalNec",
                "NonEpi", "Muscle", "MuscleAdj", "NormalAdj", "Nerve", "Vessel",
                "Blood", "LymphNode", "TLSAggregate", "TLSImmature", "TLSMature", "Adipose")

sig_order <- paste0("sig.", c("PDAC.P19.Acinar","PDAC.P19.Endocrine",
                              "PDAC.P19.Ductal_1", "PDAC.collisson.classical","PDAC.moffitt.basal",
                              # "PDAC.Elyada19.panCAF","PDAC.Elyada19.iCAF","PDAC.Elyada19.myCAF",
                              "PDAC.P19.Fibroblast","PDAC.P19.Stellate","PDAC.P19.Endothelial","PDAC.P19.Macrophage",
                              "PDAC.P19.Tcell","PDAC.P19.Bcell" # "PDAC.CosMx.Mast","PDAC.CosMx.Plasma","PDAC.U.Nervous",
                              # "BMS.Pathway.TGFB","BMS.Pathway.IFNa","BMS.Pathway.IFNg","BMS.Pathway.TNFa","BMS.CL.Hypoxia"
                              ))


# Clustering resolution and order of clusters preferred


cluster_order <- NULL
cols <- NULL

## Harmony full integration, res 0.1
if(cohort_name == "PDAC108_path14_5K_harmony"){
  cluster_order <- c(3, 6, 4, 1, 8, 0, 5, 2, 7, 10, 11, 9)
  cols <- c("#056DB5", "#153C65", "#4CB0E1", "#9A2626", "#96257D", "#E5E5AA",
            "#A0D5B5", "#FF9F2C", "#DF7126", "#AEE0EA", "#D4EEF5", "#939393")
  
  resolution <- "integrated_snn_res.0.1"
}
## Harmony-defined epi niche integrated by RPCA, res 0.4
if(cohort_name == "PDAC108_path14_harmonyepi_rpca_sct"){
  # cluster_order <- c(12,11,3,10,0,4,6,14,9,1,8,5,2,7,13)
  # cols <- c("#056DB5","#153C65",
  #           "#4CB0E1","#D08EB3","#F16666","#9A2626","#96257D","#360F2E","#681D59",
  #           "#E5E5AA","#5E5E39","#358E5B","#FF9F2C","#D66100","#939393")
  
  cluster_order <- c(3, 10, 0, 4, 6, 9)
  cols <- c("#4CB0E1","#D08EB3","#F16666","#9A2626","#96257D","#681D59")
  
  resolution <- "integrated_snn_res.0.4"
}
## Harmony-defined stroma niche integrated by RPCA, res 0.2
if(cohort_name == "PDAC108_path14_harmonystr_rpca_sct"){
  cluster_order <- c(11,7,0,3,1,4,8,5,2,6,10,9)
  resolution <- "integrated_snn_res.0.2"
}
## Harmony-defined immune niche integrated by RPCA, res 0.2
if(cohort_name == "PDAC108_path14_harmonyimm_rpca_sct"){
  cluster_order <- c(10,16,9,0,1,4,3,7,14,2,11,12,13,8,5,6,15)
  resolution <- "integrated_snn_res.0.4"
}

names(cols) <- cluster_order


if(label_rename == TRUE){
  
  path_rename <- c("Exo/Endocrine",
                   "Intestine    ",
                   "Ducts        ",
                   "Tumor        ",
                   "Other        ",
                   "Act. Stroma  ",
                   "Normal Stroma",
                   "Normal Stroma",
                   "Normal Stroma",
                   "Nerve        ",
                   "Normal Stroma",
                   "Other        ",
                   "LN & TLS     ",
                   "LN & TLS     ",
                   "LN & TLS     ",
                   "LN & TLS     ",
                   "Other        ")
  
  sig_rename <- c("Acinar","Endocrine","Ductal", "Tumor classical","Tumor basal",
                  # "panCAF","iCAF","myCAF",
                  "Fibroblast","Stellate","Endothelial","Macrophage",
                  "T cell","B cell" # "Mast","Plasma","Nerves",
                  # "TGFb Pathway","IFNa Pathway","IFNg Pathway","TNFa Pathway","Hypoxia"
                  )
  
  
  ## Harmony full integration, res 0.1
  if(cohort_name == "PDAC108_path14_5K_harmony"){
    res_rename <- c("Acinar\n","Endocrine\n","Ductal\n","Tumor classical\nTumor mixed",
                    "Tumor basal\n","Fibroblast\n","Stellate\nEndothelial",
                    "Immune 1\n","Immune 2\n","Intestine 1\n","Intestine 2\n","Blood\n")
  }
  ## Harmony-defined epi niche integrated by RPCA, res 0.4
  if(cohort_name == "PDAC108_path14_harmonyepi_rpca_sct"){
    res_rename <- c("Acinar","Endocrine",
                    "Ductal","PanIN","Tumor Classical",
                    "Tumor Mixed\nProliferative","Tumor Mixed","Tumor Basal\nTumor Mixed",
                    "Tumor Basal\nHypoxia","Fibroblast","Fibroblast\nHypoxia",
                    "Endothelial","Macrophage","Infiltrated\nfibroblast","Nerve")
  }
  ## Harmony-defined stroma niche integrated by RPCA, res 0.2
  if(cohort_name == "PDAC108_path14_harmonystr_rpca_sct"){
    res_rename <- c("Acinar","Tumor Mixed",
                    "Fibroblast 1","Fibroblast 2","Fibroblast / Hypoxia",
                    "Stellate 1","Stellate 2","Endothelial","Macrophage / Hypoxia",
                    "Early PDAC / IG","IFNa pathway / resistance","Nerve")
  }
  ## Harmony-defined immune niche integrated by RPCA, res 0.2
  if(cohort_name == "PDAC108_path14_harmonyimm_rpca_sct"){
    res_rename <- c("Acinar","Endocrine","Tumor Mixed","Fibroblast",
                    "Fibroblast / IG 1","Fibroblast / IG 2","Fibroblast / IG 3",
                    "Endothelial","Stellate",
                    "Plasma 1","Plasma 2","Plasma 3","B cell / Proliferation",
                    "Macrophage","T cell","Tumor / T cell","Nerves")
  }
  
} else{
  path_rename <- path_order
  sig_rename <- sig_order
  res_rename <- cluster_order
}


## PATHS ----


# Location of pre-processed data
input_path <- paste0("XXXX")

# Output folder
output_path_init <- "XXXX"
output_path <- paste0(output_path_init, run_name, "_", cohort_name, "_rename", label_rename[1], "/")

# Create output folders
dir.create(output_path_init, showWarnings = FALSE)
dir.create(output_path, showWarnings = FALSE)


# A file of additional barcode labels
meta_data_labels <- "XXXX"


## INGEST DATA ----


# Seurat data
data_seurat <- readRDS(input_path)

# Remove clusters if necessary
if(!is.null(cluster_order)){
  barcode_list <- rownames(data_seurat@meta.data[which(data_seurat@meta.data[[resolution]] %in% cluster_order),])
  data_seurat <- subset(data_seurat, cells=barcode_list)
}


# Read in signature meta data and make a long tibble
sig_wide_df <- readr::read_delim(meta_data_labels, delim="\t")

sig_df <- sig_wide_df %>%
  df_wide2long_my(key="Sig_Name", val="Sig_Occurrence")


## EXTRACT META DATA ----


# Extract meta data and add signature assignments
meta_df <- tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate") %>%
  dplyr::inner_join(sig_df, by="Coordinate") %>%
  dplyr::filter(!!rlang::sym(resolution) %in% cluster_order &
                  !!rlang::sym("Pathology.Group") %in% path_order &
                  !!rlang::sym("Sig_Name") %in% sig_order)

# Simplify group names
meta_df[[resolution]] <- as.character(meta_df[[resolution]])
for(i in 1:length(cluster_order)){
  meta_df[[resolution]] <- gsub(paste0("^", cluster_order[i], "$"), res_rename[i], meta_df[[resolution]])
}
for(i in 1:length(path_order)){
  meta_df[["Pathology.Group"]] <- gsub(paste0("^", path_order[i], "$"), path_rename[i], meta_df[["Pathology.Group"]])
}
for(i in 1:length(sig_order)){
  meta_df[["Sig_Name"]] <- gsub(paste0("^", sig_order[i], "$"), sig_rename[i], meta_df[["Sig_Name"]])
}


# Count the size of clusters
cluster_size_df <- meta_df %>%
  dplyr::select(Coordinate, !!rlang::sym(resolution)) %>%
  unique() %>%
  dplyr::group_by(!!rlang::sym(resolution)) %>%
  dplyr::summarise(CountClust = dplyr::n_distinct(Coordinate))

# Count the size of pathology niches within clusters
cluster_path_size_df <- meta_df %>%
  dplyr::select(Coordinate, !!rlang::sym(resolution), Pathology.Group) %>%
  unique() %>%
  dplyr::group_by(!!rlang::sym(resolution), Pathology.Group) %>%
  dplyr::summarise(CountClustPath = dplyr::n_distinct(Coordinate))

# Summarize data by de novo cluster, signature, and pathology niche
cluster_df <- meta_df %>%
  dplyr::group_by(!!rlang::sym(resolution), Pathology.Group, Sig_Name) %>%
  dplyr::summarise(PathCount = dplyr::n_distinct(Coordinate),
                   SigCount = sum(Sig_Occurrence)) %>%
  dplyr::inner_join(cluster_size_df, by=resolution) %>%
  dplyr::inner_join(cluster_path_size_df, by=c(resolution, "Pathology.Group")) %>%
  dplyr::mutate(PathPerc = round(PathCount / CountClust * 100, 1),
                SigPerc = round(SigCount / CountClustPath * 100, 1)) %>%
  dplyr::arrange(desc(PathPerc), desc(SigPerc))


# Print underlying labels and their stats to files
filename <- paste0(output_path, "barcode_cluster_path_sig.txt")
readr::write_delim(tibble::as_tibble(data_seurat@meta.data, rownames="Coordinate") %>%
                     dplyr::select(dplyr::all_of(c("Coordinate", resolution, "Pathology.Group"))) %>%
                     dplyr::inner_join(sig_wide_df, by="Coordinate"), 
                   filename, delim="\t")

filename <- paste0(output_path, "table_cluster_path_sig.txt")
readr::write_delim(cluster_df, filename, delim="\t")


# Fix order of all groups
cluster_df[[resolution]] <- factor(cluster_df[[resolution]], levels=unique(res_rename))
cluster_df[["Pathology.Group"]] <- factor(cluster_df[["Pathology.Group"]], levels=unique(path_rename))
cluster_df[["Sig_Name"]] <- factor(cluster_df[["Sig_Name"]], levels=unique(sig_rename))


# Summarize data by patient and de novo cluster
cluster_pat_df <- meta_df %>%
  dplyr::group_by(user.Block_ID) %>%
  dplyr::mutate(PatientTumorSpotNum = dplyr::n_distinct(Coordinate)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(user.Block_ID, PatientTumorSpotNum, !!rlang::sym(resolution)) %>%
  dplyr::summarise(ClustCount = dplyr::n_distinct(Coordinate),
                   SectionCount = dplyr::n_distinct(user.Sample_Name),
                   ClustCountPerSection = round(ClustCount / SectionCount)) %>%
  dplyr::mutate(ClustPerc = round(ClustCount / PatientTumorSpotNum * 100, 1)) %>%
  dplyr::arrange(desc(ClustPerc))

# Fix order of all groups
cluster_pat_df[[resolution]] <- factor(cluster_pat_df[[resolution]], levels=unique(res_rename))


## CREATE BAR PLOTS ----


# Plot clusters by patients
for(var in c("ClustCount", "ClustCountPerSection", "ClustPerc")){
  
  p <- create_bar_plot_my(cluster_pat_df, x_label="user.Block_ID", y_label=var, fill_label=resolution, 
                          position="stack", filename=NULL, 
                          labels=c("Patient", "De novo cluster membership", 
                                   "Patients by their de novo cluster composition"),
                          reorder_x = FALSE, cols=cols)
  
  filename <- paste0(output_path, "bar_patients_by_cluster_", var, "_", run_name)
  write_plot2file_my(p, filename, num_row=1, 
                     num_col=round(length(unique(cluster_pat_df[["user.Block_ID"]]))/10))
  
}


## CREATE BUBBLE HEATMAP ----


# Bubble heatmap of de novo clusters, pathology composition, and signature occurrence
if(!"bubbleHeatmap" %in% installed.packages()){
  install.packages("bubbleHeatmap")
}

library(bubbleHeatmap)


# Parameters for bubble heatmap
x_name <- "Sig_Name"
y_name <- "Pathology.Group"
z_name <- resolution
value_fill <- "PathPerc"
value_size <- "SigPerc"
fill_threshold <- 10  # Threshold on the value corresponding to the y_name var
scale_threshold <- 50

# Create a series of bubble heatmap
# https://cran.r-project.org/web/packages/bubbleHeatmap/vignettes/Using_Bubbleheatmap.html
z_list <- levels(cluster_df[[z_name]])
bubble_list <- list()
for(i in 1:length(z_list)){
  
  z <- z_list[i]
  
  cluster_loc_df <- cluster_df %>%
    dplyr::filter(!!rlang::sym(z_name) == z) %>%
    dplyr::select(dplyr::all_of(c(x_name, y_name, value_fill, value_size))) %>%
    dplyr::filter(!!rlang::sym(value_fill) >= fill_threshold) %>%
    as.data.frame()
  
  names <- list(paste0("leftLabels", 1:6), paste0("topLabels", 1:10))
  
  colorMat <- cluster_loc_df %>%
    dplyr::select(dplyr::all_of(c(x_name, y_name, value_fill))) %>%
    df_long2wide_my(rows=x_name, cols=y_name, value=value_fill) %>%
    tibble::column_to_rownames(x_name) %>%
    as.matrix()
  
  sizeMat <- cluster_loc_df %>%
    dplyr::select(dplyr::all_of(c(x_name, y_name, value_size))) %>%
    df_long2wide_my(rows=x_name, cols=y_name, value=value_size) %>%
    tibble::column_to_rownames(x_name) %>%
    as.matrix()
  
  colorMat[is.na(colorMat)] <- 0
  sizeMat[is.na(sizeMat)] <- 0
  
  
  colorSeq <- RColorBrewer::brewer.pal(n=6, name="PuRd")
  colorBreaks <- seq(0,scale_threshold,10)
  sizeBreaks <- seq(0,scale_threshold,10)
  
  # Setup bubble heatmap parameters depending on where in the list of plots we are
  if(i == length(z_list)){
    showColorLegend <- T
    showBubbleLegend <- T
  } else{
    showColorLegend <- F
    showBubbleLegend <- F
  }
  
  if(i == 1){
    yTitle <- x_name
  } else{
    yTitle <- F
    rownames(colorMat) <- letters[1:nrow(colorMat)]
    rownames(sizeMat) <- letters[1:nrow(sizeMat)]
  }
  
  # Create plot tree
  bubble_list[[i]] <- bubbleHeatmap::bubbleHeatmap(colorMat, sizeMat, treeName = z,
                                                   leftLabelsTitle = F, showRowBracket = F,
                                                   rowTitle = F, showColBracket = F, colTitle=F,
                                                   plotTitle=z, 
                                                   xTitle=F, yTitle=yTitle,
                                                   legendTitles = c(value_fill, value_size),
                                                   colorSeq = colorSeq,
                                                   colorBreaks = colorBreaks,
                                                   sizeBreaks = sizeBreaks,
                                                   colorLim = c(0, scale_threshold), sizeLim = c(0, scale_threshold),
                                                   showColorLegend = showColorLegend, 
                                                   showBubbleLegend = showBubbleLegend)
  
  # Write figure to file
  filename <- paste0(output_path, "hm_bubble_", gsub("[ /\n]+", "", z), ".png")
  png(file=filename, width=6, height=7, units="in", res=300)
  
  grid.newpage()
  grid.draw(bubble_list[[i]])
  
  dev.off()
  # grid.arrange(grobs = plist, ncol = 2) ## display plot
  # ggsave(file = OutFileName, arrangeGrob(grobs = plist, ncol = 2))  ## save plot
}


# Write figure to file
filename <- paste0(output_path, "joint_hm_bubble_", run_name, ".png")
png(file=filename, width=1.5*length(bubble_list), height=7, units="in", res=300)

grid.newpage()

# Set plot width ratios
norm_plot_ratio <- 1/(length(bubble_list) + 2*2)
big_plot_ratio <- norm_plot_ratio*3
plot_ratio_list <- c(big_plot_ratio, rep(norm_plot_ratio, length(bubble_list)-2), big_plot_ratio)

cowplot::plot_grid(plotlist = bubble_list, align = "hv", ncol = length(bubble_list), 
                   rel_widths = plot_ratio_list)

dev.off()
