
`%>%` <- magrittr::`%>%`
ab_subtypes_ccle <- readRDS("XXXX")
full_metadata_ccle <- readRDS("XXXX")
output_path <- "XXXX"
dir.create(output_path, showWarnings = FALSE)

meta_ccle <- full_metadata_ccle$meta
norm_exp <- full_metadata_ccle$log2TPM

cell_line_rows <- rownames(meta_ccle)
panc_row_names <- rownames(ab_subtypes_ccle)

row_indices <- match(panc_row_names, cell_line_rows)

panc_metadata <- meta_ccle[row_indices,]
panc_norm_exp <- norm_exp[,row_indices]


col_select <- c(1, 2, 6, 7, 11, 12, 13, 16:19, 21, 28)
select_panc_metadata <- panc_metadata[,col_select]

DepMap_ID <-panc_row_names
final_panc_metadata <- tidyr::as_tibble(cbind(DepMap_ID, select_panc_metadata, ab_subtypes_ccle))
readr::write_csv(final_panc_metadata, file = "XXXX")

# Get Specific Gene Expression

# gene_select <- c("KRAS", "GATA6", "GATA4", "HNF1A", "HNF4G", "TP63", "SMAD4", "ONECUT2")
gene_select <- c("GATA6", "TP63", "KRT5", "KRT6A")
colnames(panc_norm_exp) <- select_panc_metadata$stripped_cell_line_name
gene_subset_norm_exp <- panc_norm_exp[gene_select,] %>% t %>% scale

col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#4218F1", "#EEEEEE", "#F14232"))
subtype <- data.frame(Purist_Subtype = final_panc_metadata$Subtypes_PuristN_cohortTSP)
rownames(subtype) <- final_panc_metadata$stripped_cell_line_name
ann_colors <- list(Purist_Subtype =c(Purist_Basal ="#222222", Purist_Classical="#2ca52f"))

png(filename = paste0(output_path, "subtype_marker_heatmap.png"), width = 600, height = 800)
ComplexHeatmap::pheatmap(gene_subset_norm_exp,
                         col = col_fun,
                         clustering_method = "ward.D2",
                         annotation_row = subtype,
                         annotation_colors = ann_colors)
dev.off()
# pheatmap::pheatmap(gene_subset_norm_exp, scale = 'row', )
