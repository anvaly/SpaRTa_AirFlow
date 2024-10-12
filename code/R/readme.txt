TESTING!
TESTING!

## 10X data analysis


# INITIAL DATA PROCESSING - a_Wrangle folder


1. Sync data from NGS360 to Stash
  - Script: code/R/a_Wrangle/10X_ngs_s3_to_stash.R

2. Re-run 10X filtering procedure using CellBender
  - Script: code/R/a_Wrangle/10X_CellBender.R

3. Create QC visualizations using 10X reports
  - Script: code/R/a_Wrangle/cohort_qc_10X.R

4. Create a normalized, filtered Seurat object from 10X folder data
  - The same script has 3 versions: for SC, ST, and Perturb data
  - Script: code/R/a_Wrangle/10X_to_Seurat_ST.R
  - Script: code/R/a_Wrangle/10X_to_Seurat_SC.R
  - Script: code/R/a_Wrangle/10X_to_Seurat_Perturb.R

4. Create post-filtering QC visualizations using the final Seurat object
  - Script: code/R/a_Wrangle/cohort_qc_Seurat.R

5. Update meta data information from file in a set of Seurat RDS objects
  - Script: code/R/a_Wrangle/Seurat_meta_update.R
  
6. Extract gene lists from representative experiments for reference
  - Script: code/R/a_Wrangle/Seurat_extract_gene_lists.R
  
7. Plot a set of signatures for a set of Seurat RDS objects
  - Script: code/R/a_Wrangle/Seurat_plot_set_singatures.R


# INITIAL DATA PROCESSING - b_Cluster_Biomarkers folder


1. Perform clustering and find markers at a variety of resolutions
  - Script: code/R/b_Cluster_Biomarkers/Seurat_cluster_ST.R

2. Annotate markers with signatures at preferred resolution


# PATHOLOGY - f_Pathology_integration folder


1. Cell counts
2. Compartment integration
