# 10X data analysis with ST focus: functionality overview

## Sparta Automation Pipeline

Running of the Sparta pipeline has been automated via Airlfow, where users can trigger a run by passing in parameters to the script(s).

Directions: Access Airlfow via the URL below and then enter the Sparta airflow credentials. Next, click on "sparta_automationETL". In the upper right corner, click on the "Play" button and then click "Trigger DAG w/ config". Enter the necessary parameters explained below:
- email: Email address to receive a notifcation of when the pipeline is completed. Must be a BMS email ending in "bms.com" (required)
- download_from_s3: true or false - If data should be downloaded from the NGS360 S3 bucket to a location in stash (required)
- input_path: Path to the factor sheet that should be used. Note: if downloading from S3, then the factor sheet should reference S3 locations, otherwise it should reference stash locations (required)
- output_folder_name: Provide a name for the folder (optional) where NGS360 S3 data will be downloaded to and where other specified script(s) will output to. The time stamp will be concatenated to this name. If no name is provided, then the folder name will just be the time stamp. The pipeline automatically outputs to this location: /XXXX/<output_folder_name><time_stamp> 


To see which script(s) are currently be ran by this automated pipeline, please see the script "run.sh". There you will find an array "scripts_to_run" that will list the script(s) in the order of their execution.


Airflow server: https://XXXX

Please reach out to Anna Lyubetskaya for Sparta airflow credentials.


Note: a link to the pipeline's logs will be included in the email notifications. A user must be logged into Airflow to view the logs.



## General rules

1. All but one scripts work from Stash and write to Stash
2. All but one scripts ingest Seurat objects as RDS files from Stash
3. Scripts that make major changes to Seurat objects also write said Seurat objects to Stash to serve as input to other scripts
4. Small assist files are in this repo in the data folder, everything else is on Stash


## Data management


1. Sync data from NGS360 to Stash
  - Script: code/R/a_Wrangle/10X_ngs_s3_to_stash.R


## Early data QC


1. Re-run 10X filtering procedure using CellBender
  - Script: code/R/a_Wrangle/10X_CellBender.R

2. Create QC visualizations from 10X reports
  - Script: code/R/a_Wrangle/qc_cohort_10X.R

3. Analyze gene expression levels and gene ocurrence in spots under and outside of tissues; identify promiscuous genes
  - Script: code/R/a_Wrangle/qc_sample_data_trends.R


## Main data processing


1. Create a normalized, filtered Seurat object from 10X folder data
  - The same script has 3 versions: for SC, ST, and Perturb data
    - Script: code/R/a_Wrangle/10X_to_Seurat_ST.R
    - Script: code/R/a_Wrangle/10X_to_Seurat_SC.R
    - Script: code/R/a_Wrangle/10X_to_Seurat_Perturb.R
  
2. Perform clustering and find markers at a variety of resolutions
  - Script: code/R/b_Cluster_Biomarkers/Seurat_cluster_ST.R
  - Script: code/R/b_Cluster_Biomarkers/Seurat_cluster_SC.R
  - Script: code/R/b_Cluster_Biomarkers/Seurat_cluster_Perturb.R
  - All 3 scripts work using a utility set of functions: code/R/b_Cluster_Biomarkers/Seurat_cluster_utils.R 

3. Annotate DE markers using signatures and pathway analysis
  - Script: code/R/b_Cluster_Biomarkers/Seurat_markers_annotate.R


## Late data QC


1. Create post-filtering QC visualizations using the final Seurat object
  - Script: code/R/a_Wrangle/cohort_qc_Seurat.R


## Visualization scripts


1. Plot a set of signatures for a set of Seurat RDS objects
  - Script: code/R/a_Wrangle/Seurat_plot_set_singatures.R
2. Visualize cluster biomarkers from FindMarkers
  - Script: code/R/b_Cluster_Biomarkers/Seurat_cluster_heatmaps.R


## Assist scripts


1. Update meta data information from file in a set of Seurat RDS objects
  - Script: code/R/a_Wrangle/Seurat_meta_update.R
  
2. Extract gene lists from representative experiments for reference
  - Script: code/R/a_Wrangle/Seurat_extract_gene_lists.R


## Pathology analysis


1. Integrate cell counts with Visium data
  - Script: code/R/f_Pathology_integration/Seurat_pathology_integrate_cell_counts.R
2. Compartment integration
  - Script: code/R/f_Pathology_integration/Seurat_pathology_integrate_compartments.R
3. Export Seurat spots to an XML format compatible with HALO
  - Script: code/R/f_Pathology_integration/Seurat_to_HALO.R


## Integrate data cohort


1. Integrate a set of Seurat objects using merge, CCA, or RPCA
  - Script: code/R/b_Cluster_Biomarkers/Seurat_integrate_cohort.R
