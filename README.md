# Unravelling subclonal heterogeneity and aggressive disease states in TNBC through single-cell RNA-seq

This repository contains materials for the paper [**Unravelling subclonal heterogeneity and aggressive disease states in TNBC through single-cell RNA-seq**](https://www.nature.com/articles/s41467-018-06052-0), by Karaayvaz*, Cristea* et al, *Nature Communications* (2018).

## Scripts
* *analysis1.R*: : performs most of the analyses of the single cell transcriptomic data in the paper, including cycling cells identification, cell types identification (via markers and clustering), t-sne of all cells, clustering of epithelial cells, differential expression between clusters, basal PNAS signature, Lehman signatures, normal signatures, mammaprint signature (PS), zena werb signature (MBS), carlos artega signature (RTS), expression of genes in the metabolism and immunity pathways; it generates the following figs: 1b-e, 2a-c, 3, 4a-b,d, S2, S6-S10, S13-S16

* *correlation_map.R*: creates expression correlation maps for normal (Gao et al 2017 Nat Comm) and TNBC cells; it generatres figs. 2f and S12

* *cnv_single_cell.R*: computes copy number values for single cells (TNBC cells, as well as normal cells from Gao et al 2017 Nat Comm); it generates figs. 2d and S11

* *funcs.R*: helper functions for various aspects of cell type identification, clustering, signatures etc

* *funcs_markers.R*: helper functions to identify cell types based on markers

* *preprocess.R*: preprocesses the raw data to obtain the normalized data (*mat_norm*, and its phenodata, *pd_norm*) 

* *preprocess_normals.R*: preprocessses the normal breast single cells (Gao et al 2017 Nat Comm) with the same pipeline as the TNBC cells (script *preprocess.R*), to more accurately compare the two sets in a correlation map (as done in the script *correlation_map.R*)




## Plots
The scripts and data in this repository can generate almost all figures in the main text and supplement, as follows:

  * Main text: fig.1 (b-e); fig.2 (a-d,f); fig. 3(a-g); fig. 4(a,b,d)
  
  * Supplementary Information: figs. S2, S6-S16


## Data
The following useful data structures are included in the repository:

  * *artega_sig.txt*: list of genes comprising the 354-gene residual tumor signature (RTS) (Balko et al., 2012)
  
  * *basal_signature.txt*: list of genes (positive and negative associations), together with effect size (log fold-change) comprising the basal signature from the Normal Breast Signatures (Lim et al., 2009a)
  
  * *cell_types_tab_S9.txt*: the 1,112 classified cells and their corresponding cell type (as in Table S9 in the Supplement)
  
  * *counts_celltypes_patient_1c.txt*: the number of cells from each different cell type, per patient (as in fig. 1c)
  
  * *counts_original.txt*: counts matrix (as downloaded from GEO)
  
  * *epithelual_cells_cluster.txt*: the 868 epithelial cells and their assigned epithelial cluster (from 1 to 5)
  
  * *fd_norm.txt*: feature information on the 13,280 selected genes, such as: chromosome location (*chromosome_name*), start and end position, GC content, percent dropout (*pct_dropout*), mean expression across cells (*mean_exprs*), and others
  
  * *genes_for_basal_vs_non_basal_tnbc_PNAS.txt*: list of genes for the Intrinsic Basal signature (Sørlie et al., 2001); column *Basal.epithelial.cell.enriched.cluster* from the table
  
  * *housekeepers.txt*: list of 98 housekeeping genes compiled in Tirosh et al., 2016, to be used in data preprocessing, to remove sources of unwanted variation
  
  * *Lehman_signature.txt*: table with all genes characteristic of the 6 Lehmann signatures, together with direction of regulations and number of supporting samples from the original study. The table includes the TNBCtype4 signatures (Lehmann et al., 2016)
  
  * *LP_signature.txt*: list of genes (positive and negative associations), together with effect size (log fold-change) comprising the luminal progenitor signature from the Normal Breast Signatures (Lim et al., 2009)
  
  * *mammaprint_sig_new.txt*: list of genes comprising the 70-gene prognostic signature (PS) (Van’t Veer et al., 2002a), also known as Mammaprint
  
  * *mappings.RDS*: RData with feature information on the 13,280 selected genes, *hgnc_symbol*, *chromosome_name*, *start_position*, *end_position*, *gc*, *gene_short_name* (a subset of the *fd_norm.txt* table)
  
  * *markers_clean.txt*: table with curated cell type markers (Table S7)
  
  * *melanoma_cellcycle.txt*: sets of relevant genes used for determining the cell cycle state of cells, as done in Tirosh et al., 2016
  
  * *ML_signature.txt*: list of genes (positive and negative associations), together with effect size (log fold-change) comprising the mature luminal signature from the Normal Breast Signatures (Lim et al., 2009)

  * *order_samples_cnv.RDS*: RData used for plotting, containing the order of cells in 2d, to be kept the same for 2e
  
  * *original_clustering_epithelial.RDS*: RData with the original cluster assignment (5 clusters) of epithelial cells in the manuscript, done with Monocle
  
  * *original_clustering.RDS*: RData with the original cluster assignment (9 clusters) for cell type identification in the manuscript, done with Monocle
  
  * *pd_ct.RDS*: RData with detailed phenotypic information on the 1,112 classified cells, such as quality control metrics (*uniquely_mapped_percent* - percent of uniquely mapped reads, *insertion_length*, *deletion_length*, *total_reads*, *GC*, and others), *patient*, sequencing batch (*number_batch*), cycling scores (*mel_scores_g1s* - score for the G1S phase, *mel_scores_g2m* - score for the G2M phase), assigned cell type (*cell_types_cl_all*)
  
  * *pd_epith.RDS*: RData with detailed phenotypic information on the 868 epithelial cells, including all features from the object *pd_ct.RDS*, together with the additional features, such as assignments to the Lehman signatures *assignments_lehman_both*, Normal signatures *assignments_normsig_both*, average expression under the MBS (*zenawerb_avg_exprs*), RTS (*artega_avg_exprs*), PS (*mammaprint_avg_exprs*) and Intrinsic Basal (*basal_PNAS_avg_exprs*) signatures
  
  * *pd_norm.txt*
  
  * *qc_original.txt*
  
  * *RSEM_TPM_240_NormalCells.txt*
  
  * *tpm_original.txt*: TPM matrix (as downloaded from GEO)
  
  * *werb_49_metastasis_sig.txt* 
  

## Tables
The tables contain lists of differentially expressed genes between epithelial clusters. *allvs_i_* is a list of differentially expressed genes between the epithelial cells in cluster *i*, and all other epithelial cells.
