## this script processes the TPM data to obtain the normalized data (mat_norm, and its phenodata, pd_norm) 

library(here)
source(here::here("code", "libraries.R"))
source(here::here("code", "funcs.R"))

## function to intersect multiple vectors
intersect_all <- function(a,b,...){
  Reduce(intersect, list(a,b,...))
}

## read data
tpm.rsem <- read.table(here::here("data", "tpm_original.txt"), sep = "\t")
counts.rsem <- read.table(here::here("data", "counts_original.txt"), sep = "\t")
qc <- read.table(here::here("data", "qc_original.txt"), sep = "\t", stringsAsFactors = FALSE) # quality control
mappings <- readRDS(here::here("data", "mappings.RDS")) # genes info

## scater intialize
min_thresh_log_tpm <- 0.1
sceset_all <- SingleCellExperiment(assays = list(counts = as.matrix(counts.rsem), 
                                                 tpm = as.matrix(tpm.rsem),
                                                 exprs = log2(as.matrix(tpm.rsem) + 1), 
                                                 expressed = log2(tpm.rsem + 1) > min_thresh_log_tpm),
                                   colData = qc, 
                                   rowData = mappings)

## scater quality control
sceset_all <- calculateQCMetrics(sceset_all, exprs_values = "exprs", 
                                 detection_limit = min_thresh_log_tpm,
                                 cell_controls = list("regular" = which(qc$pool_H12 == 0), "pool" = which(qc$pool_H12 == 1)))
if (length(which(is.na(as.numeric(sceset_all$uniquely_mapped_percent)))) > 0)
  sceset_all <- sceset_all[, -which(is.na(as.numeric(sceset_all$uniquely_mapped_percent)))]

# total reads (library size)
reads.drop <- isOutlier(as.numeric(sceset_all$total_reads), nmads = 4, type = "lower", log = TRUE)
hist(log10(sceset_all$total_reads), breaks = 100, main = "", col = "grey80", xlab = expression(log[10]~"library size"), ylab = "frequency", cex.lab = 1.4, cex.axis = 1.4)
abline(v = max(log10(sceset_all$total_reads[which(reads.drop == 1)])), col = "blue", lwd = 2, lty = 2)

# number of expressed genes
feature.drop <- isOutlier(sceset_all$total_features_by_exprs, nmads = 4, type = "lower", log = TRUE)
hist(log10(sceset_all$total_features_by_exprs), breaks = 100, main = "", col = "grey80", xlab = expression(log[10]~"number of expressed genes"), ylab = "frequency", cex.lab = 1.4, cex.axis = 1.4)
abline(v = max(log10(sceset_all$total_features_by_exprs[which(feature.drop == 1)])), col = "blue", lwd = 2, lty = 2)

# mRNA
colData(sceset_all)$total_tpm <- colSums(assays(sceset_all)$exprs) # total_tpm is the sum of log2(tpm + 1)
mRNA.drop <- isOutlier(sceset_all$total_tpm, nmads = 4, type = "lower", log = TRUE)
hist(log10(sceset_all$total_tpm), breaks = 100, main = "", col = "grey80", xlab = expression(log[10]~"total mRNA"), ylab = "frequency", cex.lab = 1.4, cex.axis = 1.4)
abline(v = max(log10(sceset_all$total_tpm[which(mRNA.drop == 1)])), col = "blue", lwd = 2, lty = 2)

keep.samples <- !(reads.drop | feature.drop | mRNA.drop)
keep.samples[which(is.na(keep.samples))] <- FALSE
sceset_all <- sceset_all[ ,keep.samples]


## monocole filter unexpressed genes
pd <- new("AnnotatedDataFrame", data = as.data.frame(colData(sceset_all)))
featureNames(pd) <- rownames(pd)
fd <- new("AnnotatedDataFrame", data = as.data.frame(rowData(sceset_all)))
featureNames(fd) <- rowData(sceset_all)$gene_short_name
# expressionFamily = gaussianff()) not supported from VGAM 1.0.6 onwards
HSMM <- newCellDataSet(assays(sceset_all)$exprs, featureData = fd, phenoData = pd,
                       lowerDetectionLimit = min_thresh_log_tpm, expressionFamily = gaussianff())

# genes unexpressed in less than 5% of cells for all patients
HSMM <- detectGenes(HSMM, min_expr = min_thresh_log_tpm)
expr_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= dim(HSMM)[2]*0.05))

# genes unexpressed in less than 5% of cells for each patient
HSMM126 <- HSMM[, row.names(subset(pData(HSMM), patient == "PT126"))]
HSMM126 <- detectGenes(HSMM126, min_expr = min_thresh_log_tpm)
remove_genes_126 <- row.names(subset(fData(HSMM126), num_cells_expressed < dim(HSMM126)[2]*0.05))
HSMM126 <- HSMM126[setdiff(row.names(HSMM126), remove_genes_126), ]

HSMM39 <- HSMM[, row.names(subset(pData(HSMM), patient == "PT039"))]
HSMM39 <- detectGenes(HSMM39, min_expr = min_thresh_log_tpm)
remove_genes_39 <- row.names(subset(fData(HSMM39), num_cells_expressed < dim(HSMM39)[2]*0.05))
HSMM39 <- HSMM39[setdiff(row.names(HSMM39), remove_genes_39), ]

HSMM58 <- HSMM[, row.names(subset(pData(HSMM), patient == "PT058"))]
HSMM58 <- detectGenes(HSMM58, min_expr = min_thresh_log_tpm)
remove_genes_58 <- row.names(subset(fData(HSMM58), num_cells_expressed < dim(HSMM58)[2]*0.05))
HSMM58 <- HSMM58[setdiff(row.names(HSMM58), remove_genes_58), ]

HSMM81 <- HSMM[, row.names(subset(pData(HSMM), patient == "PT081"))]
HSMM81 <- detectGenes(HSMM81, min_expr = min_thresh_log_tpm)
remove_genes_81 <- row.names(subset(fData(HSMM81), num_cells_expressed < dim(HSMM81)[2]*0.05))
HSMM81 <- HSMM81[setdiff(row.names(HSMM81), remove_genes_81), ]

HSMM84 <- HSMM[, row.names(subset(pData(HSMM), patient == "PT084"))]
HSMM84 <- detectGenes(HSMM84, min_expr = min_thresh_log_tpm)
remove_genes_84 <- row.names(subset(fData(HSMM84), num_cells_expressed < dim(HSMM84)[2]*0.05))
HSMM84 <- HSMM84[setdiff(row.names(HSMM84), remove_genes_84), ]

HSMM89 <- HSMM[, row.names(subset(pData(HSMM), patient == "PT089"))]
HSMM89 <- detectGenes(HSMM89, min_expr = min_thresh_log_tpm)
remove_genes_89 <- row.names(subset(fData(HSMM89), num_cells_expressed < dim(HSMM89)[2]*0.05))
HSMM89 <- HSMM89[setdiff(row.names(HSMM89), remove_genes_89), ]

remove_genes_all <- intersect_all(remove_genes_126, remove_genes_39, remove_genes_58, remove_genes_81, remove_genes_84, remove_genes_89)
keep.genes <- setdiff(rownames(fData(HSMM)), remove_genes_all)
length(keep.genes)
HSMM <- HSMM[keep.genes, ]


## scater recompute quality metrics
sceset <- sceset_all[keep.genes, ]
sceset <- calculateQCMetrics(sceset, exprs_values = "exprs", 
                             detection_limit = min_thresh_log_tpm,
                             cell_controls = list("regular" = which(colData(sceset)$pool_H12 == 0), "pool" = which(colData(sceset)$pool_H12 == 1)))


## monocle relative counts
min_thresh_tpm <- 1
pd_2 <- new("AnnotatedDataFrame", data = as.data.frame(colData(sceset)))
featureNames(pd_2) <- rownames(pd_2)
fd_2 <- new("AnnotatedDataFrame", data = as.data.frame(rowData(sceset)))
featureNames(fd_2) <- rowData(sceset)$gene_short_name
HSMM <- newCellDataSet(assays(sceset)$tpm, phenoData = pd_2, featureData = fd_2, 
                       expressionFamily = tobit(), lowerDetectionLimit = min_thresh_tpm)
HSMM <- detectGenes(HSMM, min_expr = min_thresh_tpm)
colData(sceset)$num_genes_expressed <- pData(HSMM)$num_genes_expressed

# add relative counts
rpc_matrix <- relative2abs(HSMM)
HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                       phenoData = phenoData(HSMM),
                       featureData = featureData(HSMM),
                       lowerDetectionLimit = min_thresh_log_tpm,
                       expressionFamily = negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)


## scater update slots
assays(sceset)$monocle <- as.matrix(exprs(HSMM))
assays(sceset)$log2_tpm <- log2(assays(sceset)$tpm + 1)
assays(sceset)$log2_monocle <- as.matrix(log2(exprs(HSMM) + 1))
assays(sceset)$log2_tpm_median <- t(scale(t(assays(sceset)$exprs), scale = FALSE))

# normalize by dividing with monocle size factors
exprs_filtered <- t(t(exprs(HSMM)/pData(HSMM)$Size_Factor))
nz_genes <- which(exprs_filtered != 0)
exprs_filtered[nz_genes] <- log2(exprs_filtered[nz_genes] + 1)
assays(sceset)$norm_log2_monocle <- as.matrix(exprs_filtered)
rm(exprs_filtered)
rm(nz_genes)


## scran normalization on Monocle counts
sceset_mon_norm <- sceset
assays(sceset_mon_norm)$counts <- assays(sceset)$monocle

# original code (working slightly differently now with the more recent versions of quickCluster and computeSumFactors from scran)
#clusters_mon_norm <- quickCluster(sceset_mon_norm, min.size = 300)
#sceset_mon_norm <- computeSumFactors(sceset_mon_norm, cluster = clusters_mon_norm, sizes = 300, positive = TRUE)
clusters_mon_norm <- quickCluster(sceset_mon_norm, min.size = 200)
sceset_mon_norm <- computeSumFactors(sceset_mon_norm, cluster = clusters_mon_norm, sizes = 200, positive = TRUE)

if (length(which(sizeFactors(sceset_mon_norm) == 0)) > 0) {
  sceset_mon_norm <- sceset_mon_norm[,-which(sizeFactors(sceset_mon_norm) == 0)]
  assays(sceset_mon_norm)$log2_tpm_median <- t(scale(t(assays(sceset_mon_norm)$exprs), scale = FALSE))
}
sceset_mon_norm <- scater::normalize(sceset_mon_norm)
assays(sceset_mon_norm)$exprs <- assays(sceset_mon_norm)$logcounts


## RUVg remove additional sources of unwanted variation using housekeeper genes
housekeepers <- read.table(here("data", "housekeepers.txt"), sep = "\t")
colnames(housekeepers) <- "gene"
housekeepers <- unique((housekeepers))
housekeepers$index <- match(housekeepers$gene, rowData(sceset_mon_norm)$hgnc_symbol)
if (length(which(is.na(housekeepers$index)) > 0))
  housekeepers <- housekeepers[-which(is.na(housekeepers$index)), ]

# works slightly different because of the changes in quickCluster and computeSumFactors from above
ruvg1_mon_scran <- RUVg(assays(sceset_mon_norm)$exprs, housekeepers$index, k = 1, isLog = 1)
assays(sceset_mon_norm)$ruvg1_mon_scran <- ruvg1_mon_scran$normalizedCounts


## store the normalized data in new object
sceset_final <- sceset_mon_norm
assays(sceset_final)$exprs <- assays(sceset_final)$ruvg1_mon_scran
sceset_final <- calculateQCMetrics(sceset_final, 
                                   exprs_values = "exprs", 
                                   detection_limit = min_thresh_log_tpm,
                                   cell_controls = list("regular" = which(colData(sceset_final)$pool_H12 == 0), "pool" = which(colData(sceset_final)$pool_H12 == 1)))
colData(sceset_final)$total_norm_exprs <- colSums(assays(sceset_final)$exprs)

## normalized data to be used further
mat_norm <- assays(sceset_final)$exprs
pd_norm <- colData(sceset_final)

# due to updates/changes in quickCluster and computeSumFactors (at a minimum), the resulting values are slightly different
# than the ones obtained when carrying out the original analysis. to reproduce the results from the paper, we can simply
# read in the normalized data and the phenoData, as preprocessed in the original analysis

mat_norm <- read.table(here("data", "norm_data.txt"), sep = "\t", header = TRUE)
pd_norm <- read.table(here("data", "pd_norm.txt"), sep = "\t", header = TRUE)

