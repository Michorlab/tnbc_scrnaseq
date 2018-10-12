## this script preprocessses the normal breast single cells (Gao et al 2017 Nat Comm) with the same pipeline as the TNBC
## cells (script preprocess.R), to more accurately compare the two sets in a correlation map (as done in the script correlation_map.R)

library(here)
source(here("code", "libraries.R"))

mat_ct <- readRDS(here("data", "mat_ct.RDS")) # normalized expression values (as resulting from analysis1.R script)
mappings <- readRDS(here("data", "mappings.RDS")) # genes info
qc <- read.table(here("data", "qc_original.txt"), sep = "\t", stringsAsFactors = FALSE) # quality control


## read normal single cell breast data from Nick Navin
nick_normalize <- read.table(here("data", "RSEM_TPM_240_NormalCells.txt"), sep = "\t", header = TRUE)

# match gene names with the ones we have
mat_ct <- mat_ct[match(intersect(rownames(nick_normalize), rownames(mat_ct)), rownames(mat_ct)),]
nick_normalize <- nick_normalize[match(intersect(rownames(nick_normalize), rownames(mat_ct)), rownames(nick_normalize)),]
all.equal(rownames(mat_ct), rownames(nick_normalize))
mappings <- mappings[match(rownames(mat_ct), mappings$hgnc_symbol),]
all.equal(mappings$hgnc_symbol, rownames(mat_ct))

# rename normal cells
colnames(nick_normalize) <- paste("normal", c(1:ncol(nick_normalize)), sep = "")


## qc on normals
qc_normals <- matrix(NA, nrow = ncol(nick_normalize), ncol = ncol(qc))
qc_normals <- as.data.frame(qc_normals)
colnames(qc_normals) <- colnames(qc)
qc_normals$Sample <- colnames(nick_normalize)
qc_normals$experiment <- "Exp4"
qc_normals$pool_H12 <- 0
rownames(qc_normals) <- colnames(nick_normalize)


## scater normals initialization
fd_norm <- new("AnnotatedDataFrame", data = as.data.frame(mappings))
rownames(fd_norm) <- fd_norm$hgnc_symbol
pd_norm <- new("AnnotatedDataFrame", data = as.data.frame(qc_normals))
rownames(pd_norm) <- rownames(qc_normals)
sceset_norm <- SingleCellExperiment(assays = list(tpm = as.matrix(nick_normalize)), 
                                    rowData = as.data.frame(mappings),
                                    colData = as.data.frame(qc_normals))

## monocle TPM object (to transform to relative counts)
min_thresh_tpm <- 1
HSMM_norm <- newCellDataSet(tpm(sceset_norm), 
                            phenoData = pd_norm, 
                            featureData = fd_norm, 
                            expressionFamily = tobit(), lowerDetectionLimit = min_thresh_tpm)
HSMM_norm <- detectGenes(HSMM_norm, min_expr = min_thresh_tpm)

# monocle add relative counts
rpc_matrix_norm <- relative2abs(HSMM_norm)
min_thresh_log_tpm <- 0.1
HSMM_norm <- newCellDataSet(as(as.matrix(rpc_matrix_norm), "sparseMatrix"),
                            phenoData = phenoData(HSMM_norm),
                            featureData = featureData(HSMM_norm),
                            lowerDetectionLimit = min_thresh_log_tpm,
                            expressionFamily = negbinomial.size())
HSMM_norm <- estimateSizeFactors(HSMM_norm)
HSMM_norm <- estimateDispersions(HSMM_norm)


## scater: fill in other expression slots
assay(sceset_norm, "exprs") <- log2(tpm(sceset_norm) + 1)
assay(sceset_norm, "monocle") <- as.matrix(exprs(HSMM_norm))
assay(sceset_norm, "log2_tpm") <- log2(tpm(sceset_norm) + 1)
assay(sceset_norm, "log2_monocle") <- as.matrix(log2(exprs(HSMM_norm) + 1))
assay(sceset_norm, "log2_tpm_median") <- scale(assays(sceset_norm)$exprs, scale = FALSE)

# monocle: normalize by dividing with size factors estimated by Monocle and logged2, and add to scater
exprs_filtered <- t(t(exprs(HSMM_norm)/pData(HSMM_norm)$Size_Factor))
nz_genes <- which(exprs_filtered != 0)
exprs_filtered[nz_genes] <- log2(exprs_filtered[nz_genes] + 1)
assay(sceset_norm, "norm_log2_monocle") <- as.matrix(exprs_filtered)
rm(exprs_filtered)
rm(nz_genes)


## scran normalization done on the Monocle counts
sceset_norm_mon <- sceset_norm
counts(sceset_norm_mon) <- assays(sceset_norm)$monocle
sceset_norm_mon <- computeSumFactors(sceset_norm_mon, positive = TRUE)
if (length(which(sizeFactors(sceset_norm_mon) == 0)) > 0) {
  sceset_norm_mon <- sceset_norm_mon[,-which(sizeFactors(sceset_norm_mon) == 0)]
  set_exprs(sceset_norm_mon, "log2_tpm_median") <- scale(assays(sceset_norm_mon)$exprs, scale = FALSE)
}
summary(sizeFactors(sceset_norm_mon))
sceset_norm_mon <- scater::normalize(sceset_norm_mon)


## remove additional sources of unwanted variation
housekeepers <- read.table(here("data", "housekeepers.txt"), sep = "\t")
colnames(housekeepers) <- "gene"
housekeepers <- unique((housekeepers))
housekeepers$index <- match(housekeepers$gene, rowData(sceset_norm_mon)$hgnc_symbol)
if (length(which(is.na(housekeepers$index)) > 0))
  housekeepers <- housekeepers[-which(is.na(housekeepers$index)), ]

# RUVg on scran Monocle counts
ruvg1_norm_mon_scran <- RUVg(exprs(sceset_norm_mon), housekeepers$index, k = 1, isLog = 1)
assay(sceset_norm_mon, "ruvg1_mon_scran") <- ruvg1_norm_mon_scran$normalizedCounts


## combine the two matrices (normals and TNBC)
mat_norm <- assays(sceset_norm_mon)$ruvg1_mon_scran
mat_complete <- cbind(mat_norm, mat_ct)


