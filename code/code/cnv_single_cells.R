## this script computes copy number values for single cells (TNBC cells as well as normal cells from Gao et al 2017 Nat Comm)
## this script generates the following figures
##    fig. 2d
##    fig. S11

library(here)
source(here("code", "libraries.R"))

## read tpm data
tpm.rsem <- read.table(here::here("data", "tpm_original.txt"), sep = "\t")

## read the normalized data (outputs of analysis1.R)
mat_ct <- readRDS(here::here("data","mat_ct.RDS"))
pd_ct <- readRDS(here::here("data", "pd_ct.RDS")) # phenoData
fd_ct <- readRDS(here::here("data", "fd_ct.RDS")) # featureData

## read data of normal single breast cells
nick_normalize <- read.table(here::here("data", "RSEM_TPM_240_NormalCells.txt"), sep = "\t", header = TRUE)

# match gene names with the ones we have
tpm.cnv <- tpm.rsem[match(intersect(rownames(tpm.rsem), rownames(nick_normalize)), rownames(tpm.rsem)),]
nick_normalize <- nick_normalize[match(intersect(rownames(tpm.rsem), rownames(nick_normalize)), rownames(nick_normalize)),]
all.equal(rownames(tpm.cnv), rownames(nick_normalize))

# keep only the cells of interest
tpm.cnv <- tpm.cnv[,match(colnames(mat_ct), colnames(tpm.cnv))]

# transform to log2(TPM+1)
log.tpm.cnv <- log2(tpm.cnv + 1)
log.nick_normalize <- log2(nick_normalize + 1)

# divide by 10 to mimic the library complexity of 1 million
log.tpm.cnv <- log.tpm.cnv/10
log.nick_normalize <- log.nick_normalize/10

# rename normal cells
colnames(log.nick_normalize) <- paste("normal", c(1:ncol(log.nick_normalize)), sep = "")

# keep all genes with expression higher than threshold
threshold_to_keep <- 0.1
log.tpm.cnv <- log.tpm.cnv[(apply(log.tpm.cnv, 1, mean) >= threshold_to_keep),]

# get the mapping coordinates for the filtered genes from the sceset object
if (length(which(is.na(match(rownames(log.tpm.cnv), rownames(mat_ct))))) == 0)
  fd.cnas <- fd_ct[match(rownames(log.tpm.cnv), rownames(mat_ct)),]
all.equal(rownames(fd.cnas), rownames(log.tpm.cnv))
mappings.cnas <- fd.cnas[,c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position")]

# order genes in mappings according to their positions on the chromosomes
chr.st.cnas <- paste(mappings.cnas$chromosome_name, mappings.cnas$start_position, sep = " ")
order.map.cnas <- match(gtools::mixedsort(chr.st.cnas), chr.st.cnas)
mappings.cnas <- mappings.cnas[order.map.cnas, ]
if (length(which(duplicated(mappings.cnas))) > 0)
  mappings.cnas <- mappings.cnas[-which(duplicated(mappings.cnas)), ]
mappings.cnas <- as.data.frame(mappings.cnas) %>% dplyr::filter(chromosome_name %in% c(1:24, "X", "Y", "x", "y"))
mappings.cnas$chromosome_name <- revalue(mappings.cnas$chromosome_name, c("X" = "23", "x" = "23", "Y" = "24", "y" = "24"))
rownames(mappings.cnas) <- mappings.cnas$hgnc_symbol

# order genes in the tpm data according to the odering in mappings
log.tpm.cnv <- log.tpm.cnv[match(mappings.cnas$hgnc_symbol, rownames(log.tpm.cnv)),]
all.equal(rownames(log.tpm.cnv), rownames(mappings.cnas))

# restrict the normal dataset to the same genes
log.nick_normalize <- log.nick_normalize[match(rownames(log.tpm.cnv), rownames(log.nick_normalize)),]
all.equal(rownames(log.nick_normalize), rownames(log.tpm.cnv))

# compute average of normal cells
mean.genes <- apply(log.nick_normalize, 1, mean)

# remove mean normal expression from tumor data
cnv.data <- sweep(log.tpm.cnv, 1, mean.genes)

# add normal epithelial cells to the tumor data
cnv.data <- cbind(log.nick_normalize, cnv.data)

# level off extreme values in data
if (length(which(cnv.data > 3) > 0))
  cnv.data <- cnv.data[-which(cnv.data > 3)]
if (length(which(cnv.data < -3) > 0))
  cnv.data <- cnv.data[-which(cnv.data < -3)]

# compute sliding window
window.size <- 101
sl.avg <- apply(cnv.data, 2, caTools::runmean, k = window.size)
rownames(sl.avg) <- rownames(cnv.data)

# center CNA values
scaled.sl.avg <- scale(sl.avg, scale = FALSE)

# clustering method and metric
clustering_method_columns <- "ward.D"
clustering_distance_columns <- "pearson"


## clustering of all epithelial cells to be used as annotation
pd_epith <- readRDS(here::here("data", "pd_epith.RDS"))
# same as HSMM_allepith_clustering$Cluster in analysis1.R; also same as original_clustering_epithelial.RDS
clustering_allepith <- pd_epith$Cluster

patients_now <- sort(unique(pd_epith$patient))
clusterings_sep_allepith <- list()
for (i in patients_now) {
  clusterings_sep_allepith[[i]] <- clustering_allepith[which(pd_epith$patient == i)]
  names(clusterings_sep_allepith[[i]]) <- rownames(pd_epith)[which(pd_epith$patient == i)]
}


## separate data into normals and epithelials per patient for separate heatmaps
norm_oth_idx <- which(pd_ct$cell_types_cl_all != "epithelial")

mat_normepith_cna <- scaled.sl.avg[,1:ncol(log.nick_normalize)]
scaled.sl.avg <- scaled.sl.avg[,-c(1:ncol(log.nick_normalize))]

all.equal(rownames(pd_ct), colnames(scaled.sl.avg))
mat_normoth_cna <- scaled.sl.avg[,norm_oth_idx]
pd_normoth_ct <- pd_ct[norm_oth_idx,]

mats_ct <- list()
pds_ct <- list()
for (i in 1:length(patients_now)) {
  mats_ct[[i]] <- mat_ct[,pd_ct$patient == patients_now[i]]
  pds_ct[[i]] <- pd_ct[pd_ct$patient == patients_now[i],]
}
names(mats_ct) <- patients_now
names(pds_ct) <- patients_now

# data to plot
mats_epith_cna <- list()
pds_epith_ct <- list()
for (i in 1:length(patients_now)) {
  mats_epith_cna[[i]] <- scaled.sl.avg[,intersect(which(pd_ct$cell_types_cl_all == "epithelial"), which(pd_ct$patient == patients_now[i]))]
  pds_epith_ct[[i]] <- pds_ct[[i]][which(pds_ct[[i]]$cell_types_cl_all == "epithelial"),]
}
names(mats_epith_cna) <- patients_now
names(pds_epith_ct) <- patients_now

# heatmap annotation for normal epithelial cells
ha_cnas_normepith <- HeatmapAnnotation(data.frame(celltype = "normal"), 
                                       col = list(celltype = c("normal" = "violet")),
                                       show_annotation_name = TRUE,
                                       annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 8),
                                       annotation_legend_param = list(list(title_position = "topcenter",
                                                                           title = c("cell type"))),
                                       show_legend = TRUE,
                                       gap = unit(c(2), "mm"))

# heatmap annotation for all epithelial cells with separate clusters
ha_cnas_epith_sep <- list()
for (i in 1:length(patients_now)) {
  
  if (i == 1)
    ha_cnas_epith_sep[[i]] <- HeatmapAnnotation(data.frame(cluster_all = clusterings_sep_allepith[[i]]),
                                                col = list(cluster_all = c("1" = "#ee204d", "2" = "#17806d", "3" = "#b2ec5d", "4" = "#cda4de", "5" = "#1974d2")),
                                                annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 12),
                                                annotation_legend_param = list(list(title_position = "topcenter", title = c(paste("clust_all")))),
                                                show_annotation_name = FALSE,
                                                gap = unit(c(2), "mm"),
                                                show_legend = FALSE)
  
  if (i > 1 && i != 5)
    ha_cnas_epith_sep[[i]] <- HeatmapAnnotation(data.frame(cluster_all = clusterings_sep_allepith[[i]]),
                                                col = list(cluster_all = c("1" = "#ee204d", "2" = "#17806d", "3" = "#b2ec5d", "4" = "#cda4de", "5" = "#1974d2")),
                                                annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 12),
                                                annotation_legend_param = list(list(title_position = "topcenter", title = c(paste("clust_all")))),
                                                show_annotation_name = FALSE,
                                                gap = unit(c(2), "mm"),
                                                show_legend = FALSE)
  
  if (i == 5)
    ha_cnas_epith_sep[[i]] <- HeatmapAnnotation(data.frame(cluster_all = clusterings_sep_allepith[[i]]),
                                                col = list(cluster_all = c("1" = "#ee204d", "2" = "#17806d", "3" = "#b2ec5d", "4" = "#cda4de", "5" = "#1974d2")),
                                                annotation_name_side = "right", annotation_name_gp = gpar(fontsize = 12),
                                                annotation_legend_param = list(list(title_position = "topcenter", title = c(paste("clust_all")))),
                                                show_annotation_name = FALSE,
                                                gap = unit(c(2), "mm"),
                                                show_legend = TRUE)
}

# do separate heatmap for normals and each epithelial patient
all.equal(rownames(mats_epith_cna[[1]]), mappings.cnas$hgnc_symbol)

splits_chroms <- as.factor(mappings.cnas$chromosome_name)
splits_chroms <- factor(splits_chroms, levels(splits_chroms)[c(1, 12, 18:24, 2:11, 13:17)])

cols_chroms <- rep(c("black", "gray"), 12)
names(cols_chroms) <- c(1:24)

ha_rows_chroms <- HeatmapAnnotation(df = data.frame(chromosome = mappings.cnas$chromosome_name),
                                    annotation_legend_param = list(chromosome = list(ncol = 2, title = "chromosome", title_position = "topcenter")),
                                    which = "row", col = list("chromosome" = cols_chroms), annotation_width = unit(3, "mm"))

# heatmap without clustering annotation
ht_cnas_chroms_high <- ha_rows_chroms + 
  Heatmap(mat_normepith_cna, 
          col = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red")),
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE, show_row_names = FALSE, 
          name = "norm epith", 
          clustering_distance_columns = clustering_distance_columns, clustering_method_columns = clustering_method_columns,
          split = splits_chroms, 
          column_title = "norm epith", column_title_gp = gpar(fontsize = 10),
          #top_annotation = ha_cnas_normepith, top_annotation_height = unit(8, "mm"),
          heatmap_legend_param = list(title = paste("norm epith"), title_position = "topcenter", color_bar = "continuous",
                                      title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
          show_heatmap_legend = TRUE,
          gap = unit(1, "mm"),
          use_raster = TRUE, raster_device = "CairoPNG") + 
  Heatmap(mats_epith_cna[[1]],
          col = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red")), 
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE, show_row_names = FALSE, 
          name = patients_now[1], 
          clustering_distance_columns = clustering_distance_columns, clustering_method_columns = clustering_method_columns,
          column_title = patients_now[1], column_title_gp = gpar(fontsize = 10),
          #top_annotation = ha_cnas_epith_sep[[1]], top_annotation_height = unit(8, "mm"),
          heatmap_legend_param = list(title = paste("epithelial", patients_now[1]), title_position = "topcenter", color_bar = "continuous",
                                      title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
          show_heatmap_legend = TRUE,
          gap = unit(1, "mm"),
          use_raster = TRUE, raster_device = "CairoPNG") + 
  Heatmap(mats_epith_cna[[2]],
          col = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red")),
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE, show_row_names = FALSE, 
          name = patients_now[2], 
          clustering_distance_columns = clustering_distance_columns, clustering_method_columns = clustering_method_columns,
          column_title = patients_now[2], column_title_gp = gpar(fontsize = 10),
          #top_annotation = ha_cnas_epith_sep[[2]], top_annotation_height = unit(8, "mm"),
          heatmap_legend_param = list(title = paste(" epithelial", patients_now[2]), title_position = "topcenter", color_bar = "continuous",
                                      title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
          show_heatmap_legend = TRUE,
          gap = unit(1, "mm"),
          use_raster = TRUE, raster_device = "CairoPNG") +
  Heatmap(mats_epith_cna[[3]],
          col = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red")),
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE, show_row_names = FALSE, 
          name = patients_now[3], 
          clustering_distance_columns = clustering_distance_columns, clustering_method_columns = clustering_method_columns,
          column_title = patients_now[3], column_title_gp = gpar(fontsize = 10),
          #top_annotation = ha_cnas_epith_sep[[3]], top_annotation_height = unit(8, "mm"),
          heatmap_legend_param = list(title = paste("epithelial", patients_now[3]), title_position = "topcenter", color_bar = "continuous",
                                      title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
          show_heatmap_legend = TRUE,
          gap = unit(1, "mm"),
          use_raster = TRUE, raster_device = "CairoPNG") +
  Heatmap(mats_epith_cna[[4]],
          col = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red")),
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE, show_row_names = FALSE, 
          name = patients_now[4],
          clustering_distance_columns = clustering_distance_columns, clustering_method_columns = clustering_method_columns,
          column_title = patients_now[4], column_title_gp = gpar(fontsize = 10),
          #top_annotation = ha_cnas_epith_sep[[4]], top_annotation_height = unit(8, "mm"),
          heatmap_legend_param = list(title = paste("epithelial", patients_now[4]), title_position = "topcenter", color_bar = "continuous",
                                      title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
          show_heatmap_legend = TRUE, 
          gap = unit(1, "mm"),
          use_raster = TRUE, raster_device = "CairoPNG") +
  Heatmap(mats_epith_cna[[5]],
          col = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red")),
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE, show_row_names = FALSE, 
          name = patients_now[5], 
          clustering_distance_columns = clustering_distance_columns, clustering_method_columns = clustering_method_columns,
          column_title = patients_now[5], column_title_gp = gpar(fontsize = 10),
          #top_annotation = ha_cnas_epith_sep[[5]], top_annotation_height = unit(8, "mm"),
          heatmap_legend_param = list(title = paste("epithelial", patients_now[5]), title_position = "topcenter", color_bar = "continuous",
                                      title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
          show_heatmap_legend = TRUE, 
          gap = unit(1, "mm"),
          use_raster = TRUE, raster_device = "CairoPNG") +
  Heatmap(mats_epith_cna[[6]],
          col = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red")),
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE, show_row_names = FALSE, 
          name = patients_now[6],
          clustering_distance_columns = clustering_distance_columns, clustering_method_columns = clustering_method_columns,
          column_title = patients_now[6], column_title_gp = gpar(fontsize = 10),
          #top_annotation = ha_cnas_epith_sep[[6]], top_annotation_height = unit(8, "mm"),
          heatmap_legend_param = list(title = paste("epithelial", patients_now[6]), title_position = "topcenter", color_bar = "continuous",
                                      title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
          show_heatmap_legend = TRUE, 
          gap = unit(1, "mm"),
          use_raster = TRUE, raster_device = "CairoPNG")
#pdf(here::here("plots", "fig2d.pdf"), onefile = FALSE, height = 12, width = 12)
#print(draw(ht_cnas_chroms_high, gap = unit(0.2, "cm"), heatmap_legend_side = "bottom"))
#dev.off()

# heatmap with clustering annotation
ht_cnas_wclust_chroms_high <- ha_rows_chroms + 
  Heatmap(mat_normepith_cna, 
          col = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red")),
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE, show_row_names = FALSE, 
          name = "norm epith", 
          clustering_distance_columns = clustering_distance_columns, clustering_method_columns = clustering_method_columns,
          split = splits_chroms, 
          column_title = "norm epith", column_title_gp = gpar(fontsize = 10),
          top_annotation = ha_cnas_normepith, top_annotation_height = unit(8, "mm"),
          heatmap_legend_param = list(title = paste("norm epith"), title_position = "topcenter", color_bar = "continuous",
                                      title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
          show_heatmap_legend = TRUE,
          gap = unit(1, "mm"),
          use_raster = TRUE, raster_device = "CairoPNG") + 
  Heatmap(mats_epith_cna[[1]],
          col = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red")), 
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE, show_row_names = FALSE, 
          name = patients_now[1], 
          clustering_distance_columns = clustering_distance_columns, clustering_method_columns = clustering_method_columns,
          column_title = patients_now[1], column_title_gp = gpar(fontsize = 10),
          top_annotation = ha_cnas_epith_sep[[1]], top_annotation_height = unit(8, "mm"),
          heatmap_legend_param = list(title = paste("epithelial", patients_now[1]), title_position = "topcenter", color_bar = "continuous",
                                      title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
          show_heatmap_legend = TRUE,
          gap = unit(1, "mm"),
          use_raster = TRUE, raster_device = "CairoPNG") + 
  Heatmap(mats_epith_cna[[2]],
          col = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red")),
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE, show_row_names = FALSE, 
          name = patients_now[2], 
          clustering_distance_columns = clustering_distance_columns, clustering_method_columns = clustering_method_columns,
          column_title = patients_now[2], column_title_gp = gpar(fontsize = 10),
          top_annotation = ha_cnas_epith_sep[[2]], top_annotation_height = unit(8, "mm"),
          heatmap_legend_param = list(title = paste(" epithelial", patients_now[2]), title_position = "topcenter", color_bar = "continuous",
                                      title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
          show_heatmap_legend = TRUE,
          gap = unit(1, "mm"),
          use_raster = TRUE, raster_device = "CairoPNG") +
  Heatmap(mats_epith_cna[[3]],
          col = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red")),
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE, show_row_names = FALSE, 
          name = patients_now[3], 
          clustering_distance_columns = clustering_distance_columns, clustering_method_columns = clustering_method_columns,
          column_title = patients_now[3], column_title_gp = gpar(fontsize = 10),
          top_annotation = ha_cnas_epith_sep[[3]], top_annotation_height = unit(8, "mm"),
          heatmap_legend_param = list(title = paste("epithelial", patients_now[3]), title_position = "topcenter", color_bar = "continuous",
                                      title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
          show_heatmap_legend = TRUE,
          gap = unit(1, "mm"),
          use_raster = TRUE, raster_device = "CairoPNG") +
  Heatmap(mats_epith_cna[[4]],
          col = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red")),
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE, show_row_names = FALSE, 
          name = patients_now[4],
          clustering_distance_columns = clustering_distance_columns, clustering_method_columns = clustering_method_columns,
          column_title = patients_now[4], column_title_gp = gpar(fontsize = 10),
          top_annotation = ha_cnas_epith_sep[[4]], top_annotation_height = unit(8, "mm"),
          heatmap_legend_param = list(title = paste("epithelial", patients_now[4]), title_position = "topcenter", color_bar = "continuous",
                                      title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
          show_heatmap_legend = TRUE, 
          gap = unit(1, "mm"),
          use_raster = TRUE, raster_device = "CairoPNG") +
  Heatmap(mats_epith_cna[[5]],
          col = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red")),
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE, show_row_names = FALSE, 
          name = patients_now[5], 
          clustering_distance_columns = clustering_distance_columns, clustering_method_columns = clustering_method_columns,
          column_title = patients_now[5], column_title_gp = gpar(fontsize = 10),
          top_annotation = ha_cnas_epith_sep[[5]], top_annotation_height = unit(8, "mm"),
          heatmap_legend_param = list(title = paste("epithelial", patients_now[5]), title_position = "topcenter", color_bar = "continuous",
                                      title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
          show_heatmap_legend = TRUE, 
          gap = unit(1, "mm"),
          use_raster = TRUE, raster_device = "CairoPNG") +
  Heatmap(mats_epith_cna[[6]],
          col = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red")),
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE, show_row_names = FALSE, 
          name = patients_now[6],
          clustering_distance_columns = clustering_distance_columns, clustering_method_columns = clustering_method_columns,
          column_title = patients_now[6], column_title_gp = gpar(fontsize = 10),
          top_annotation = ha_cnas_epith_sep[[6]], top_annotation_height = unit(8, "mm"),
          heatmap_legend_param = list(title = paste("epithelial", patients_now[6]), title_position = "topcenter", color_bar = "continuous",
                                      title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
          show_heatmap_legend = TRUE, 
          gap = unit(1, "mm"),
          use_raster = TRUE, raster_device = "CairoPNG")
#pdf(here::here("plots", "figS11.pdf"), onefile = FALSE, height = 12, width = 12)
#print(draw(ht_cnas_wclust_chroms_high, gap = unit(0.2, "cm"), heatmap_legend_side = "bottom"))
#dev.off()
    

order_samples_cnv <- column_order(ht_cnas_chroms_high)
