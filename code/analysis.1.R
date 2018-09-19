# this script performs the following analyses from the manuscript: 
#   cycling cells identification
#   cell types identification (via markers and clustering)
#   t-sne of all cells
#   clustering of epithelial cells
#   differential expression between clusters
#   basal PNAS signature
#   Lehman signatures
#   normal signatures
#   mammaprint signature (PS)
#   zena werb signature (MBS)
#   carlos artega signature (RTS)
#   expression of genes in the metabolsims and immunity pathways
# this script generates the following figures from the manuscript: 
#   fig. 1b (heatmap)
#   fig. 1c (barplot)
#   fig. 1e and S2 (plots per patient)
#   fig. 1d (barplot)
#   fig. 2a, 2b, 2c (t-sne)
#   fig. S9 (per patient and per cluster plots)
#   fig. 3a (t-sne)
#   fig. S6 (barplot)
#   fig. 3d (heatmap)
#   fig. 3g and S8 (per patient and per cluster plots)
#   fig. 3e (dot plot)
#   fig. 3b (heatmap)
#   fig. 3f and S7 (per patient and per cluster plots)
#   fig. S13 (per patient and per cluster plots)
#   fig. S14 (per patient and per cluster plots)
#   fig. S15 (per patient and per cluster plots)
#   fig. 4a (heatmap)
#   fig. 4b (violin plot)
#   fig. S16 (Venn diagram)
#   fig. S10 (t-sne plots)
#   fig. S4d (heatmap)

library(here)
source(here("code", "libraries.R"))
source(here("code", "funcs.R"))
source(here("code", "funcs_markers.R"))

## read in normalized data and phenotypic information
mat_norm <- read.table(here("data", "norm_data.txt"), sep = "\t", header = TRUE)
pd_norm <- read.table(here("data", "pd_norm.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
fd_norm <- read.table(here("data", "fd_norm.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sceset_final <- SingleCellExperiment(assays = list(exprs = as.matrix(mat_norm)),
                                     colData = pd_norm, 
                                     rowData = fd_norm)


## colors for plotting
epithelial_col <- brocolors("crayons")["Maroon"]
basal_epithelial_col <- brocolors("crayons")["Red"]
luminal_epithelial_col <- brocolors("crayons")["Sunset Orange"]
luminal_progenitor_col <- brocolors("crayons")["Salmon"]

stroma_col <- brocolors("crayons")["Aquamarine"]
endothelial_col <- brocolors("crayons")["Wisteria"]

PTPRC_col <- brocolors("crayons")["Inchworm"]
t_cell_col <- brocolors("crayons")["Screamin' Green"]
b_cell_col <- brocolors("crayons")["Fern"]
macrophage_col <- brocolors("crayons")["Tropical Rain Forest"]

marker_cols <- c("epithelial" = unname(epithelial_col), "basal epithelial" = unname(basal_epithelial_col), 
                 "luminal epithelial" = unname(luminal_epithelial_col), "luminal progenitor" = unname(luminal_progenitor_col),
                 "stroma" = unname(stroma_col), "endothelial" = unname(endothelial_col), 
                 "immune" = unname(PTPRC_col), "T cell" = unname(t_cell_col), "B cell" = unname(b_cell_col), 
                 "macrophage" = unname(macrophage_col))
cycling_mel_cols <- c("non-cycling" = "gainsboro",
                      "cycling" = unname(brocolors("crayons")["Mulberry"]))
depletion_cols <- c("depleted" = unname(brocolors("crayons")["White"]), "not depleted" = unname(brocolors("crayons")["Red"]))
pats_cols <- c("PT039" = unname(brocolors("crayons")["Orange Red"]), "PT058" = unname(brocolors("crayons")["Orange"]), 
               "PT081" = unname(brocolors("crayons")["Pink Flamingo"]), "PT084" = unname(brocolors("crayons")["Fern"]), 
               "PT089" = unname(brocolors("crayons")["Blue Violet"]), "PT126" = unname(brocolors("crayons")["Sky Blue"]))
tsne_cols <- c("epithelial" = unname(basal_epithelial_col), "stroma" = unname(stroma_col), "endothelial" = unname(endothelial_col),
               "Tcell" = unname(t_cell_col), "Bcell" = unname(b_cell_col), "macrophage" = unname(macrophage_col))
#tsne_cols_unknown <- c(tsne_cols, "unknown" = "black", "undecided" = "black")
anno_colors <- list("marker" = marker_cols, "cycling" = cycling_mel_cols, "immune depletion" = depletion_cols, 
                    "patient" = pats_cols, "tsne" = tsne_cols)


## cell cycle assignment as done in melanoma in Tirosch et al 2016
melanoma_cellcycle <- read.table(here("data", "melanoma_cellcycle.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
melanoma_g1s <- melanoma_cellcycle$G1S
melanoma_g1s <- match_clean_vector_genes(mat_norm, melanoma_g1s)
scores_g1s <- avg_expr_genes(mat_norm, melanoma_g1s$index)
melanoma_g2m <- melanoma_cellcycle$G2M
melanoma_g2m <- match_clean_vector_genes(mat_norm, melanoma_g2m)
scores_g2m <- avg_expr_genes(mat_norm, melanoma_g2m$index)

cycling_mel <- rep(NA, length(scores_g1s))
for (i in 1:length(cycling_mel)) {
  if (scores_g1s[i] >= (median(scores_g1s) + 2 * mad(scores_g1s)) && scores_g2m[i] < (median(scores_g2m) + 2 * mad(scores_g2m)))
    cycling_mel[i] <- "cycling"
  if (scores_g1s[i] < (median(scores_g1s) + 2 * mad(scores_g1s)) && scores_g2m[i] >= (median(scores_g2m) + 2 * mad(scores_g2m)))
    cycling_mel[i] <- "cycling"
  if (scores_g1s[i] < (median(scores_g1s) + 2 * mad(scores_g1s)) && scores_g2m[i] < (median(scores_g2m) + 2 * mad(scores_g2m)))
    cycling_mel[i] <- "non-cycling"
  if (scores_g1s[i] >= (median(scores_g1s) + 2 * mad(scores_g1s)) && scores_g2m[i] >= (median(scores_g2m) + 2 * mad(scores_g2m)))
    cycling_mel[i] <- "cycling"
}


## update colData and pd_norm
colData(sceset_final)$mel_scores_g1s <- scores_g1s
colData(sceset_final)$mel_scores_g2m <- scores_g2m
colData(sceset_final)$cycling_mel <- cycling_mel
pd_norm <- colData(sceset_final)
pd_norm <- as.data.frame(pd_norm)


## initialize partial matrices
patients_now <- c()
mats_now <- list()
pds_now <- list()
for (i in 1:length(unique(pd_norm$patient))) {
  patients_now[i] <- sort(unique(pd_norm$patient))[i]
  mats_now[[i]] <- mat_norm[, pd_norm$patient == patients_now[i]]
  pds_now[[i]] <- pd_norm[pd_norm$patient == patients_now[i],]
}
names(mats_now) <- patients_now
names(pds_now) <- patients_now


## cell type markers
all_markers <- read.table(here("data", "markers_clean.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
markers <- unique(all_markers[all_markers$gene %in% rowData(sceset_final)$hgnc_symbol, ])

# marker names
markers$type_heatmap <- markers$type
markers$type_heatmap[which(markers$type_heatmap == "luminalprogenitor")] <- "luminal progenitor"
markers$type_heatmap[which(markers$type_heatmap == "luminalepithelial")] <- "luminal epithelial"
markers$type_heatmap[which(markers$type_heatmap == "basalepithelial")] <- "basal epithelial"
markers$type_heatmap[which(markers$type_heatmap %in% c("EPCAM", "EGFR", "CDH1"))] <- "epithelial"
markers$type_heatmap[which(markers$type_heatmap == "Bcell")] <- "B cell"
markers$type_heatmap[which(markers$type_heatmap == "Tcell")] <- "T cell"

markers$type_long_heatmap <- markers$type_long
markers$type_long_heatmap[which(markers$type == "stroma")] <- "stroma"
markers$type_long_heatmap[which(markers$type == "endothelial")] <- "endothelial"


# fig 1b
colors_markers_ch <- markers$type_heatmap
for (i in c(1:length(names(anno_colors$marker)))) {
  colors_markers_ch <- replace(colors_markers_ch, colors_markers_ch == names(anno_colors$marker)[i], anno_colors$marker[i])
}

splits_ch <- as.factor(markers$type_long_heatmap)
splits_ch <- factor(splits_ch, levels(splits_ch)[c(2,4,1,3)])

colors_anno_markers_ch <- as.factor(markers$type_heatmap)
colors_anno_markers_ch <- factor(colors_anno_markers_ch, levels(colors_anno_markers_ch)[c(4,2,5,6,9,3,8,10,1,7)])
ha_rows <- HeatmapAnnotation(df = data.frame(type = colors_anno_markers_ch),
                             annotation_legend_param = list(type = list(ncol = 2, title = "cell type", title_position = "topcenter")),
                             which = "row", col = list("type" = anno_colors$marker), annotation_width = unit(3, "mm"))

cycling_now <- list()
depletion_now <- list()
ha_cols_up <- list()
ha_cols_bottom <- list()
for (i in 1:length(patients_now)) {
  cycling_now[[i]] <- pds_now[[i]][,"cycling_mel"]
  
  if (i == 1)
    ha_cols_up[[i]] <- HeatmapAnnotation(data.frame(cycling = cycling_now[[i]]), 
                                         col = list(cycling = anno_colors$cycling), 
                                         show_annotation_name = TRUE, 
                                         annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 11),
                                         annotation_legend_param = list(list(title_position = "topcenter",
                                                                             title = c("cycling status"))),
                                         gap = unit(c(1, 1), "mm"))
  if (i > 1)
    ha_cols_up[[i]] <- HeatmapAnnotation(data.frame(cycling = cycling_now[[i]]), 
                                         col = list(cycling = anno_colors$cycling), 
                                         show_legend = FALSE,
                                         gap = unit(c(1, 1), "mm"))
}

for (i in 1:length(patients_now)) {
  depletion_now[[i]] <- pds_now[[i]][,"depletion_batch"]
  
  if (i == 1)
    ha_cols_bottom[[i]] <- HeatmapAnnotation(data.frame('CD45' = depletion_now[[i]]), 
                                             col = list('CD45' = c("depleted_yes" = "gainsboro", "depleted_no" = "gray54")),
                                             show_annotation_name = TRUE,
                                             annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 11),
                                             annotation_legend_param = list(title = "CD45 status", 
                                                                            title_position = "topcenter", 
                                                                            at = c("depleted_yes", "depleted_no"),
                                                                            labels = c("CD45 depleted","CD45 unselected")),
                                             show_legend = TRUE,
                                             gap = unit(c(1), "mm"))
  
  if (i == 2)
    ha_cols_bottom[[i]] <- HeatmapAnnotation(data.frame('CD45' = depletion_now[[i]]), 
                                             col = list('CD45' = c("depleted_yes" = "gainsboro", "depleted_no" = "gray54")),
                                             show_annotation_name = FALSE, 
                                             show_legend = FALSE,
                                             annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 11),
                                             annotation_legend_param = list(title = "CD45 status", 
                                                                            title_position = "topcenter", 
                                                                            at = c("depleted_yes", "depleted_no"),
                                                                            labels = c("CD45 depleted","CD45 unselected")),
                                             gap = unit(c(1), "mm"))
  if (i > 2)
    ha_cols_bottom[[i]] <- HeatmapAnnotation(data.frame('CD45' = depletion_now[[i]]), 
                                             col = list('CD45' = c("depleted_yes" = "gainsboro", "depleted_no" = "gray54")),
                                             show_legend = FALSE,
                                             gap = unit(c(1), "mm"))
}
names(cycling_now) <- patients_now
names(depletion_now) <- patients_now
names(ha_cols_up) <- patients_now
names(ha_cols_bottom) <- patients_now

ht_list <- ha_rows + 
  Heatmap(mats_now[[1]][match(markers$gene,rownames(mats_now[[1]])),], 
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE,
          name = patients_now[1], clustering_distance_columns = "euclidean", row_names_side = "left", 
          row_names_gp = gpar(fontsize = 10, col = colors_markers_ch),
          split = splits_ch, gap = unit(1, "mm"), column_title = patients_now[1], 
          column_title_gp = gpar(fontsize = 11),
          row_title_gp = gpar(font = 11), top_annotation = ha_cols_up[[1]], top_annotation_height = unit(1, "mm"),
          heatmap_legend_param = list(title = "expression", title_position = "topcenter", color_bar = "continuous"),
          bottom_annotation = ha_cols_bottom[[1]], bottom_annotation_height = unit(1, "mm")) + 
  Heatmap(mats_now[[2]][match(markers$gene,rownames(mats_now[[1]])),],
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE,
          name = patients_now[2], clustering_distance_columns = "euclidean",
          show_row_names = FALSE, column_title = patients_now[2],
          column_title_gp = gpar(fontsize = 11),
          show_heatmap_legend = FALSE,
          row_title_gp = gpar(font = 11), top_annotation = ha_cols_up[[2]], bottom_annotation = ha_cols_bottom[[2]],
          gap = unit(1, "mm")) +
  Heatmap(mats_now[[3]][match(markers$gene,rownames(mats_now[[3]])),],
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE,
          name = patients_now[3], clustering_distance_columns = "euclidean",
          show_row_names = FALSE, column_title = patients_now[3],
          column_title_gp = gpar(fontsize = 11),
          show_heatmap_legend = FALSE,
          row_title_gp = gpar(font = 11), top_annotation = ha_cols_up[[3]], bottom_annotation = ha_cols_bottom[[3]],
          gap = unit(1, "mm")) +
  Heatmap(mats_now[[4]][match(markers$gene,rownames(mats_now[[4]])),],
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE,
          name = patients_now[4], clustering_distance_columns = "euclidean",
          show_row_names = FALSE, column_title = patients_now[4], 
          column_title_gp = gpar(fontsize = 11),
          show_heatmap_legend = FALSE,
          row_title_gp = gpar(font = 11), top_annotation = ha_cols_up[[4]], bottom_annotation = ha_cols_bottom[[4]],
          gap = unit(1, "mm")) +
  Heatmap(mats_now[[5]][match(markers$gene,rownames(mats_now[[5]])),], 
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE,
          name = patients_now[5], clustering_distance_columns = "euclidean",
          show_row_names = FALSE, column_title = patients_now[5], 
          column_title_gp = gpar(fontsize = 11), 
          show_heatmap_legend = FALSE,
          row_title_gp = gpar(font = 11), top_annotation = ha_cols_up[[5]], bottom_annotation = ha_cols_bottom[[5]],
          gap = unit(1, "mm")) + 
  Heatmap(mats_now[[6]][match(markers$gene,rownames(mats_now[[6]])),], 
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE,
          name = patients_now[6], clustering_distance_columns = "euclidean",
          row_names_gp = gpar(fontsize = 10, col = colors_markers_ch), 
          split = splits_ch, column_title = patients_now[6], 
          column_title_gp = gpar(fontsize = 11),
          row_title_gp = gpar(font = 11), top_annotation = ha_cols_up[[6]], bottom_annotation = ha_cols_bottom[[6]],
          show_heatmap_legend = FALSE,
          gap = unit(1, "mm"))
#pdf(here("plots", "fig1a.pdf"), onefile = FALSE, width = 11, height = 10)
#draw(ht_list, gap = unit(0.1, "cm"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
#dev.off()


## cell types by markers (done as described in SI of paper)
thresh <- 1
cells_markers <- lists_markers(mat_norm, thresh, markers)
epithelial_markers <- cells_markers$epithelial_cells
is_epithelial <- decide_is_epithelial(epithelial_markers)
immune_markers <- cells_markers$immune_cells
is_immune <- decide_is_immune(immune_markers)
other_markers <- cells_markers$other_cells
is_other <- decide_is_other(other_markers)

one_epithelial_marker <- expression_one_epithelial_marker(mat_norm, pd_norm, is_epithelial, epithelial_markers, "pats", 0.5)
is_epithelial[which(one_epithelial_marker$is_epithelial_extra == 1)] <- 1

is_epithelial_simple <- is_epithelial
is_epithelial_simple[which(is_epithelial == 1)] <- "epithelial"
is_immune_simple <- is_immune
is_immune_simple[which(is_immune == "immune_mix")] <- 0
is_other_simple <- is_other
is_other_simple[which(is_other == "other_mix")] <- 0

cells_types <- paste(is_epithelial_simple, is_immune_simple, is_other_simple, sep = "_")
names(cells_types) <- names(is_epithelial)
cell_types <- sapply(strsplit(cells_types, "_"), function(x){
  # none of the cell types (epithelial, immune, other)
  if (sum(x == 0) == 3) return("unknown") else 
    if (sum(x == 0) == 2) return(setdiff(x, "0")) else
      if (sum(c("epithelial", "stroma", "0") %in% x) == 3) return("epithelial") else
        return(paste(setdiff(x, "0"),collapse = "_"))})
cell_types_simple <- cell_types
cell_types_simple[which(sapply(strsplit(cell_types, "_"), length) > 1)] <- "undecided"
table(cell_types_simple)

# update colData and pd_norm
colData(sceset_final)$cell_types_markers <- cell_types_simple
pd_norm <- colData(sceset_final)


## cell types by unsupervised clustering (done as described in SI of paper)
HSMM_clustering_ct <- monocle_unsup_clust_plots(sceset_obj = sceset_final, mat_to_cluster = mat_norm, anno_colors = anno_colors, 
                                                name_in_phenodata = "cluster_allregr_disp", disp_extra = 1, save_plots = 0,
                                                path_plots = NULL, type_pats = "allpats", regress_pat = 1)
HSMM_clustering_ct$Cluster <- HSMM_clustering_ct$cluster_allregr_disp
table(HSMM_clustering_ct$Cluster)

# due to changes in Monocle's functions (reduceDimension and clusterCells), the resulting clustering is slightly different from
# the original clustering from the paper. for reproducibility, we read in the original cluster assignment
original_clustering <- readRDS(file = here("data", "original_clustering.RDS"))
HSMM_clustering_ct$Cluster <- original_clustering
table(HSMM_clustering_ct$Cluster)

# update the sceset_final object
colData(sceset_final) <- cbind(colData(sceset_final), pData(HSMM_clustering_ct)[,c(96:104)])
colData(sceset_final)$cell_types_cl_all <- colData(sceset_final)$cell_types_markers
pd_norm <- colData(sceset_final)


## assign unknown and undecided cell types from clusters
cells_to_assign <- list()
cells_to_reassign <- list()
mean_exprs <- list()
mean_reassign_exprs <- list()

# cluster1: only epithelial cells, so nothing to be done
cluster_here <- 1
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
mean_exprs[[cluster_here]] <- NULL
mean_reassign_exprs[[cluster_here]] <- NULL
cells_to_assign[[cluster_here]] <- NULL
cells_to_reassign[[cluster_here]] <- NULL


# cluster2: try to assign the unknown and undecided cells to macrophages
cluster_here <- 2
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
to_assign_cluster <- c("undecided", "unknown")
mean_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                      cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster)),
                                                      epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)

# only assign the cells whose mean immune expression is highest
cells_to_assign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster))[which(apply(mean_exprs[[cluster_here]], 1, function(x){if (x[2] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_assign[[cluster_here]]] <- "macrophage"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])

# also check whether the epithelial cells have high immune expression
mean_reassign_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                               cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "epithelial")),
                                                               epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only re-assign the epithelial cells whose mean immune expression is highest
cells_to_reassign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "epithelial"))[which(apply(mean_reassign_exprs[[cluster_here]], 1, function(x){if (x[2] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_reassign[[cluster_here]]] <- "macrophage"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])


# cluster3: mixed, so nothing to be done
cluster_here <- 3
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
mean_exprs[[cluster_here]] <- NULL
mean_reassign_exprs[[cluster_here]] <- NULL
cells_to_assign[[cluster_here]] <- NULL
cells_to_reassign[[cluster_here]] <- NULL


# cluster 4: try to assign the unknown and undecided cells to epithelial
cluster_here <- 4
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
to_assign_cluster <- c("undecided", "unknown")
mean_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                      cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster)),
                                                      epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only assign the cells whose epithelial expression is highest
cells_to_assign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster))[which(apply(mean_exprs[[cluster_here]], 1, function(x){if (x[1] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_assign[[cluster_here]]] <- "epithelial"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])

# also check whether the stroma cells have higher epithelial expression
mean_reassign_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                               cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "stroma")),
                                                               epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only re-assign the stroma cells if their epithelial expression is higher
cells_to_reassign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "stroma"))[which(apply(mean_reassign_exprs[[cluster_here]], 1, function(x){if (x[1] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_reassign[[cluster_here]]] <- "epithelial"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])


# cluster 5: mixed, so nothing to be done
cluster_here <- 5
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
mean_exprs[[cluster_here]] <- NULL
mean_reassign_exprs[[cluster_here]] <- NULL
cells_to_assign[[cluster_here]] <- NULL
cells_to_reassign[[cluster_here]] <- NULL


# cluster 6: try to assign the unknown cells to epithelial
cluster_here <- 6
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
to_assign_cluster <- c("unknown")
mean_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                      cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster)),
                                                      epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)

# only assign the cells whose epithelial expression is highest
cells_to_assign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster))[which(apply(mean_exprs[[cluster_here]], 1, function(x){if (x[1] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_assign[[cluster_here]]] <- "epithelial"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])

# also check whether the Bcell has higher epithelial expression
mean_reassign_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                               cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "Bcell")),
                                                               epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only re-assign the Bcell if its epithelial expression is higher
cells_to_reassign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "Bcell"))[which(apply(mean_reassign_exprs[[cluster_here]], 1, function(x){if (x[1] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_reassign[[cluster_here]]] <- "epithelial"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])


# cluster 7: try to assign the unknown cell to epithelial
cluster_here <- 7
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
to_assign_cluster <- c("unknown")
mean_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                      cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster)),
                                                      epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)

# only assign the cells whose epithelial expression is highest
cells_to_assign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster))[which(apply(mean_exprs[[cluster_here]], 1, function(x){if (x[1] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_assign[[cluster_here]]] <- "epithelial"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])


# cluster 8: try to assign the unknown and undecided cells to epithelial
cluster_here <- 8
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
to_assign_cluster <- c("undecided", "unknown")
mean_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                      cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster)),
                                                      epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only assign the cells whose epithelial expression is highest
cells_to_assign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster))[which(apply(mean_exprs[[cluster_here]], 1, function(x){if (x[1] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_assign[[cluster_here]]] <- "epithelial"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])

# also check whether the stroma cells have higher epithelial expression
mean_reassign_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                               cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "stroma")),
                                                               epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only re-assign the stroma cells if their epithelial expression is higher
cells_to_reassign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "stroma"))[which(apply(mean_reassign_exprs[[cluster_here]], 1, function(x){if (x[1] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_reassign[[cluster_here]]] <- "epithelial"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])

# also check whether the Bcells or the Tcell have higher epithelial expression
mean_reassign_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                               cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% c("Bcell", "Tcell"))),
                                                               epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only re-assign the Bcells or Tcell if their epithelial expression is higher
cells_to_reassign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% c("Bcell", "Tcell")))[which(apply(mean_reassign_exprs[[cluster_here]], 1, function(x){if (x[1] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_reassign[[cluster_here]]] <- "epithelial"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])


# cluster 9: try to assign the unknown and undecided cells to macrophage
cluster_here <- 9
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
to_assign_cluster <- c("undecided", "unknown")
mean_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                      cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster)),
                                                      epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)

# only assign the cells whose mean immune expression is highest
cells_to_assign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster))[which(apply(mean_exprs[[cluster_here]], 1, function(x){if (x[2] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_assign[[cluster_here]]] <- "macrophage"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])

# also check whether the epithelial cells have high immune expression
mean_reassign_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                               cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "epithelial")),
                                                               epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only re-assign the epithelial cells whose mean immune expression is highest
cells_to_reassign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "epithelial"))[which(apply(mean_reassign_exprs[[cluster_here]], 1, function(x){if (x[2] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_reassign[[cluster_here]]] <- "macrophage"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])

colData(sceset_final)$cell_types_cl_all <- pd_norm$cell_types_cl_all


## remove the unknown and undecided cells
unkund <- which(pd_norm$cell_types_cl_all %in% c("undecided", "unknown"))

# update all the matrices
sceset_ct <- sceset_final[,-unkund]
pd_ct <- colData(sceset_ct)
mat_ct <- assays(sceset_ct)$exprs
mats_ct <- list()
pds_ct <- list()
for (i in 1:length(patients_now)) {
  mats_ct[[i]] <- mat_ct[,pd_ct$patient == patients_now[i]]
  pds_ct[[i]] <- pd_ct[pd_ct$patient == patients_now[i],]
}
names(mats_ct) <- patients_now
names(pds_ct) <- patients_now


## barplot cell types
match_celltype_levels <- c("epithelial", "stroma", "endothelial", "Tcell", "Bcell", "macrophage")
tbl_pd_ct <- tbl_df(pd_ct)
tbl_pd_ct <- tbl_pd_ct %>%
  group_by(patient) %>%
  mutate(cell_types_cl_all = factor(cell_types_cl_all, levels = match_celltype_levels)) %>%
  arrange(cell_types_cl_all)

#pdf(here("plots/", "fig1c.pdf"))
ggplot() +
  geom_bar(data = tbl_pd_ct, aes(x = patient, fill = factor(cell_types_cl_all)), position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = anno_colors$tsne) +
  labs(fill = "cell type", y = "fraction of cells")
#dev.off()


## patient cycling plots
for (i in 1:length(patients_now)) {
  #pdf(here("plots", paste("cycling_patient", patients_now[i], ".pdf", sep = "")))
  percent_epith <- length(intersect_all(which(pd_ct$patient == patients_now[i]), 
                                        which(pd_ct$cell_types_cl_all == "epithelial"), 
                                        which(pd_ct$cycling_mel == "cycling")))/length(intersect_all(
                                          which(pd_ct$patient == patients_now[i]), 
                                          which(pd_ct$cell_types_cl_all == "epithelial")))*100
  percent_all <- length(intersect_all(which(pd_ct$patient == patients_now[i]), 
                                      which(pd_ct$cycling_mel == "cycling")))/length(which(pd_ct$patient == patients_now[i]))*100
  print(ggplot(as.data.frame(pds_ct[[i]]), aes(x = mel_scores_g1s, y = mel_scores_g2m)) +
          geom_rect(ggplot2::aes(xmin = median(pd_ct$mel_scores_g1s) + 2 * mad(pd_ct$mel_scores_g1s),
                        xmax = Inf,
                        ymin = -Inf,
                        ymax = Inf),
                    fill = "gainsboro", alpha = 0.05) +
          geom_rect(aes(ymin = median(pd_ct$mel_scores_g2m) + 2 * mad(pd_ct$mel_scores_g2m),
                        ymax = Inf,
                        xmin = -Inf,
                        xmax = Inf),
                    fill = "gainsboro", alpha = 0.05) +
          geom_point(aes(col = factor(cell_types_cl_all, levels = names(anno_colors$tsne)), 
                         shape = factor(cycling_mel)), size  = 5) +
          xlim(-0.15, 2) +
          ylim(-0.15, 2.8) +
          labs(col = "cell type", shape = "cycling", x = "G1S score", y = "G2M score", 
               title = paste("patient ", patients_now[i], " (", round(percent_all), "% cycling cells)", sep = "")) + 
          scale_color_manual(values = anno_colors$tsne))
  #dev.off()
}


## barplot cell cycle phase for all patients
tbl_pd_cycle <- tbl_pd_ct %>%
  group_by(patient) %>%
  mutate(cycling_mel = factor(cycling_mel, levels = c("cycling", "non-cycling"))) %>%
  arrange(cycling_mel)

#pdf(here("plots", "fig1d.pdf"))
ggplot() +
  geom_bar(data = tbl_pd_cycle, aes(x = patient, fill = factor(cycling_mel)), position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = anno_colors$cycling) +
  labs(fill = "cycling status", y = "fraction of cells")
#dev.off()


## tsne on cell types
# all cells
to_plot_ct <- unique(pd_ct$cell_types_cl_all)
mat_short_ct <- mat_ct[, which(pd_ct$cell_types_cl_all %in% to_plot_ct)]
pd_short_ct <- pd_ct[which(pd_ct$cell_types_cl_all %in% to_plot_ct), ]
tsne_short_ct <- Rtsne(t(mat_short_ct), perplexity = 30)
colnames(tsne_short_ct$Y) <- c("col1", "col2")
tsne_short_ct$Y <- as.data.frame(tsne_short_ct$Y)
tsne_short_ct$Y$cell_types_cl_all <- pd_short_ct$cell_types_cl_all
tsne_short_ct$Y$cell_types_markers <- pd_short_ct$cell_types_markers
tsne_short_ct$Y$patient <- pd_short_ct$patient

#pdf(here("plots", "fig2a.pdf"))
ggplot(tsne_short_ct$Y, aes(x = col1, y = col2, color = factor(cell_types_cl_all, levels = names(anno_colors$tsne)), 
                            shape = patient)) + 
  geom_point(size = 4) + 
  scale_color_manual(values = anno_colors$tsne) +
  labs(col = "patient", x = "tSNE dimension 1", y = "tSNE dimension 2", shape = "patient")
#dev.off()

# epithelial cells
to_plot_ct <- c("epithelial")
mat_short_ct <- mat_ct[, which(pd_ct$cell_types_cl_all %in% to_plot_ct)]
pd_short_ct <- pd_ct[which(pd_ct$cell_types_cl_all %in% to_plot_ct), ]
tsne_short_ct <- Rtsne(t(mat_short_ct), perplexity = 30)
colnames(tsne_short_ct$Y) <- c("col1", "col2")
tsne_short_ct$Y <- as.data.frame(tsne_short_ct$Y)
tsne_short_ct$Y$cell_types_cl_all <- pd_short_ct$cell_types_cl_all
tsne_short_ct$Y$cell_types_markers <- pd_short_ct$cell_types_markers
tsne_short_ct$Y$patient <- pd_short_ct$patient

#pdf(here("plots", "fig2b.pdf"))
ggplot(tsne_short_ct$Y, aes(x = col1, y = col2, color = factor(patient, levels = names(anno_colors$patient)), 
                            shape = cell_types_cl_all)) + 
  geom_point(size = 4) + 
  scale_color_manual(values = anno_colors$patient) +
  labs(col = "patient", x = "tSNE dimension 1", y = "tSNE dimension 2", shape = "cell type")
#dev.off()

# normal cells
to_plot_ct <- c("Bcell", "macrophage", "Tcell", "stroma", "endothelial")
mat_short_ct <- mat_ct[, which(pd_ct$cell_types_cl_all %in% to_plot_ct)]
pd_short_ct <- pd_ct[which(pd_ct$cell_types_cl_all %in% to_plot_ct), ]
tsne_short_ct <- Rtsne(t(mat_short_ct), perplexity = 30)
colnames(tsne_short_ct$Y) <- c("col1", "col2")
tsne_short_ct$Y <- as.data.frame(tsne_short_ct$Y)
tsne_short_ct$Y$cell_types_cl_all <- pd_short_ct$cell_types_cl_all
tsne_short_ct$Y$cell_types_markers <- pd_short_ct$cell_types_markers
tsne_short_ct$Y$patient <- pd_short_ct$patient

#pdf(here("plots", "fig2c.pdf"))
ggplot(tsne_short_ct$Y, aes(x = col1, y = col2, color = factor(cell_types_cl_all, levels = names(anno_colors$tsne)), 
                            shape = patient)) + 
  geom_point(size = 4) + 
  labs(col = "cell type", x = "tSNE dimension 1", y = "tSNE dimension 2") + 
  scale_color_manual(values = anno_colors$tsne)
#dev.off()


## clustering of epithelial cells
HSMM_allepith_clustering <- monocle_unsup_clust_plots(sceset_obj = sceset_ct[,which(colData(sceset_ct)$cell_types_cl_all == "epithelial")], 
                                                      mat_to_cluster = mat_ct[,which(colData(sceset_ct)$cell_types_cl_all == "epithelial")], 
                                                      anno_colors = anno_colors, name_in_phenodata = "cluster_allepith_regr_disp", 
                                                      disp_extra = 1, save_plots = 0, path_plots = NULL, 
                                                      type_pats = "allpats", regress_pat = 1, use_known_colors = 1, use_only_known_celltypes = 1)
table(HSMM_allepith_clustering$Cluster)


# due to changes in Monocle's functions (reduceDimension and clusterCells), the resulting clustering of epithelial cells
# is slightly different from the original clustering from the paper. for reproducibility, we read in the original 
# clustering of epithelial cells
original_clustering_epithelial <- readRDS(file = here("data", "original_clustering_epithelial.RDS"))
table(original_clustering_epithelial)

HSMM_allepith_clustering$Cluster <- original_clustering_epithelial
clustering_allepith <- HSMM_allepith_clustering$Cluster

#pdf(here("plots", "fig3a.pdf"))
plot_cell_clusters(HSMM_allepith_clustering, 1, 2, color = "Cluster", cell_size = 2) + 
  scale_color_manual(values = c("1" = "#ee204d", "2" = "#17806d", "3" = "#b2ec5d", "4" = "#cda4de", "5" = "#1974d2"))
#dev.off()

clusterings_sep_allepith <- list()
for (i in patients_now) {
  clusterings_sep_allepith[[i]] <- clustering_allepith[which(HSMM_allepith_clustering$patient == i)]
  names(clusterings_sep_allepith[[i]]) <- colnames(HSMM_allepith_clustering)[which(HSMM_allepith_clustering$patient == i)]
}


## barplot cycling/non-cycling per Cluster
tbl_pd_cluster <- tbl_df(pData(HSMM_allepith_clustering))
tbl_pd_cluster <- tbl_pd_cluster %>%
  group_by(Cluster) %>%
  mutate(cycling_mel = factor(cycling_mel, levels = c("cycling", "non-cycling"))) %>%
  arrange(cycling_mel)

#pdf(here("plots", "figS6.pdf"))
ggplot() +
  geom_bar(data = tbl_pd_cluster, aes(x = Cluster, fill = factor(cycling_mel)), position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = anno_colors$cycling) +
  labs(fill = "cycling status", y = "fraction of cells")
#dev.off()


## differential expression between clusters, each cluster vs. all the others
HSMM_for_DE <- HSMM_allepith_clustering
diff_test_res <- list()

HSMM_for_DE$allvs1 <- clustering_allepith
HSMM_for_DE$allvs1 <- as.numeric(HSMM_for_DE$allvs1)
HSMM_for_DE$allvs1[which(HSMM_for_DE$allvs1 != 1)] <- 2
diff_test_res$allvs1 <- differentialGeneTest(HSMM_for_DE, fullModelFormulaStr = "~allvs1", cores = 3)
diff_test_res$allvs1 <- diff_test_res$allvs1[order(diff_test_res$allvs1$qval),]
diff_test_res$allvs1 <- diff_test_res$allvs1[which(diff_test_res$allvs1$qval <= 0.1),]
head(diff_test_res$allvs1[,1:5], n = 10)
#write.table(diff_test_res$allvs1, here("tables", "diffexp_genes_allvs1.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

HSMM_for_DE$allvs2 <- clustering_allepith
HSMM_for_DE$allvs2 <- as.numeric(HSMM_for_DE$allvs2)
HSMM_for_DE$allvs2[which(HSMM_for_DE$allvs2 != 2)] <- 3
diff_test_res$allvs2 <- differentialGeneTest(HSMM_for_DE, fullModelFormulaStr = "~allvs2", cores = 3)
diff_test_res$allvs2 <- diff_test_res$allvs2[order(diff_test_res$allvs2$qval),]
diff_test_res$allvs2 <- diff_test_res$allvs2[which(diff_test_res$allvs2$qval <= 0.1),]
head(diff_test_res$allvs2[,1:5], n = 10)
#write.table(diff_test_res$allvs2, here("tables", "diffexp_genes_allvs2.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

HSMM_for_DE$allvs3 <- clustering_allepith
HSMM_for_DE$allvs3 <- as.numeric(HSMM_for_DE$allvs3)
HSMM_for_DE$allvs3[which(HSMM_for_DE$allvs3 != 3)] <- 4
diff_test_res$allvs3 <- differentialGeneTest(HSMM_for_DE, fullModelFormulaStr = "~allvs3", cores = 3)
diff_test_res$allvs3 <- diff_test_res$allvs3[order(diff_test_res$allvs3$qval),]
diff_test_res$allvs3 <- diff_test_res$allvs3[which(diff_test_res$allvs3$qval <= 0.1),]
head(diff_test_res$allvs3[,1:5], n = 10)
#write.table(diff_test_res$allvs3, here("tables", "diffexp_genes_allvs3.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

HSMM_for_DE$allvs4 <- clustering_allepith
HSMM_for_DE$allvs4 <- as.numeric(HSMM_for_DE$allvs4)
HSMM_for_DE$allvs4[which(HSMM_for_DE$allvs4 != 4)] <- 5
diff_test_res$allvs4 <- differentialGeneTest(HSMM_for_DE, fullModelFormulaStr = "~allvs4", cores = 3)
diff_test_res$allvs4 <- diff_test_res$allvs4[order(diff_test_res$allvs4$qval),]
diff_test_res$allvs4 <- diff_test_res$allvs4[which(diff_test_res$allvs4$qval <= 0.1),]
head(diff_test_res$allvs4[,1:5], n = 10)
#write.table(diff_test_res$allvs4, here("tables", "diffexp_genes_allvs4.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

HSMM_for_DE$allvs5 <- clustering_allepith
HSMM_for_DE$allvs5 <- as.numeric(HSMM_for_DE$allvs5)
HSMM_for_DE$allvs5[which(HSMM_for_DE$allvs5 != 5)] <- 6
diff_test_res$allvs5 <- differentialGeneTest(HSMM_for_DE, fullModelFormulaStr = "~allvs5", cores = 3)
diff_test_res$allvs5 <- diff_test_res$allvs5[order(diff_test_res$allvs5$qval),]
diff_test_res$allvs5 <- diff_test_res$allvs5[which(diff_test_res$allvs5$qval <= 0.1),]
head(diff_test_res$allvs5[,1:5], n = 10)
#write.table(diff_test_res$allvs5, here("tables", "diffexp_genes_allvs5.txt"), sep = "\t", quote = FALSE, row.names = FALSE)


## basal vs. nonbasal PNAS signature
basal_PNAS_all <- read.table(here("data", "genes_for_basal_vs_non_basal_tnbc_PNAS.txt"), header = TRUE, sep = "\t")
basal_PNAS_long <- basal_PNAS_all$Basal.epithelial.cell.enriched.cluster
basal_PNAS <- intersect(basal_PNAS_long, rownames(mat_ct))
basal_PNAS_avg_exprs <- apply(mat_ct[match(basal_PNAS, rownames(mat_ct)),], 2, mean)
all.equal(names(basal_PNAS_avg_exprs), colnames(mat_ct))
basal_PNAS_avg_exprs <- basal_PNAS_avg_exprs[which(pd_ct$cell_types_cl_all == "epithelial")]

# per patient and per cluster PAM50 signature plots
all.equal(colnames(HSMM_allepith_clustering), names(basal_PNAS_avg_exprs))
pData(HSMM_allepith_clustering)$basal_PNAS_avg_exprs <- basal_PNAS_avg_exprs

#pdf(here("plots", "figS9b.pdf"), width = 10)
plot_cell_clusters(HSMM_allepith_clustering, 1, 2, color = "basal_PNAS_avg_exprs", cell_size = 2) + facet_wrap(~patient) +
  scale_color_continuous(low = "yellow", high = "blue")
#dev.off()

#pdf(here("plots", "figS9a.pdf"), width = 10)
plot_cell_clusters(HSMM_allepith_clustering, 1, 2, color = "basal_PNAS_avg_exprs", cell_size = 2) + facet_wrap(~Cluster) +
  scale_color_continuous(low = "yellow", high = "blue")
#dev.off()


## Lehman signature
lehman_long <- read.table(here("data", "Lehman_signature.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
for (i in 0:5) {
  
  gene <- "gene"
  regulation <- "regulation"
  no_samples <- "no_samples"
  signature <- "signature"
  
  if (i == 0) {
    lehman <- lehman_long[, 1:4]
    lehman <- lehman[-which(lehman$signature == ""),]
  }
  
  if (i > 0) {
    gene <- paste("gene", i, sep = ".")
    regulation <- paste("regulation", i, sep = ".")
    no_samples <- paste("no_samples", i, sep = ".")
    signature <- paste("signature", i, sep = ".")
    
    mat_to_bind <- lehman_long[, c(gene, regulation, no_samples, signature)]
    colnames(mat_to_bind) <- c("gene", "regulation", "no_samples", "signature")
    if (length(which(is.na(mat_to_bind$no_samples))) > 0 )
      mat_to_bind <- mat_to_bind[-which(mat_to_bind$signature == ""),]
    lehman <- rbind(lehman, mat_to_bind)
  }
}

lehman <- tbl_df(lehman) %>%
  group_by(signature)
lehman <- lehman[which(!is.na(match(lehman$gene, rownames(mat_ct)))),]

lehman_signatures <- unique(lehman$signature)

lehman_avg_exps <- apply(mat_ct, 2, function(x){
  
  mns <- matrix(NA, nrow = length(lehman_signatures), ncol = 2)
  rownames(mns) <- lehman_signatures
  for (s in 1:length(lehman_signatures)) {
    sign <- lehman_signatures[s] # current signature
    lehman_here <- lehman %>%
      dplyr::filter(signature == sign)
    lehman_here_up <- lehman_here %>%
      dplyr::filter(regulation == "UP")
    lehman_here_down <- lehman_here %>%
      dplyr::filter(regulation == "DOWN")
    
    # indices of genes in the expression matrix
    idx_genes_up <- match(lehman_here_up$gene, rownames(mat_ct)) 
    idx_genes_down <- match(lehman_here_down$gene, rownames(mat_ct))
    
    mns[s,] <- c(mean(x[idx_genes_up]), mean(x[idx_genes_down]))
  }
  return(mns)
})
all.equal(colnames(lehman_avg_exps), rownames(pd_ct))
lehman_avg_exprs_epithelial <- lehman_avg_exps[,which(pd_ct$cell_types_cl_all == "epithelial")]

lehman_avg_ups <- lehman_avg_exps[c(1:6), ]
rownames(lehman_avg_ups) <- lehman_signatures
all.equal(colnames(lehman_avg_ups), rownames(pd_ct))
lehman_avg_ups_epithelial <- lehman_avg_ups[,which(pd_ct$cell_types_cl_all == "epithelial")]

lehman_avg_downs <- lehman_avg_exps[c(7:12),]
rownames(lehman_avg_downs) <- lehman_signatures
all.equal(colnames(lehman_avg_downs), rownames(pd_ct))
lehman_avg_downs_epithelial <- lehman_avg_downs[,which(pd_ct$cell_types_cl_all == "epithelial")]

lehman_avg_both <- lehman_avg_ups - lehman_avg_downs
all.equal(colnames(lehman_avg_both), rownames(pd_ct))
lehman_avg_both_epithelial <- lehman_avg_both[,which(pd_ct$cell_types_cl_all == "epithelial")]

assignments_lehman_both <- apply(lehman_avg_both, 2, function(x){rownames(lehman_avg_both)[which.max(x)]})
assignments_lehman_both_epithelial <- assignments_lehman_both[which(pd_ct$cell_types_cl_all == "epithelial")]

# update lehman signatures by removing the immunomodulatory and mesenchymal_stem_like signatures
lehman_avg_both_epithelial_new <- lehman_avg_both_epithelial[-which(rownames(lehman_avg_both_epithelial) %in% c("immunomodulatory", "mesenchymal_stem_like")),]
assignments_lehman_both_epithelial_new <- apply(lehman_avg_both_epithelial_new, 2, function(x){rownames(lehman_avg_both_epithelial_new)[which.max(x)]})

# Heatmap on lehman expression per patient
pd_ct_epith <- pd_ct[which(pd_ct$cell_types_cl_all == "epithelial"),]
lehmans_epith_pat_both <- list()
lehmans_epith_pat_ups <- list()
pds_epith_ct <- list()
for (i in 1:length(patients_now)) {
  
  lehmans_epith_pat_both[[i]] <- lehman_avg_both_epithelial[,which(pd_ct_epith$patient == patients_now[i])]
  lehmans_epith_pat_ups[[i]] <- lehman_avg_ups_epithelial[,which(pd_ct_epith$patient == patients_now[i])]
  pds_epith_ct[[i]] <- pds_ct[[i]][which(pds_ct[[i]]$cell_types_cl_all == "epithelial"),]
  
}
names(lehmans_epith_pat_both) <- patients_now
names(lehmans_epith_pat_ups) <- patients_now
names(pds_epith_ct) <- patients_now

# lehmann expression per patient
lehmans_epith_pat_both_new <- list()
for (i in 1:length(patients_now)) {
  lehmans_epith_pat_both_new[[i]] <- lehman_avg_both_epithelial_new[,which(pd_ct_epith$patient == patients_now[i])]
}
names(lehmans_epith_pat_both_new) <- patients_now

# annotations for epithelial cells separate per patiet
ha_lehman_epith_pat <- list()
for (i in 1:length(patients_now)) {
  
  if (i == 1)
    ha_lehman_epith_pat[[i]] <- HeatmapAnnotation(data.frame(cluster_all = clusterings_sep_allepith[[i]]), 
                                                  col = list(cluster_all = c("1" = "#ee204d", "2" = "#17806d", "3" = "#b2ec5d", "4" = "#cda4de", "5" = "#1974d2")),
                                                  annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 12),
                                                  annotation_legend_param = list(list(title_position = "topcenter", title = "cluster")),
                                                  show_annotation_name = FALSE,
                                                  gap = unit(c(2), "mm"),
                                                  show_legend = FALSE)
  
  if (i > 1 && i != 5 )
    ha_lehman_epith_pat[[i]] <- HeatmapAnnotation(data.frame(cluster_all = clusterings_sep_allepith[[i]]), 
                                                  col = list(cluster_all = c("1" = "#ee204d", "2" = "#17806d", "3" = "#b2ec5d", "4" = "#cda4de", "5" = "#1974d2")),
                                                  annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 12),
                                                  annotation_legend_param = list(list(title_position = "topcenter", title = "cluster")),
                                                  show_annotation_name = FALSE,
                                                  gap = unit(c(2), "mm"),
                                                  show_legend = FALSE)
  
  if (i == 5)
    ha_lehman_epith_pat[[i]] <- HeatmapAnnotation(data.frame(cluster_all = clusterings_sep_allepith[[i]]), 
                                                  col = list(cluster_all = c("1" = "#ee204d", "2" = "#17806d", "3" = "#b2ec5d", "4" = "#cda4de", "5" = "#1974d2")),
                                                  annotation_name_side = "right", annotation_name_gp = gpar(fontsize = 12),
                                                  annotation_legend_param = list(list(title_position = "topcenter",title = "cluster")),
                                                  show_annotation_name = FALSE,
                                                  gap = unit(c(2), "mm"),
                                                  show_legend = TRUE)
}
all.equal(names(lehmans_epith_pat_both), patients_now)


# add basal signature
lehmans_epith_pat_both_wbasal_new <- lehmans_epith_pat_both_new
for (i in 1:length(patients_now)) {
  lehmans_epith_pat_both_wbasal_new[[i]] <- rbind(lehmans_epith_pat_both_new[[i]], pData(HSMM_allepith_clustering)$basal_PNAS_avg_exprs[which(HSMM_allepith_clustering$patient == patients_now[i])])
  rownames(lehmans_epith_pat_both_wbasal_new[[i]])[5] <- "intrinsic_basal"
}

# new heatmap with basal
ht_sep_lehmans_both_wbasal_new <-
  Heatmap(lehmans_epith_pat_both_wbasal_new[[1]],
          col = colorRamp2(c(-0.7, 0, 1), c("blue","white", "red")),
          cluster_rows = FALSE,
          show_column_names = FALSE,
          column_title = patients_now[1],
          column_title_gp = gpar(fontsize = 12),
          top_annotation = ha_lehman_epith_pat[[1]],
          name = patients_now[1], 
          show_row_names = FALSE,
          top_annotation_height = unit(c(1), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9))) +
  Heatmap(lehmans_epith_pat_both_wbasal_new[[2]],
          col = colorRamp2(c(-0.7, 0, 1), c("blue","white", "red")),
          cluster_rows = FALSE,
          show_column_names = FALSE,
          column_title = patients_now[2],
          column_title_gp = gpar(fontsize = 12),
          top_annotation = ha_lehman_epith_pat[[2]],
          name = patients_now[2], 
          show_row_names = FALSE,
          top_annotation_height = unit(c(1), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9))) +
  Heatmap(lehmans_epith_pat_both_wbasal_new[[3]],
          col = colorRamp2(c(-0.7, 0, 1), c("blue","white", "red")),
          cluster_rows = FALSE,
          show_column_names = FALSE,
          column_title = patients_now[3],
          column_title_gp = gpar(fontsize = 12),
          top_annotation = ha_lehman_epith_pat[[3]],
          name = patients_now[3], 
          show_row_names = FALSE,
          top_annotation_height = unit(c(1), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9))) +
  Heatmap(lehmans_epith_pat_both_wbasal_new[[4]],
          col = colorRamp2(c(-0.7, 0, 1), c("blue","white", "red")),
          cluster_rows = FALSE,
          show_column_names = FALSE,
          column_title = patients_now[4],
          column_title_gp = gpar(fontsize = 12),
          top_annotation = ha_lehman_epith_pat[[4]],
          name = patients_now[4], 
          show_row_names = FALSE,
          top_annotation_height = unit(c(1), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9))) +
  Heatmap(lehmans_epith_pat_both_wbasal_new[[5]],
          col = colorRamp2(c(-0.7, 0, 1), c("blue","white", "red")),
          cluster_rows = FALSE,
          show_column_names = FALSE,
          column_title = patients_now[5],
          column_title_gp = gpar(fontsize = 12),
          top_annotation = ha_lehman_epith_pat[[5]],
          name = patients_now[5], 
          show_row_names = FALSE,
          top_annotation_height = unit(c(1), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9))) +
  Heatmap(lehmans_epith_pat_both_wbasal_new[[6]],
          col = colorRamp2(c(-0.7, 0, 1), c("blue","white", "red")),
          cluster_rows = FALSE,
          row_names_side = "right",
          column_title = patients_now[6],
          column_title_gp = gpar(fontsize = 12),
          top_annotation = ha_lehman_epith_pat[[6]],
          name = patients_now[6], 
          show_column_names = FALSE,
          top_annotation_height = unit(c(1), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9)))
#pdf(here("plots", "fig3d.pdf"), onefile = FALSE, width = 20)
print(draw(ht_sep_lehmans_both_wbasal_new, annotation_legend_side = "bottom"))
#dev.off()


# per patient and per cluster Lehman plots
all.equal(colnames(HSMM_allepith_clustering), names(assignments_lehman_both_epithelial_new))
pData(HSMM_allepith_clustering)$assignments_lehman_both_new <- assignments_lehman_both_epithelial_new

#pdf(here("plots", "fig3g.pdf"), width = 10)
plot_cell_clusters(HSMM_allepith_clustering, 1, 2, color = "assignments_lehman_both_new", cell_size = 2) + facet_wrap(~patient)
#dev.off()

#pdf(here("plots", "figS8.pdf"), width = 10)
plot_cell_clusters(HSMM_allepith_clustering, 1, 2, color = "assignments_lehman_both_new", cell_size = 2) + facet_wrap(~Cluster)
#dev.off()


# dot plots new lehman signatures
clust_avg_lehman_both_new <- matrix(NA, nrow = length(unique(HSMM_allepith_clustering$Cluster)), ncol = nrow(lehman_avg_both_epithelial_new))
rownames(clust_avg_lehman_both_new) <- paste("clust", c(1:length(unique(HSMM_allepith_clustering$Cluster))), sep = "")
colnames(clust_avg_lehman_both_new) <- rownames(lehman_avg_both_epithelial_new)
for (c in 1:length(unique(HSMM_allepith_clustering$Cluster))) {
  clust_avg_lehman_both_new[c,] <- apply(lehman_avg_both_epithelial_new[,which(HSMM_allepith_clustering$Cluster == c)], 1, mean)
}

clust_avg_lehman_both_new <- as.data.frame(clust_avg_lehman_both_new)
clust_avg_lehman_both_new$Cluster <- rownames(clust_avg_lehman_both_new)
clust_avg_lehman_melt_new <- melt(clust_avg_lehman_both_new, "Cluster")

pdf(here("plots", "fig3e.pdf"), width = 7)
ggplot(clust_avg_lehman_melt_new, aes(Cluster, value, fill = factor(variable), color = factor(variable), 
                                      shape = factor(variable))) + 
  geom_point(size = 3, stroke = 1) +
  scale_shape_discrete(solid = T) + 
  #guides(colour = guide_legend(override.aes = list(size=3))) + 
  ylab("average expression of signature in cluster") +
  xlab("cluster") +
  ylim(c(-0.35, 0.5))
dev.off()


## normal signatures
ml_signature_long <- read.table(here("data", "ML_signature.txt"), sep = "\t", header = TRUE)
if (length(which(ml_signature_long$Symbol == "")) > 0)
  ml_signature_long <- ml_signature_long[-which(ml_signature_long$Symbol == ""),]
ml_signature_long <- ml_signature_long[order(ml_signature_long$Symbol, -abs(ml_signature_long$Average.log.fold.change) ), ]
ml_signature_long <- ml_signature_long[ !duplicated(ml_signature_long$Symbol), ]
ml_signature <- ml_signature_long[which(!is.na(match(ml_signature_long$Symbol, rownames(mat_ct)))), ]
ml_up <- ml_signature[which(ml_signature$Average.log.fold.change > 0), ]
ml_down <- ml_signature[which(ml_signature$Average.log.fold.change < 0), ]
idx_ml_up <- match(ml_up$Symbol, rownames(mat_ct))
idx_ml_down <- match(ml_down$Symbol, rownames(mat_ct))

basal_signature_long <- read.table(here("data", "basal_signature.txt"), sep = "\t", header = TRUE)
if (length(which(basal_signature_long$Symbol == "")) > 0)
  basal_signature_long <- basal_signature_long[-which(basal_signature_long$Symbol == ""),]
basal_signature_long <- basal_signature_long[order(basal_signature_long$Symbol, -abs(basal_signature_long$Average.log.fold.change) ), ]
basal_signature_long <- basal_signature_long[ !duplicated(basal_signature_long$Symbol), ]
basal_signature <- basal_signature_long[which(!is.na(match(basal_signature_long$Symbol, rownames(mat_ct)))), ]
basal_up <- basal_signature[which(basal_signature$Average.log.fold.change > 0), ]
basal_down <- basal_signature[which(basal_signature$Average.log.fold.change < 0), ]
idx_basal_up <- match(basal_up$Symbol, rownames(mat_ct))
idx_basal_down <- match(basal_down$Symbol, rownames(mat_ct))

lp_signature_long <- read.table(here("data", "lp_signature.txt"), sep = "\t", header = TRUE)
if (length(which(lp_signature_long$Symbol == "")) > 0)
  lp_signature_long <- lp_signature_long[-which(lp_signature_long$Symbol == ""),]
lp_signature_long <- lp_signature_long[order(lp_signature_long$Symbol, -abs(lp_signature_long$Average.log.fold.change) ), ]
lp_signature_long <- lp_signature_long[ !duplicated(lp_signature_long$Symbol), ]
lp_signature <- lp_signature_long[which(!is.na(match(lp_signature_long$Symbol, rownames(mat_ct)))), ]
lp_up <- lp_signature[which(lp_signature$Average.log.fold.change > 0), ]
lp_down <- lp_signature[which(lp_signature$Average.log.fold.change < 0), ]
idx_lp_up <- match(lp_up$Symbol, rownames(mat_ct))
idx_lp_down <- match(lp_down$Symbol, rownames(mat_ct))

normsig_avg_exprs <- apply(mat_ct, 2, function(x){
  
  avg_ml_up <- mean(x[idx_ml_up])
  avg_ml_down <- mean(x[idx_ml_down])
  avg_ml_both <- avg_ml_up - avg_ml_down
  
  avg_basal_up <- mean(x[idx_basal_up])
  avg_basal_down <- mean(x[idx_basal_down])
  avg_basal_both <- avg_basal_up - avg_basal_down
  
  avg_lp_up <- mean(x[idx_lp_up])
  avg_lp_down <- mean(x[idx_lp_down])
  avg_lp_both <- avg_lp_up - avg_lp_down
  
  return(c(avg_ml_up, avg_basal_up, avg_lp_up, avg_ml_both, avg_basal_both, avg_lp_both))
})
rownames(normsig_avg_exprs) <- c("avg_ml_up", "avg_basal_up", "avg_lp_up", "avg_ml_both", "avg_basal_both", "avg_lp_both")
all.equal(colnames(normsig_avg_exprs), rownames(pd_ct))
normsig_avg_exprs_epithelial <- normsig_avg_exprs[,which(pd_ct$cell_types_cl_all == "epithelial")]

normsig_avg_ups <- normsig_avg_exprs[c(1:3), ]
all.equal(colnames(normsig_avg_ups), rownames(pd_ct))
normsig_avg_ups_epithelial <- normsig_avg_ups[,which(pd_ct$cell_types_cl_all == "epithelial")]

normsig_avg_both <- normsig_avg_exprs[c(4:6),]
all.equal(colnames(normsig_avg_both), rownames(pd_ct))
normsig_avg_both_epithelial <- normsig_avg_both[,which(pd_ct$cell_types_cl_all == "epithelial")]

assignments_normsig_ups <- apply(normsig_avg_ups, 2, function(x){rownames(normsig_avg_ups)[which.max(x)]})
assignments_normsig_ups_epithelial <- assignments_normsig_ups[which(pd_ct$cell_types_cl_all == "epithelial")]
assignments_normsig_both <- apply(normsig_avg_both, 2, function(x){rownames(normsig_avg_both)[which.max(x)]})
assignments_normsig_both_epithelial <- assignments_normsig_both[which(pd_ct$cell_types_cl_all == "epithelial")]

# heatmaps on normal signatures per patient
pd_ct_epith <- pd_ct[which(pd_ct$cell_types_cl_all == "epithelial"),]
normsig_epith_pat_both <- list()
normsig_epith_pat_ups <- list()
pds_epith_ct <- list()
for (i in 1:length(patients_now)) {
  normsig_epith_pat_both[[i]] <- normsig_avg_both_epithelial[,which(pd_ct_epith$patient == patients_now[i])]
  normsig_epith_pat_ups[[i]] <- normsig_avg_ups_epithelial[,which(pd_ct_epith$patient == patients_now[i])]
  pds_epith_ct[[i]] <- pds_ct[[i]][which(pds_ct[[i]]$cell_types_cl_all == "epithelial"),]
}
names(normsig_epith_pat_both) <- patients_now
names(normsig_epith_pat_ups) <- patients_now
names(pds_epith_ct) <- patients_now

ht_sep_normsig_both <-
  Heatmap(normsig_epith_pat_both[[1]],
          col = colorRamp2(c(-0.7, -0.2, 0.7), c("blue","white", "red")),
          cluster_rows = FALSE,
          show_column_names = FALSE,
          column_title = patients_now[1],
          top_annotation = ha_lehman_epith_pat[[1]],
          column_title_gp = gpar(fontsize = 12),
          show_row_names = FALSE,
          name = patients_now[1], 
          top_annotation_height = unit(c(2), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9))) +
  Heatmap(normsig_epith_pat_both[[2]],
          col = colorRamp2(c(-0.7, -0.2, 0.7), c("blue","white", "red")),
          cluster_rows = FALSE,
          show_column_names = FALSE,
          column_title = patients_now[2],
          column_title_gp = gpar(fontsize = 12),
          top_annotation = ha_lehman_epith_pat[[2]],
          name = patients_now[2], 
          show_row_names = FALSE,
          top_annotation_height = unit(c(2), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9))) +
  Heatmap(normsig_epith_pat_both[[3]],
          col = colorRamp2(c(-0.7, -0.2, 0.7), c("blue","white", "red")),
          cluster_rows = FALSE,
          show_column_names = FALSE,
          column_title = patients_now[3],
          column_title_gp = gpar(fontsize = 12),
          top_annotation = ha_lehman_epith_pat[[3]],
          name = patients_now[3], 
          show_row_names = FALSE,
          top_annotation_height = unit(c(2), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9))) +
  Heatmap(normsig_epith_pat_both[[4]],
          col = colorRamp2(c(-0.7, -0.2, 0.7), c("blue","white", "red")),
          cluster_rows = FALSE,
          show_column_names = FALSE,
          column_title = patients_now[4],
          column_title_gp = gpar(fontsize = 12),
          top_annotation = ha_lehman_epith_pat[[4]],
          name = patients_now[4], 
          show_row_names = FALSE,
          top_annotation_height = unit(c(2), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9))) +
  Heatmap(normsig_epith_pat_both[[5]],
          col = colorRamp2(c(-0.7, -0.2, 0.7), c("blue","white", "red")),
          cluster_rows = FALSE,
          show_column_names = FALSE,
          column_title = patients_now[5],
          column_title_gp = gpar(fontsize = 12),
          top_annotation = ha_lehman_epith_pat[[5]],
          name = patients_now[5], 
          show_row_names = FALSE,
          top_annotation_height = unit(c(2), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9))) +
  Heatmap(normsig_epith_pat_both[[6]],
          col = colorRamp2(c(-0.7, -0.2, 0.7), c("blue","white", "red")),
          cluster_rows = FALSE,
          row_names_side = "right",
          column_title = patients_now[6],
          column_title_gp = gpar(fontsize = 12),
          top_annotation = ha_lehman_epith_pat[[6]],
          name = patients_now[6], 
          show_column_names = FALSE,
          top_annotation_height = unit(c(2), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9)))
#pdf(here("plots", "fig3b.pdf"), onefile = FALSE, width = 20)
print(draw(ht_sep_normsig_both, annotation_legend_side = "bottom"))
#dev.off()


# per patient and per cluster normal signatures plots
all.equal(colnames(HSMM_allepith_clustering), names(assignments_normsig_both_epithelial))
pData(HSMM_allepith_clustering)$assignments_normsig_both <- assignments_normsig_both_epithelial
pData(HSMM_allepith_clustering)$assignments_normsig_ups <- assignments_normsig_ups_epithelial

#pdf(here("plots", "fig3c.pdf"), width = 10)
plot_cell_clusters(HSMM_allepith_clustering, 1, 2, color = "assignments_normsig_both", cell_size = 2) + facet_wrap(~patient)
#dev.off()

#pdf(here("plots", "figS7.pdf"), width = 10)
plot_cell_clusters(HSMM_allepith_clustering, 1, 2, color = "assignments_normsig_both", cell_size = 2) + facet_wrap(~Cluster)
#dev.off()


# dot plots normal signatures
all.equal(HSMM_allepith_clustering$Cluster, clustering_allepith)
all.equal(colnames(normsig_avg_both_epithelial), colnames(HSMM_allepith_clustering))

clust_avg_normsig_both <- matrix(NA, nrow = length(unique(HSMM_allepith_clustering$Cluster)), ncol = nrow(normsig_avg_both_epithelial))
rownames(clust_avg_normsig_both) <- paste("clust", c(1:length(unique(HSMM_allepith_clustering$Cluster))), sep = "")
colnames(clust_avg_normsig_both) <- rownames(normsig_avg_both_epithelial)
for (c in 1:length(unique(HSMM_allepith_clustering$Cluster))) {
  clust_avg_normsig_both[c,] <- apply(normsig_avg_both_epithelial[,which(HSMM_allepith_clustering$Cluster == c)], 1, mean)
}

clust_avg_normsig_both <- as.data.frame(clust_avg_normsig_both)
clust_avg_normsig_both$Cluster <- rownames(clust_avg_normsig_both)
clust_avg_normsig_melt <- melt(clust_avg_normsig_both, "Cluster")

#pdf(here("plots", "fig3c.pdf"), width = 6.5)
ggplot(clust_avg_normsig_melt, aes(Cluster, value, fill = factor(variable), color = factor(variable), 
                                   shape = factor(variable))) + 
  geom_point(size = 3, stroke = 1) +
  scale_shape_discrete(solid = T) + 
  ylab("average expression of signature in cluster") +
  xlab("cluster") +
  ylim(c(-0.35, 0.5))
#dev.off()


## Mammaprint signature
mammaprint_long <- read.table(here("data", "mammaprint_sig_new.txt"), header = TRUE, sep = "\t")
mammaprint <- apply(mammaprint_long, 2, function(x){return(intersect(x, rownames(mat_ct)))})[,1]
mammaprint_avg_exprs <- apply(mat_ct[match(mammaprint, rownames(mat_ct)),], 2, mean)
all.equal(names(mammaprint_avg_exprs), colnames(mat_ct))
mammaprint_avg_exprs <- mammaprint_avg_exprs[which(pd_ct$cell_types_cl_all == "epithelial")]

# per patient and per cluster Mammaprint plots
all.equal(colnames(HSMM_allepith_clustering), names(mammaprint_avg_exprs))
pData(HSMM_allepith_clustering)$mammaprint_avg_exprs <- mammaprint_avg_exprs

#pdf(here("plots", "figS13b.pdf"), width = 10)
plot_cell_clusters(HSMM_allepith_clustering, 1, 2, color = "mammaprint_avg_exprs", cell_size = 2) + facet_wrap(~patient) +
  scale_color_continuous(low = "yellow", high = "blue")
#dev.off()

#pdf(here("plots", "figS13a.pdf"), width = 10)
plot_cell_clusters(HSMM_allepith_clustering, 1, 2, color = "mammaprint_avg_exprs", cell_size = 2) + facet_wrap(~Cluster) +
  scale_color_continuous(low = "yellow", high = "blue")
#dev.off()


## zena werb signature
zenawerb_long <- read.table(here("data", "werb_49_metastasis_sig.txt"), header = TRUE, sep = "\t")
zenawerb <- apply(zenawerb_long, 2, function(x){return(intersect(x, rownames(mat_ct)))})[,1]
zenawerb_avg_exprs <- apply(mat_ct[match(zenawerb, rownames(mat_ct)),], 2, mean)
all.equal(names(zenawerb_avg_exprs), colnames(mat_ct))
zenawerb_avg_exprs <- zenawerb_avg_exprs[which(pd_ct$cell_types_cl_all == "epithelial")]

# per patient and per cluster Zena Werb Signature plots
all.equal(colnames(HSMM_allepith_clustering), names(zenawerb_avg_exprs))
pData(HSMM_allepith_clustering)$zenawerb_avg_exprs <- zenawerb_avg_exprs

#pdf(here("plots", "figS14b.pdf"), width = 10)
plot_cell_clusters(HSMM_allepith_clustering, 1, 2, color = "zenawerb_avg_exprs", cell_size = 2) + facet_wrap(~patient) +
  scale_color_continuous(low = "yellow", high = "blue")
#dev.off()

#pdf(here("plots", "figS14a.pdf"), width = 10)
plot_cell_clusters(HSMM_allepith_clustering, 1, 2, color = "zenawerb_avg_exprs", cell_size = 2) + facet_wrap(~Cluster) +
  scale_color_continuous(low = "yellow", high = "blue")
#dev.off()


## carlos artega signature
artega_long <- read.table(here("data", "artega_sig.txt"), header = TRUE, sep = "\t")
artega <- apply(artega_long, 2, function(x){return(intersect(x, rownames(mat_ct)))})[,1]
artega_avg_exprs <- apply(mat_ct[match(artega, rownames(mat_ct)),], 2, mean)
all.equal(names(artega_avg_exprs), colnames(mat_ct))
artega_avg_exprs <- artega_avg_exprs[which(pd_ct$cell_types_cl_all == "epithelial")]

# per patient and per cluster Artega signature plots
all.equal(colnames(HSMM_allepith_clustering), names(artega_avg_exprs))
pData(HSMM_allepith_clustering)$artega_avg_exprs <- artega_avg_exprs

#pdf(here("plots", "figS15a.pdf"), width = 10)
plot_cell_clusters(HSMM_allepith_clustering, 1, 2, color = "artega_avg_exprs", cell_size = 2) + facet_wrap(~patient) +
  scale_color_continuous(low = "yellow", high = "blue")
#dev.off()

#pdf(here("plots", "figS15b.pdf"), width = 10)
plot_cell_clusters(HSMM_allepith_clustering, 1, 2, color = "artega_avg_exprs", cell_size = 2) + facet_wrap(~Cluster) +
  scale_color_continuous(low = "yellow", high = "blue")
#dev.off()


## heatmap on prognosis signatures
prognosis_sig <- cbind(mammaprint_avg_exprs, zenawerb_avg_exprs, artega_avg_exprs)
colnames(prognosis_sig) <- c("mammaprint", "zenawerb", "artega")

# mammaprint, zenawerb and carlos artega signatures separate per patient
prognosis_epith_pat <- list()
for (i in 1:length(patients_now)) {
  prognosis_epith_pat[[i]] <- t(prognosis_sig)[,which(pd_ct_epith$patient == patients_now[i])]
}
names(prognosis_epith_pat) <- patients_now
for (i in 1:length(patients_now)) {
  print(all.equal(colnames(prognosis_epith_pat[[1]]), rownames(pds_epith_ct[[1]])))
  print(all.equal(names(clusterings_sep_allepith[[1]]), colnames(prognosis_epith_pat[[1]])))
}

ht_sep_prognosis <-
  Heatmap(prognosis_epith_pat[[1]],
          cluster_rows = FALSE,
          col = colorRamp2(c(-0.2, 0.2, 1), c("blue","white", "red")),
          show_column_names = FALSE,
          column_title = patients_now[1],
          top_annotation = ha_lehman_epith_pat[[1]],
          column_title_gp = gpar(fontsize = 12),
          show_row_names = FALSE,
          name = patients_now[1], 
          top_annotation_height = unit(c(2), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9))) +
  Heatmap(prognosis_epith_pat[[2]],
          col = colorRamp2(c(-0.2, 0.2, 1), c("blue","white", "red")),
          cluster_rows = FALSE,
          show_column_names = FALSE,
          column_title = patients_now[2],
          column_title_gp = gpar(fontsize = 12),
          top_annotation = ha_lehman_epith_pat[[2]],
          name = patients_now[2], 
          show_row_names = FALSE,
          top_annotation_height = unit(c(2), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9))) +
  Heatmap(prognosis_epith_pat[[3]],
          col = colorRamp2(c(-0.2, 0.2, 1), c("blue","white", "red")),
          cluster_rows = FALSE,
          show_column_names = FALSE,
          column_title = patients_now[3],
          column_title_gp = gpar(fontsize = 12),
          top_annotation = ha_lehman_epith_pat[[3]],
          name = patients_now[3], 
          show_row_names = FALSE,
          top_annotation_height = unit(c(2), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9))) +
  Heatmap(prognosis_epith_pat[[4]],
          col = colorRamp2(c(-0.2, 0.2, 1), c("blue","white", "red")),
          cluster_rows = FALSE,
          show_column_names = FALSE,
          column_title = patients_now[4],
          column_title_gp = gpar(fontsize = 12),
          top_annotation = ha_lehman_epith_pat[[4]],
          name = patients_now[4], 
          show_row_names = FALSE,
          top_annotation_height = unit(c(2), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9))) +
  Heatmap(prognosis_epith_pat[[5]],
          col = colorRamp2(c(-0.2, 0.2, 1), c("blue","white", "red")),
          cluster_rows = FALSE,
          show_column_names = FALSE,
          column_title = patients_now[5],
          column_title_gp = gpar(fontsize = 12),
          top_annotation = ha_lehman_epith_pat[[5]],
          name = patients_now[5], 
          show_row_names = FALSE,
          top_annotation_height = unit(c(2), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9))) +
  Heatmap(prognosis_epith_pat[[6]],
          col = colorRamp2(c(-0.2, 0.2, 1), c("blue","white", "red")),
          cluster_rows = FALSE,
          row_names_side = "right",
          column_title = patients_now[6],
          column_title_gp = gpar(fontsize = 12),
          top_annotation = ha_lehman_epith_pat[[6]],
          name = patients_now[6], 
          show_column_names = FALSE,
          top_annotation_height = unit(c(2), "cm"),
          heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 9)))
#pdf(here("plots", "fig4a.pdf"), onefile = FALSE, width = 20)
print(draw(ht_sep_prognosis, annotation_legend_side = "bottom"))
#dev.off()


# violin ranking plot for prognosis signatures
clust_avg_prognosis <- matrix(NA, nrow = length(unique(HSMM_allepith_clustering$Cluster)), ncol = ncol(prognosis_sig))
rownames(clust_avg_prognosis) <- paste("clust", c(1:length(unique(HSMM_allepith_clustering$Cluster))), sep = "")
colnames(clust_avg_prognosis) <- colnames(prognosis_sig)
for (c in 1:length(unique(HSMM_allepith_clustering$Cluster))) {
  clust_avg_prognosis[c,] <- apply(prognosis_sig[which(HSMM_allepith_clustering$Cluster == c),], 2, mean)}

prognosis_sig <- as.data.frame(prognosis_sig)
all.equal(rownames(prognosis_sig), colnames(HSMM_allepith_clustering))
prognosis_sig$Cluster <- as.numeric(HSMM_allepith_clustering$Cluster)

prognosis_melt <- melt(prognosis_sig, id.vars = c("Cluster"))
prognosis_melt$value <- as.numeric(prognosis_melt$value)
prognosis_melt <- prognosis_melt %>%
  filter(Cluster %in% c(2,3,4))

#pdf(here("plots", "fig4b.pdf"), width = 9)
p <- ggplot(prognosis_melt, aes(factor(Cluster), value, fill = factor(Cluster))) +
  scale_fill_manual(values = c("1" = "#ee204d", "2" = "#17806d", "3" = "#b2ec5d", "4" = "#cda4de", "5" = "#1974d2"))
p + geom_violin(adjust = .5) + facet_wrap(~variable) + stat_summary(fun.y = mean, geom = "point", shape = 22, size = 3)
#dev.off()


## venn diagram signatures intersection
#pdf(here("plots", "figS16.pdf"))
venn(list("PS" = mammaprint, "MBS" = zenawerb, "RTS" = artega))
#dev.off()


## t-sne clustering from 3a, with different colors of dots
# repeat plot 3a with normal signatures
#pdf(here("plots", "figS10b.pdf"), width = 10, height = 10)
plot_cell_clusters(HSMM_allepith_clustering, 1, 2, color = "assignments_normsig_both", cell_size = 3)
#dev.off()

# repeat plot 3a with TNBCtype4 Lehman signatures
lehman_avg_both_epithelial_new <- lehman_avg_both_epithelial[-which(rownames(lehman_avg_both_epithelial) %in% c("immunomodulatory", "mesenchymal_stem_like")),]
assignments_lehman_both_epithelial_new <- apply(lehman_avg_both_epithelial_new, 2, function(x){rownames(lehman_avg_both_epithelial_new)[which.max(x)]})
all.equal(colnames(HSMM_allepith_clustering), names(assignments_lehman_both_epithelial_new))
pData(HSMM_allepith_clustering)$assignments_lehman_both_new <- assignments_lehman_both_epithelial_new

#pdf(here("plots", "figS10c.pdf"), width = 10, height = 10)
plot_cell_clusters(HSMM_allepith_clustering, 1, 2, color = "assignments_lehman_both_new", cell_size = 3)
#dev.off()

# repeat plot 3a with mammaprint signature (PS)
#pdf(here("plots", "figS10d.pdf"), width = 10, height = 10)
plot_cell_clusters(HSMM_allepith_clustering, 1, 2, color = "mammaprint_avg_exprs", cell_size = 3) +
  scale_color_continuous(low = "yellow", high = "blue")
#dev.off()


## heatmap for pathways
path1 <- c("ST3GAL4", "ST3GAL6", "ST8SIA1", "FUT1", "FUT2", "FUT3", "FUT4", "FUT6", "FUT9", "GLTP", "SPTLC1", "KDSR", "SPTLC2", "CERK", "ARSA")
idx_path1 <- match(path1, rownames(HSMM_allepith_clustering))

path2 <- c("CCL20", "CCL22", "CCL4", "CCR6", "IL11", "IL12RB1", "IL6R", "CSF2RA", "BMP7", "GLMN", "GPI", "INHBA", "TNF", 
           "TNFSF15", "GHR", "LEPR", "TLR1", "TLR2", "TLR5", "TLR7", "TLR10", "F11R")
idx_path2 <- match(path2, rownames(HSMM_allepith_clustering))

path3 <- c("ERBB2", "AKT1", "AKT3", "PIK3R3", "PIK3R4", "RPS6KB2", "TRIB3", "BTK", "GRB10", "GRB2", "ILK", "PAK1", "PRKCZ", "CSNK2A1",
           "IRS1", "IRAK1", "MYD88", "MAP2K1", "MAPK8", "MAPK1", "PTPN11", "EIF4E", "EIF4EBP1", "EIF4G1", "FKBP1A", "PDK1", "RHEB", "RPS6KB1")
idx_path3 <- match(path3, rownames(HSMM_allepith_clustering))

all_idx_paths <- c(idx_path1, idx_path2)

names_path <- c(rep("glycosphigolipid metabolism", length(idx_path1)), rep("innate immunity", length(idx_path2)))
anno_cols_path <- c("glycosphigolipid metabolism" = "#ff9baa", "innate immunity" = "#17806d")
ha_path_rows <- HeatmapAnnotation(df = data.frame(pathway = names_path),
                                  annotation_legend_param = list(pathway = list(ncol = 1, title = "pathway", title_position = "topcenter")),
                                  which = "row", col = list("pathway" = anno_cols_path), annotation_width = unit(3, "mm"))

# separate matrices per cluster and only extract the relevant genes
mat_epith_allpats <- exprs(HSMM_allepith_clustering)
mats_epith_patient <- list()
pds_epith_patient <- list()
mats_epith_patient_clusters <- list()
for (i in 1:length(patients_now)) {
  pds_epith_patient[[i]] <- pData(HSMM_allepith_clustering)[which(pData(HSMM_allepith_clustering)$patient == patients_now[i]), ]
  mats_epith_patient[[i]] <- mat_epith_allpats[all_idx_paths, which(pData(HSMM_allepith_clustering)$patient == patients_now[i])]
  mats_epith_patient_clusters[[i]] <- list()
  for (c in 1:length(unique(HSMM_allepith_clustering$Cluster))) {
    mats_epith_patient_clusters[[i]][[c]] <- mats_epith_patient[[i]][, which(pds_epith_patient[[i]]$Cluster == c)]
  }
  names(mats_epith_patient_clusters[[i]]) <- paste0("clust", c(1:5))
}
names(mats_epith_patient) <- patients_now
names(pds_epith_patient) <- patients_now

# heatmap
ht_paths <-
  ha_path_rows + 
  Heatmap(mats_epith_patient_clusters[[1]][[1]], 
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          show_column_dend = TRUE, show_column_names = FALSE,
          name = "cluster1", column_title = "cluster1",
          row_names_gp = gpar(fontsize = 9, col = anno_cols_path),
          split = names_path) +
  Heatmap(mats_epith_patient_clusters[[1]][[2]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE,
          name = "cluster2", column_title = "cluster2",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[1]][[3]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster3", column_title = "cluster3",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[1]][[4]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster4", column_title = "cluster4",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[1]][[5]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster5", column_title = "cluster5",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[2]][[1]], 
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          show_column_dend = TRUE, show_column_names = FALSE,
          name = "cluster1_2", column_title = "cluster1",
          row_names_gp = gpar(fontsize = 9, col = anno_cols_path),
          split = names_path) +
  Heatmap(mats_epith_patient_clusters[[2]][[2]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE,
          name = "cluster2_2", column_title = "cluster2",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[2]][[3]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster3_2", column_title = "cluster3",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[2]][[4]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster4_2", column_title = "cluster4",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[2]][[5]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster5_2", column_title = "cluster5",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) + 
  Heatmap(mats_epith_patient_clusters[[3]][[1]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE,
          name = "cluster1_3", column_title = "cluster1",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[3]][[2]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE,
          name = "cluster2_3", column_title = "cluster2",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[3]][[3]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster3_3", column_title = "cluster3",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[3]][[4]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster4_3", column_title = "cluster4",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[3]][[5]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster5_3", column_title = "cluster5",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[4]][[1]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE,
          name = "cluster1_4", column_title = "cluster1",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[4]][[2]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE,
          name = "cluster2_4", column_title = "cluster2",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[4]][[3]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster3_4", column_title = "cluster3",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[4]][[4]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster4_4", column_title = "cluster4",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[4]][[5]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster5_4", column_title = "cluster5",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[5]][[1]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE,
          name = "cluster1_5", column_title = "cluster1",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[5]][[2]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE,
          name = "cluster2_5", column_title = "cluster2",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[5]][[3]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster3_5", column_title = "cluster3",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[5]][[4]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster4_5", column_title = "cluster4",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[5]][[5]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster5_5", column_title = "cluster5",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[6]][[1]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE,
          name = "cluster1_6", column_title = "cluster1",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[6]][[2]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE,
          name = "cluster2_6", column_title = "cluster2",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[6]][[3]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster3_6", column_title = "cluster3",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[6]][[4]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster4_6", column_title = "cluster4",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE) +
  Heatmap(mats_epith_patient_clusters[[6]][[5]],
          col = colorRamp2(c(-0.5, 1, 3), c("blue", "white", "red")),
          cluster_rows = TRUE, show_row_dend = FALSE, 
          name = "cluster5_6", column_title = "cluster5",
          show_row_names = FALSE,
          show_column_dend = TRUE, show_column_names = FALSE)
#pdf(here("plots", "fig4d.pdf"), onefile = FALSE, width = 25)
draw(ht_paths)
#dev.off()

