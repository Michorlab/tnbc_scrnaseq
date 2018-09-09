# this script performs the following analyses from the manuscript: 
#   cycling cells identification
#   cell types identification (via markers and clustering)
# this script plots the following figures from the manuscript: 
#   fig. 1b (heatmap)
#   fig. 1c (barplot)
#   fig 1d (barplot)

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
# the original cluster from the paper. for reproducibility, we read in the original clustering assignment
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
for (i in 1:length(patients_now)){
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
  pdf(here("plots", paste("cycling_patient", patients_now[i], ".pdf", sep = "")))
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
  dev.off()
}


#####
# barplot cell cycle phase for all patients
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






