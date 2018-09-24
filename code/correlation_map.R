## this script creates expression correlation maps for normal (Gao et al 2017 Nat Comm) and TNBC cells
## this script generates the following figures
##    fig. 2f
##    fig. S12

library(here)
source(here::here("code", "libraries.R"))

mat_complete <- readRDS(here::here("data","mat_complete.RDS"))
pd_ct <- readRDS(here::here("data", "pd_ct.RDS"))
mat_ct <- readRDS(here::here("data", "mat_ct.RDS"))
order_samples_cnv <- readRDS(here::here("data", "order_samples_cnv.RDS"))
pd_epith <- readRDS(here::here("data", "pd_epith.RDS"))

patients_now <- sort(unique(pd_ct$patient))

# order epithelial cells by patient, in the same order as they appear in the CNV plot in 2d, and add normal cells at the top
ord_ct <- c(intersect(which(pd_ct$cell_types_cl_all == "epithelial"), which(pd_ct$patient == patients_now[1]))[order_samples_cnv[[3]]],
            intersect(which(pd_ct$cell_types_cl_all == "epithelial"), which(pd_ct$patient == patients_now[2]))[order_samples_cnv[[4]]],
            intersect(which(pd_ct$cell_types_cl_all == "epithelial"), which(pd_ct$patient == patients_now[3]))[order_samples_cnv[[5]]],
            intersect(which(pd_ct$cell_types_cl_all == "epithelial"), which(pd_ct$patient == patients_now[4]))[order_samples_cnv[[6]]],
            intersect(which(pd_ct$cell_types_cl_all == "epithelial"), which(pd_ct$patient == patients_now[5]))[order_samples_cnv[[7]]],
            intersect(which(pd_ct$cell_types_cl_all == "epithelial"), which(pd_ct$patient == patients_now[6]))[order_samples_cnv[[8]]])
ord_norm <- order_samples_cnv[[2]]
ord_ct <- ord_ct + length(ord_norm)
ord_ct <- c(ord_norm, ord_ct)
mat_cor <- mat_complete[, ord_ct]

ct_ord_per_pat <- pd_ct$cell_types_cl_all
ct_ord_per_pat[intersect(which(pd_ct$cell_types_cl_all == "epithelial"), which(pd_ct$patient == patients_now[1]))][order_samples_cnv[[3]]] <- paste("epith", patients_now[1], sep = "")
ct_ord_per_pat[intersect(which(pd_ct$cell_types_cl_all == "epithelial"), which(pd_ct$patient == patients_now[2]))][order_samples_cnv[[4]]] <- paste("epith", patients_now[2], sep = "")
ct_ord_per_pat[intersect(which(pd_ct$cell_types_cl_all == "epithelial"), which(pd_ct$patient == patients_now[3]))][order_samples_cnv[[5]]] <- paste("epith", patients_now[3], sep = "")
ct_ord_per_pat[intersect(which(pd_ct$cell_types_cl_all == "epithelial"), which(pd_ct$patient == patients_now[4]))][order_samples_cnv[[6]]] <- paste("epith", patients_now[4], sep = "")
ct_ord_per_pat[intersect(which(pd_ct$cell_types_cl_all == "epithelial"), which(pd_ct$patient == patients_now[5]))][order_samples_cnv[[7]]] <- paste("epith", patients_now[5], sep = "")
ct_ord_per_pat[intersect(which(pd_ct$cell_types_cl_all == "epithelial"), which(pd_ct$patient == patients_now[6]))][order_samples_cnv[[8]]] <- paste("epith", patients_now[6], sep = "")
ct_ord_per_pat <- c(rep("normal", length(ord_norm)), ct_ord_per_pat) 
ct_ord <- ct_ord_per_pat[ord_ct]

# compute correlation map
cors_ct <- cor(mat_cor, method = "pearson")

# prepare for plotting
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

tsne_cols <- c("epithelial" = unname(basal_epithelial_col), "stroma" = unname(stroma_col), "endothelial" = unname(endothelial_col),
               "Tcell" = unname(t_cell_col), "Bcell" = unname(b_cell_col), "macrophage" = unname(macrophage_col))

pats_cols <- c("PT039" = unname(brocolors("crayons")["Orange Red"]), "PT058" = unname(brocolors("crayons")["Orange"]), 
               "PT081" = unname(brocolors("crayons")["Pink Flamingo"]), "PT084" = unname(brocolors("crayons")["Fern"]), 
               "PT089" = unname(brocolors("crayons")["Blue Violet"]), "PT126" = unname(brocolors("crayons")["Sky Blue"]))

anno_colors <- list("tsne" = tsne_cols, "patient" = pats_cols)

cols_ct_pats <- c(anno_colors$tsne, anno_colors$patient)
names(cols_ct_pats)[(length(cols_ct_pats) - length(patients_now) + 1):length(cols_ct_pats)] <- paste("epith",patients_now, sep = "")
cols_ct_pats <- c("normal" = "gray", cols_ct_pats)
if (length((which(names(cols_ct_pats) == "epithelial")) > 0))
  cols_ct_pats <- cols_ct_pats[-(which(names(cols_ct_pats) == "epithelial"))]
cols_ct_pats <- cols_ct_pats[match(c("normal", "stroma", "Tcell", "Bcell", "macrophage", "endothelial", paste("epith", patients_now, sep = "")), names(cols_ct_pats))]


# cluster annotation
cluster_order_allepith <- pd_epith$Cluster[match(colnames(mat_cor)[(length(ord_norm) + 1):length(colnames(mat_cor))], pd_epith$Sample)]

ha_cors_top <- HeatmapAnnotation(data.frame(celltype = ct_ord, cluster_all = c(rep(0, length(ord_norm)), cluster_order_allepith)), 
                                 col = list(celltype = cols_ct_pats, cluster_all = c("0" = "gray", "1" = "#ee204d", "2" = "#17806d", "3" = "#b2ec5d", "4" = "#cda4de", "5" = "#1974d2")),
                                 show_annotation_name = TRUE,
                                 annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 8),
                                 annotation_legend_param = list(list(title_position = "topcenter", title = c("cell type")),
                                                                list(title_position = "topcenter", title = "cluster_all")),
                                 show_legend = TRUE,
                                 gap = unit(c(2), "mm"))

ha_cors_rows <- HeatmapAnnotation(data.frame(celltype = ct_ord, cluster_all = c(rep(0, length(ord_norm)), cluster_order_allepith)), 
                                  col = list(celltype = cols_ct_pats, cluster_all = c("0" = "gray", "1" = "#ee204d", "2" = "#17806d", "3" = "#b2ec5d", "4" = "#cda4de", "5" = "#1974d2")),
                                  show_annotation_name = FALSE,
                                  annotation_legend_param = list(list(title_position = "topcenter", title = c("cell type")),
                                                                 list(title_position = "topcenter", title = "cluster_all")),
                                  show_legend = FALSE,
                                  which = "row",
                                  gap = unit(c(2), "mm"))


ht_cors <- ha_cors_rows + 
  Heatmap(cors_ct, 
          name = "expr cor",
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          show_row_names = FALSE,
          show_column_names = FALSE,
          top_annotation = ha_cors_top,
          col = colorRamp2(c(0, 0.2, 0.3, 0.4), c("lightgray", "pink", "mediumorchid3", "purple")),
          use_raster = TRUE, raster_device = "CairoPNG")
#pdf(here::here("plots", "fig2f.pdf"), onefile = FALSE, width = 12, height = 10)
#pdf(here::here("plots", "figS12.pdf"), onefile = FALSE, width = 12, height = 10)
#print(draw(ht_cors))
#dev.off()
