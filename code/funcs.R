## functions for various aspects of cell type identification, clustering, signatures etc

# function to match and clean a vector of genes
match_clean_vector_genes <- function(mat_to_fit, vec_genes){
  
  vec_genes <- unique(vec_genes)
  vec_genes <- trimws(vec_genes)
  vec_genes <- data.frame("gene" = vec_genes)
  vec_genes$index <- match(vec_genes$gene, rownames(mat_to_fit))
  if (length(which(is.na(vec_genes$index)) > 0))
    vec_genes <- vec_genes[-which(is.na(vec_genes$index)), ]
  return(vec_genes)
}



# function to average the expression over a vector of gene indices, for all cells
avg_expr_genes <- function(mat_to_fit, indices){
  scores_cells <- apply(mat_to_fit, 2, function(x){mean(x[indices])})
  return(scores_cells)
}



# function to compute monocle clustering and produce plots
monocle_unsup_clust_plots <- function(sceset_obj, slot = NULL, mat_to_cluster, anno_colors, name_in_phenodata, 
                                      disp_extra, save_plots, path_plots = NULL, type_pats, regress_pat = 0, 
                                      num_clusters = NULL, perplexity = NULL, use_known_colors = NULL, use_only_known_celltypes = NULL){
  
  if (!is.null(slot) && slot == "exprs")
    mat_to_cluster <- assays(sceset_obj)$exprs
  if (!is.null(slot) && slot != "exprs")
    mat_to_cluster <- get_exprs(sceset_obj, slot)
  
  pd <- new("AnnotatedDataFrame", data = as.data.frame(colData(sceset_obj)))
  featureNames(pd) <- rownames(pd)
  fd <- new("AnnotatedDataFrame", data = as.data.frame(rowData(sceset_obj)))
  featureNames(fd) <- rowData(sceset_obj)$gene_short_name
  
  
  HSMM_clustering <- newCellDataSet(data.matrix(mat_to_cluster), phenoData = pd, featureData = fd, 
                                    expressionFamily = negbinomial.size(),
                                    lowerDetectionLimit = 0.1)
  if (length(grep(name_in_phenodata, colnames(pData(HSMM_clustering)))) > 0) {
    print("clustering already done")
    return(HSMM_clustering)
  }
  
  if (is.null(num_clusters))
    suffix_clust_plot <- "auto"
  else
    suffix_clust_plot <- num_clusters
  
  if (is.null(use_only_known_celltypes))
    pData(HSMM_clustering)$CellType <- pData(HSMM_clustering)$cell_types_markers
  else
    pData(HSMM_clustering)$CellType <- pData(HSMM_clustering)$cell_types_cl_all
  
  HSMM_clustering <- estimateSizeFactors(HSMM_clustering)
  HSMM_clustering <- estimateDispersions(HSMM_clustering)
  disp_table <- dispersionTable(HSMM_clustering)
  
  if (disp_extra == 0) {
    unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
    suffix_disp_plot <- "nodisp"
  }
  if (disp_extra == 1) {
    unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)
    suffix_disp_plot <- "disp"
  }
  
  HSMM_clustering <- setOrderingFilter(HSMM_clustering, unsup_clustering_genes$gene_id)
  if (save_plots == 1) {
    pdf(paste(path_plots, "/monocle_", type_pats,"_genes_", suffix_disp_plot, ".pdf", sep = ""))
    print(plot_ordering_genes(HSMM_clustering))
    dev.off()
  }
  
  if (regress_pat == 1)
    HSMM_clustering <- reduceDimension(HSMM_clustering, max_components = 2, num_dim = 6,
                                       reduction_method =  "tSNE" , verbose = T, norm_method = "none",
                                       residualModelFormulaStr = "~patient", pseudo_expr = 0)
  else
    if (!is.null(perplexity))
      HSMM_clustering <- reduceDimension(HSMM_clustering, max_components = 2, num_dim = 6,
                                         reduction_method =  "tSNE" , verbose = T, norm_method = "none",
                                         pseudo_expr = 0, perplexity = perplexity)
  else
    HSMM_clustering <- reduceDimension(HSMM_clustering, max_components = 2, num_dim = 6,
                                       reduction_method =  "tSNE" , verbose = T, norm_method = "none",
                                       pseudo_expr = 0)
  
  #plot_pc_variance_explained(HSMM_clustering, return_all = F, norm_method = "none", pseudo_expr = 0)
  HSMM_clustering <- clusterCells(HSMM_clustering, num_clusters = num_clusters)
  pData(HSMM_clustering)$new_cluster <- HSMM_clustering$Cluster
  colnames(pData(HSMM_clustering))[grep("new_cluster", colnames(pData(HSMM_clustering)))] <- name_in_phenodata
  
  
  if (!is.null(use_known_colors))
    cols_cell_types <- anno_colors$tsne
  else
    cols_cell_types <- anno_colors$tsne_unknown
  
  if (save_plots) {
    p1 <- plot_cell_clusters(HSMM_clustering, 1, 2, color = "patient", cell_size = 3) + 
      scale_color_manual(values = anno_colors$patient)
    p2 <- plot_cell_clusters(HSMM_clustering, 1, 2, color = "Cluster", cell_size = 3)
    p3 <- plot_cell_clusters(HSMM_clustering, 1, 2, color = "CellType", cell_size = 3) + 
      scale_color_manual(values = cols_cell_types)
    pdf(paste(path_plots, "/monocle_", type_pats,"_", "regr", regress_pat, "_", suffix_disp_plot, "_clust", suffix_clust_plot, ".pdf", sep = ""), width = 20)
    multiplot(p1, p2, p3, cols = 3)
    dev.off()
  }
  
  plotDF_ct <- data.frame("CellType" = pData(HSMM_clustering)$CellType)
  plotDF_ct$Cluster <- pData(HSMM_clustering)$Cluster
  
  if (save_plots) {
    pdf(paste(path_plots, "/monocle_", type_pats,"_", "regr", regress_pat, "_", suffix_disp_plot, "_clust", suffix_clust_plot, "_perpats.pdf", sep = ""))
    print(plot_cell_clusters(HSMM_clustering, 1, 2, color = "Cluster", cell_size = 3) + facet_wrap(~patient))
    dev.off()
    
    pdf(paste(path_plots, "/monocle_", type_pats,"_", "regr", regress_pat, "_", suffix_disp_plot, "_clust", suffix_clust_plot, "_percelltype.pdf", sep = ""))
    print(plot_cell_clusters(HSMM_clustering, 1, 2, color = "Cluster", cell_size = 3) + facet_wrap(~CellType))
    dev.off()
    
    pdf(paste(path_plots, "/monocle_", type_pats,"_", "regr", regress_pat, "_", suffix_disp_plot, "_clust", suffix_clust_plot, "_percluster.pdf", sep = ""))
    print(plot_cell_clusters(HSMM_clustering, 1, 2, color = "CellType", cell_size = 3) + facet_wrap(~Cluster) + scale_color_manual(values = cols_cell_types))
    dev.off()
    
    pdf(paste(path_plots, "/table_monocle_", type_pats, "_clust", suffix_clust_plot, "_", "regr", regress_pat, "_", suffix_disp_plot, ".pdf", sep = ""))
    print(ggplot(plotDF_ct) + 
            stat_sum(aes(x = factor(Cluster), y = factor(CellType), 
                         size = factor(..n..), color = factor(CellType)), geom = "point") + 
            scale_size_discrete() + 
            labs(x = "Cluster", y = "Cell Type", size = "no. genes") +
            guides(color = FALSE) +
            scale_color_manual(values = cols_cell_types))
    dev.off()
  }
  
  return(HSMM_clustering)
}



# function to compute mean expression of specific markers of cell types
compute_mean_expr_types <- function(types, mat_expr, cells_pos, epithelial_markers, immune_markers, other_markers){
  
  sum_expr <- matrix(nrow = length(cells_pos), ncol = length(types))
  for (i in 1:length(types)) {
    if (types[i] == "epithelial")
      sum_expr[,i] <- t(sapply(c(1:length(epithelial_markers[cells_pos])),function(y){
        x <- epithelial_markers[cells_pos][[y]]
        sum(mat_expr[match(x, rownames(mat_expr)),cells_pos[y]])}))
    if (types[i] == "immune")
      sum_expr[,i] <- t(sapply(c(1:length(immune_markers[cells_pos])),function(y){
        x <- immune_markers[cells_pos][[y]]
        sum(mat_expr[match(x, rownames(mat_expr)),cells_pos[y]])}))
    if (types[i] == "other")
      sum_expr[,i] <- t(sapply(c(1:length(other_markers[cells_pos])),function(y){
        x <- other_markers[cells_pos][[y]]
        sum(mat_expr[match(x, rownames(mat_expr)),cells_pos[y]])}))
  }
  
  colnames(sum_expr) <- types
  return(sum_expr)
  
}



# function to intersect multiple vectors
intersect_all <- function(a,b,...){
  Reduce(intersect, list(a,b,...))
}
