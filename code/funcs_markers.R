## functions to identify cell types based on markers

# function to find which markers are expressed by cells in an expression matrix
# returns list of lists, length of each list = the number of cells in the matrix
lists_markers <- function(mat_here, thresh, markers){
  
  epithelial_cells <- list()
  immune_cells <- list()
  other_cells <- list()
  
  for (j in 1:ncol(mat_here)) {
    
    x <- mat_here[, j]
    
    vec_epithelial <- c()
    vec_immune <- c()
    vec_other <- c()
    
    for (i in 1:nrow(markers)) {
      
      name_marker <- markers[i, "gene"]
      eval(parse(text = paste(name_marker, "_in <- which(rownames(mat_here) == ", 'name_marker', ")", sep = "")))
      
      eval(parse(text = paste("if (x[", name_marker,"_in] > thresh) {vec_", markers[i, "type_long"], " <- c(vec_", markers[i, "type_long"], ", ", 'name_marker', 
                              "); names(vec_", markers[i, "type_long"], ")[length(vec_", markers[i, "type_long"], ")] <- ", 
                              'markers[i, "type"]', "}", sep = "")))
    }
    
    epithelial_cells[[length(epithelial_cells) + 1]] <- if (length(vec_epithelial)) vec_epithelial else list(NULL)
    immune_cells[[length(immune_cells) + 1]] <- if (length(vec_immune)) vec_immune else list(NULL)
    other_cells[[length(other_cells) + 1]] <- if (length(vec_other)) vec_other else list(NULL)
    
  }
  
  names(epithelial_cells) <- colnames(mat_here)
  epithelial_cells <- lapply(epithelial_cells, function(x){unlist(x)})
  
  names(immune_cells) <- colnames(mat_here)
  immune_cells <- lapply(immune_cells, function(x){unlist(x)})
  
  names(other_cells) <- colnames(mat_here)
  other_cells <- lapply(other_cells, function(x){unlist(x)})
  
  return(lists_cells = list("epithelial_cells" = epithelial_cells, "immune_cells" = immune_cells, "other_cells" = other_cells))
}


# decide whether the cells with markers are epithelial by the rules
# returns vector binary, length number of cells in the list
decide_is_epithelial <- function(list_epithelial_markers){
  
  mat_is_epithelial <- sapply(list_epithelial_markers, function(x){
    if (length(x) >= 2)
      return(1)
    return(0)
  })
  
  return(mat_is_epithelial)
}


# decide whether the cells with markers are immune by the rules
# returns vector binary, length number of cells in the list
decide_is_immune <- function(list_immune_markers){
  
  mat_is_immune <- sapply(list_immune_markers, function(x){
    
    # no marker was found
    if (length(x) == 0)
      return(0)
    # only one single marker (either specific or unspecific) --> not enough evidence
    if (length(unique(x)) == 1)
      return(0)
    # only one type with at least 2 markers (without PTPRC) --> that type
    if (length(x) > 1 & length(unique(names(x))) == 1)
      return(unique(names(x)))
    # only one type with at least 1 marker, and PTPRC --> that type
    if (length(unique(names(x))) == 2 && "PTPRC" %in% unique(names(x)))
      return(setdiff(unique(names(x)), "PTPRC"))
    # more markers of different types
    if (length(unique(names(x))) > 1) {
      names_cells <- sort(table(names(x)), decreasing = TRUE)
      # one marker is prevalent
      if (names_cells[1] >= 3 && names_cells[2] < 2)
        return(names(names_cells)[1])
      else
        return("immune_mix")
    }
  })
  
  return(mat_is_immune)
}


# decide whether the cells with markers are stroma, endothelial or adipocytes by the rules
# returns vector binary, length number of cells in the list
decide_is_other <- function(list_other_markers){
  
  mat_is_other <- sapply(list_other_markers, function(x){
    # no marker was found
    if (length(x) == 0)
      return(0)
    # only one marker was found --> not enough evidence
    if (length(x) == 1)
      return(0)
    # more than one marker of the same type --> that type
    if (length(x) > 1 & length(unique(names(x))) == 1)
      return(unique(names(x)))
    # both types
    if (length(unique(names(x))) > 1) {
      names_cells <- sort(table(names(x)), decreasing = TRUE)
      # one marker is prevalent
      if (names_cells[1] >= 3 && names_cells[2] < 2)
        return(names(names_cells)[1])
      else
        return("other_mix")
    }
  })
  return(mat_is_other)
}


# function to evaluate the distribution of expression for the cases in which cells only have one epithelial marker expressed
# type_test indicates whether we should look at the distribution of that markers in the cells from only one patient or from all patients
# thrsh_percent is the how high the expression should be, as a fraction
expression_one_epithelial_marker <- function(mat_here, pd_here, is_epithelial, epithelial_markers, type_test, thresh_percent){
  
  cdfspats_is_epithelial <- list()
  cdfsall_is_epithelial <- list()
  exprs_is_epithelial <- list()
  for (i in 1:length(is_epithelial)) {
    # if the cell is not already epithelial
    if (is_epithelial[i] == 0) {
      # but has some epithelial markers
      if (!is.null(epithelial_markers[[i]])) {
        name_cell <- names(epithelial_markers)[i]
        patient <- pd_here$patient[grep(name_cell, colnames(mat_here))]
        all_cells_pat <- colnames(mat_here)[which(pd_here$patient == patient)]
        # indices of the current patient among the cells that are epithelial by markers
        idx_cells_pat_here <- na.omit(match(all_cells_pat, colnames(mat_here)[which(is_epithelial == 1)]))
        
        cdfs_pats <- rep(NA, length(epithelial_markers[[i]]))
        names(cdfs_pats) <- epithelial_markers[[i]]
        cdfs_all <- rep(NA, length(epithelial_markers[[i]]))
        names(cdfs_all) <- epithelial_markers[[i]]
        exprs_j_cell <- rep(NA, length(epithelial_markers[[i]]))
        names(exprs_j_cell) <- epithelial_markers[[i]]
        
        for (j in 1:length(epithelial_markers[[i]])) {
          # expression of the current cell for the current epithelial marker
          exprs_j_cell[j] <- mat_here[which(rownames(mat_here) == epithelial_markers[[i]][j]), i]
          exprs_j_pats <- mat_here[,which(is_epithelial == 1)][which(rownames(mat_here) == epithelial_markers[[i]][j]), idx_cells_pat_here]
          w_pats <- ecdf(exprs_j_pats)
          cdfs_pats[j] <- w_pats(exprs_j_pats[j])
          
          exprs_j_all <- mat_here[,which(is_epithelial == 1)][which(rownames(mat_here) == epithelial_markers[[i]][j]), ]
          w_all <- ecdf(exprs_j_all)
          cdfs_all[j] <- w_all(exprs_j_all[j])
        }
        
        exprs_is_epithelial[[i]] <- exprs_j_cell
        cdfspats_is_epithelial[[i]] <- cdfs_pats
        cdfsall_is_epithelial[[i]] <- cdfs_all
      }
    }
  }
  
  for (i in 1:length(cdfsall_is_epithelial)) {
    if (type_test == "all") {
      is_epithelial_extra <- sapply(cdfsall_is_epithelial, function(x){if (!is.null(x) && names(x) %in% c("EPCAM", markers$gene[grep("KRT", markers$gene)]) && x > thresh_percent) return(1) else return(0)})
    }
    
    if (type_test == "pats") {
      is_epithelial_extra <- sapply(cdfspats_is_epithelial, function(x){if (!is.null(x) && names(x) %in% c("EPCAM", markers$gene[grep("KRT", markers$gene)]) && x > thresh_percent) return(1) else return(0)})
    }
  }
  
  return(list("exprs_is_epithelial" = exprs_is_epithelial, "cdfspats_is_epithelial" = cdfspats_is_epithelial, 
              "cdfsall_is_epithelial" = cdfsall_is_epithelial, "is_epithelial_extra" = is_epithelial_extra, "type_test" = type_test, "thresh_percent" = thresh_percent))
  
}

