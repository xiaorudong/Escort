
#' Find DE genes between disconnected clusters
#'
#' @param rawcounts A raw count data matrix: row:genes, column:cells
#' @param cls Identified disconnected clusters
#'
#' @import Seurat
#' @import dplyr
#' @export

DE_seurat <- function(rawcounts, cls) {

  a <- match(colnames(rawcounts), names(cls))
  cls <- cls[a]

  clusters <- as.factor(cls)

  scdata <- CreateSeuratObject(counts = rawcounts)
  scdata <- NormalizeData(scdata)
  scdata <- AddMetaData(scdata, metadata=clusters, col.name = "myclusters")
  Idents(scdata) <- scdata$myclusters

  de_ls <- list()
  if(length(unique(cls))==2) {
    de.markers <- FindMarkers(scdata, ident.1 = unique(cls)[1])
    de.markers_sig <- de.markers[de.markers$p_val_adj<0.05, ]
    de.markers_sig$comp <- rep(paste(unique(cls)[1], unique(cls)[2], sep = " vs "), nrow(de.markers_sig))
    de_ls[[unique(cls)[1]]] <- as.data.frame(de.markers_sig)
  }
  if(length(unique(cls))>2) {
    for (cl_i in unique(cls)) {
      de.markers <- FindMarkers(scdata, ident.1 = cl_i)
      de.markers_sig <- de.markers[de.markers$p_val_adj<0.05, ]
      de.markers_sig$comp <- rep(paste(cl_i, "others", sep = " vs "), nrow(de.markers_sig))
      de_ls[[cl_i]] <- as.data.frame(de.markers_sig)
    }
  }
  de_df <- do.call("rbind", de_ls)
  de_df <- de_df %>%
    mutate(across(1:5, round, 3)) %>%
    mutate(across(1:5, function(x) ifelse(x<0.001, "<0.0005", x)))

  return(de_df)
}



#' GSEA genes between disconnected clusters
#'
#' @param norm_counts A normalized count data matrix: row:genes, column:cells
#' @param OrgDb OrgDb
#'
#' @import clusterProfiler
#' @import dplyr
#' @import scran
#' @export

HVGs_GO <- function(norm_counts, OrgDb="org.Hs.eg.db") {
  gene.var <- modelGeneVar(norm_counts)
  HVGs <- getTopHVGs(gene.var, fdr.threshold=0.05)

  ego <- enrichGO(gene = HVGs, OrgDb = OrgDb, ont = "ALL", keyType = "GENENAME")
  if (is.null(ego)) {
    return(NULL)
  }

  if (! is.null(ego)){
    df <- ego@result[, c(1, 3 ,6, 7, 10, 9)]
    df <- as.data.frame(df) %>%
      mutate(across(3:4, round, 3)) %>%
      mutate(across(3:4, function(x) ifelse(x<0.001, "<0.0005", x)))

    colnames(df)[1] <- "ont"
    colnames(df)[4] <- "p.adj"
    colnames(df)[5] <- "N"

    return(df)
  }
}





#' Find HVGs by using scran package
#'
#' @param norm_counts A normalized count data matrix: row:genes, column:cells
#'
#' @import dplyr
#' @import scran
#' @export

HVGs_scran <- function(norm_counts) {
  gene.var <- modelGeneVar(norm_counts)

  df <- as.data.frame(gene.var) %>%
    mutate(across(1:6, round, 3)) %>%
    mutate(across(5:6, function(x) ifelse(x<0.001, "<0.0005", x)))

  # HVGs <- getTopHVGs(gene.var, fdr.threshold=0.05)
  return(df)
}
