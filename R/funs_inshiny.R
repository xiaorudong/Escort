
#' Find top 30 DE genes between disconnected clusters
#'
#' @param rawcounts A raw count data matrix: row:genes, column:cells
#' @param cls Identified disconnected clusters
#'
#' @import edgeR
#' @import dplyr
#' @export

DE_edgeR <- function(rawcounts, cls) {

  a <- match(colnames(rawcounts), names(cls))
  cls <- cls[a]
  meta <- data.frame(cell=names(cls), clusternum=cls)
  meta$cluster <- factor(paste("Cluster", meta$clusternum, sep = ""))
  mm <- model.matrix(~0+cluster, data=meta)
  colnames(mm) <- gsub("cluster", "", colnames(mm))

  combn_tb <- t(combn(levels(meta$cluster), 2))
  comp <- apply(combn_tb, 1, function(x) paste(x,collapse = "-"))
  contr <- makeContrasts(contrasts=comp, levels=mm)

  d <- DGEList(rawcounts)
  d <- calcNormFactors(d)
  d <- estimateDisp(d, mm)
  fit <- glmFit(d, mm)

  de_ls <- list()
  for (i in 1:ncol(contr)) {
    lrt <- glmLRT(fit, contrast = contr[,i])
    top.table <- topTags(lrt, adjust.method = "BH", n = 30)
    top.table <- as.data.frame(top.table)
    colnames(top.table)[which(colnames(top.table)=="PValue")] <- "p"
    colnames(top.table)[which(colnames(top.table)=="FDR")] <- "adj.p"
    top.table$Gene <- rownames(top.table)
    top.table$Comp <- gsub("-", " vs ", colnames(contr)[i])

    rownames(top.table) <- NULL
    de_ls[[i]] <- top.table[, c("Gene","logCPM", "logFC", "p", "adj.p", "Comp")]
  }

  de_df <- do.call("rbind", de_ls)
  de_df <- de_df %>%
    mutate(across(2:5, round, 3)) %>%
    mutate(across(2:5, function(x) ifelse(x<0.001, "<0.0005", x)))

  return(de_df)
}



#' GSEA genes between disconnected clusters
#'
#' @param normcounts A normalized count data matrix: row:genes, column:cells
#' @param OrgDb OrgDb
#'
#' @import clusterProfiler
#' @import dplyr
#' @export

HVGs_GO <- function(normcounts, OrgDb="org.Hs.eg.db") {
  gene.var <- quick_model_gene_var(normcounts)
  HVGs <- rownames(subset(gene.var, FDR < 0.05))

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





#' Find HVGs by using quick scran
#'
#' @param normcounts A normalized count data matrix: row:genes, column:cells
#'
#' @import dplyr
#' @export

HVGs_quick <- function(normcounts) {
  gene.var <- quick_model_gene_var(normcounts)

  df <- as.data.frame(gene.var) %>%
    mutate(across(1:6, round, 3)) %>%
    mutate(across(5:6, function(x) ifelse(x<0.001, "<0.0005", x)))

  return(df)
}
