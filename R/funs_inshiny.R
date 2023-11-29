
#' Find DE genes between disconnected clusters
#'
#' @param rawcounts A raw count data matrix: row:genes, column:cells
#' @param cls Identified disconnected clusters
#'
#' @import limma
#' @import edgeR
#' @import dplyr
#' @export

DE_limma <- function(rawcounts, cls) {

  a <- match(colnames(rawcounts), names(cls))
  cls <- paste("Cluster", cls[a], sep = "")
  clusters <- as.factor(cls)
  combn_tb <- t(combn(levels(clusters), 2))

  d <- DGEList(rawcounts)
  d <- calcNormFactors(d)
  mm <- model.matrix(~0+clusters)
  colnames(mm) <- gsub("clusters", "", colnames(mm))
  y <- voom(d, mm)
  fit <- lmFit(y, mm)
  comp <- apply(combn_tb, 1, function(x) paste(x,collapse = "-"))
  contr <- makeContrasts(contrasts=comp, levels=mm)
  fit <- contrasts.fit(fit, contr)
  res <- eBayes(fit)
  de_ls <- list()
  for (i in 1:ncol(fit$coefficients)) {
    top.table <- topTable(res, adjust="BH",coef=1, p.value=0.05, n=Inf)
    top.table <- as.data.frame(top.table)
    top.table$Comp <- colnames(fit$coefficients)[i]
    de_ls[[i]] <- top.table[, c("AveExpr", "logFC", "P.Value", "adj.P.Val", "Comp")]
  }

  de_df <- do.call("rbind", de_ls)
  de_df <- de_df %>%
    mutate(across(1:4, round, 3)) %>%
    mutate(across(1:4, function(x) ifelse(x<0.001, "<0.0005", x)))

  return(de_df)
}



#' GSEA genes between disconnected clusters
#'
#' @param norm_counts A normalized count data matrix: row:genes, column:cells
#' @param OrgDb OrgDb
#'
#' @import clusterProfiler
#' @import dplyr
#' @export

HVGs_GO <- function(norm_counts, OrgDb="org.Hs.eg.db") {
  gene.var <- quick_model_gene_var(norm_counts)
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
#' @param norm_counts A normalized count data matrix: row:genes, column:cells
#'
#' @import dplyr
#' @export

HVGs_quick <- function(norm_counts) {
  gene.var <- quick_model_gene_var(norm_counts)

  df <- as.data.frame(gene.var) %>%
    mutate(across(1:6, round, 3)) %>%
    mutate(across(5:6, function(x) ifelse(x<0.001, "<0.0005", x)))

  return(df)
}
