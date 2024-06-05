
#' Find top 30 DE genes between disconnected clusters
#'
#' @param rawcounts A raw count data matrix: row:genes, column:cells
#' @param cls Identified disconnected clusters
#'
#' @importFrom limma makeContrasts 
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmFit glmLRT topTags
#' @importFrom shinyjs show
#' @importFrom stats lag filter
#' @import magrittr
#' @importFrom dplyr mutate across
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
    mutate(across(2:5, function(x) ifelse(x<0.001, "<0.001", x)))

  return(de_df)
}



#' GSEA genes between disconnected clusters
#'
#' @param normcounts A normalized count data matrix: row:genes, column:cells
#' @param OrgDb OrgDb
#'
#' @importFrom  stats lag filter
#' @import magrittr
#' @importFrom dplyr mutate across
#' @import clusterProfiler
#' @export

HVGs_GO <- function(normcounts, OrgDb="org.Hs.eg.db", geneVar=NULL) {
  if(is.null(geneVar)) {
    geneVar <- Escort::quick_model_gene_var(normcounts)
  } 
  HVGs <- rownames(subset(geneVar, FDR < 0.05))

  ego <- clusterProfiler::enrichGO(gene = HVGs, OrgDb = OrgDb, ont = "ALL", keyType = "SYMBOL")
  if (is.null(ego)) {
    return(NULL)
  }

  if (!is.null(ego)){
    df <- ego@result[, c(1, 3 , 4, 6, 7)]
    head(df)
    df <- as.data.frame(df) %>%
      dplyr::mutate(across(4:5, round, 3)) %>%
      dplyr::mutate(across(4:5, function(x) ifelse(x<0.001, "<0.001", x)))
    head(df)
    colnames(df)[1] <- "OntologyType"
    colnames(df)[4] <- "Pval"
    colnames(df)[5] <- "AdjPval"

    return(df)
  } 
}





#' Find HVGs by using quick scran
#'
#' @param normcounts A normalized count data matrix: row:genes, column:cells
#'
#' @importFrom  stats lag filter
#' @import magrittr
#' @importFrom dplyr mutate across
#' @export

HVGs_quick <- function(normcounts) {
  gene.var <- Escort::quick_model_gene_var(normcounts)

  df <- as.data.frame(gene.var)
  df <- df[order(df$FDR, df$mean),]
  df <- df %>% dplyr::mutate(across(1:6, round, 3)) %>%
    dplyr::mutate(across(5:6, function(x) ifelse(x<0.001, "<0.001", x)))

  return(df)
}
