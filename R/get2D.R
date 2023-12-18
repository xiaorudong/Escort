
#' Generate 2D Embeddings by different methods
#'
#' @param normcounts A normalized count data matrix with genes as rows and cells in columns.
#' @param method options for dimention reduction: MDS, UMAP, PCA, TSNE, ICA
#'
#' @importFrom SCORPIUS reduce_dimensionality
#' @importFrom umap umap
#' @importFrom stats prcomp
#' @importFrom Rtsne Rtsne
#' @importFrom fastICA fastICA
#'
#' @return A data frame containing coordinates of 2-dimensional embedding.
#' @export

getDR_2D <- function(normcounts, method) {

  if (method=="MDS") {
    if (!require("SCORPIUS")) install.packages("SCORPIUS")
    dimred <- SCORPIUS::reduce_dimensionality(t(normcounts), dist="spearman", ndim = 2)
  }

  if(method=="UMAP") {
    if (!require("umap")) install.packages("umap")
    dimred <- umap::umap(t(normcounts))$layout
  }

  if(method=="PCA") {
    if (!require("stats")) install.packages("stats")
    dimred <- stats::prcomp(t(normcounts), rank. = 2)$x
    rownames(dimred) <- rownames(t(normcounts))
  }

  if(method=="TSNE") {
    if (!require("Rtsne")) install.packages("Rtsne")
    dimred <- Rtsne::Rtsne(t(normcounts), dims = 2)$Y
    rownames(dimred) <- rownames(t(normcounts))
  }

  if(method=="ICA") {
    if (!require("fastICA")) install.packages("fastICA")
    dimred <- fastICA::fastICA(t(normcounts), 2)$S
  }
  dimred <- as.data.frame(dimred)
  return(dimred)
}

