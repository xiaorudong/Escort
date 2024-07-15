#' Test Homogeneous
#'
#' @param normcounts A normalized count data matrix: row:genes, column:cells
#' @param HVGs Highly variable genes. If there is no input of the highly variable genes, they will be calcualted by `modelGeneVar()` function.
#' @param g The number of highly variable genes to test
#' @param num.sim The number of simulations generated for permutation test on spearman correlation. Default is 20000. Reduce this to improve speed.
#'
#' @importFrom jmuOutlier perm.cor.test
#' @importFrom dplyr first setequal intersect union setdiff
#' @importFrom stats cor cov start window complete.cases end var xtabs aggregate na.omit sd na.exclude
#' @importFrom clusterProfiler rename
#' 
#' @return A list about the results of trajectory signal detection
#' \itemize{
#'   \item signal_pct - the percentage of trajectory signal detected. If the percentage exceeds 50%, we conclude that the trajectory signal detected.
#'   \item decision - the decision sentence.
#' }
#' @export

step1_testHomogeneous <- function(normcounts, HVGs=NULL, g=100, num.sim = 20000) {
  if(is.null(HVGs)) {
    gene.var <- quick_model_gene_var(normcounts)
    HVGs <- rownames(gene.var)[1:g]
  }
  if(length(HVGs)<20) return("Please input more highly variable genes.")

  pcDat <- prcomp(t(normcounts))
  pca_cells <- pcDat$x
  est_pt <- sort(pca_cells[,1])

  time_df <- data.frame(Cell=names(est_pt), PT=est_pt)
  time_df$Cell <- factor(time_df$Cell, levels=colnames(normcounts))
  time_df <- time_df[order(time_df$Cell),]
  comb_df <- cbind(time_df, t(normcounts[HVGs,]))
  comb_df <- comb_df[order(comb_df$PT),]

  withProgress(message = "Testing for homogeneous cells", value = 0, {
    incProgress(1/2)
    p_vec_perm <- apply(comb_df[,3:ncol(comb_df)], 2,
                      function(x) jmuOutlier::perm.cor.test(comb_df[,2], x, method = "spearman", num.sim = num.sim)$p.value)
    incProgress(1)
  })
  padj <- p.adjust(p_vec_perm, method = "fdr")
  # res_perm_p <- mean(p_vec_perm<0.05)
  res_perm_padj <- mean(padj<0.05)
  decision <- ifelse(res_perm_padj>=0.46,
                     "The trajectory signal is detected.",
                     "No trajectory signal is detected.")

  return(list(signal_pct = res_perm_padj, decision=decision))
}

