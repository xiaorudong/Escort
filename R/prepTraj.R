#' Trajectory Evaluation Object
#'
#' @param dimred A data frame. a 2D embedding from a DR method and the row name is the cell name.
#' @param PT A data frame or a vector. estimated pseudotime from TI. Each column contains the PT for one lineage.
#' @param fitLine the fitted curve line got from TI containing line segments between pairs of points.
#'
#' @return a list containing related information used for future analysis.
#' \itemize{
#'   \item Embedding - A data frame. a 2D embedding from a DR method
#'   \item pse - A data frame. a scaled pseudotime
#'   \item fitLine - A data frame. line segments between pairs of points.
#' }
#' @export

prepTraj <- function(dimred, PT, fitLine) {

  rawpse <- as.data.frame(PT)
  pse <- rawpse / max(rawpse, na.rm=TRUE)
  pse <- pse[match(rownames(dimred), rownames(pse)),]
  pse <- as.data.frame(pse)
  rownames(pse) <- rownames(dimred)
  fitLine <- as.data.frame(fitLine)

  obj <- list(Embedding=dimred, pse=pse, fitLine=fitLine)
  return(obj)
}
