#' Trajectory Evaluation Object
#'
#' @param dimred A data frame. A 2D embedding from a DR method and the row name is the cell name.
#' @param PT A data frame or a vector of estimated pseudotime from TI. Each column contains the pseudotime for one lineage.
#' @param fitLine the fitted curve line from trajectory containing line segments between pairs of points.
#'
#' @return a list containing related information used for future analysis.
#' \itemize{
#'   \item Embedding - A data frame. A 2D embedding from a DR method.
#'   \item pse - A data frame. Scaled pseudotime.
#'   \item fitLine - A data frame. Contains line segments between pairs of points.
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

#' Format line segments (when using Slingshot)
#'
#' @param fittedLines principle curves obtained form the slingCurves() function in Slingshot
#'
#' @export

segFormat <- function(fittedLines) {

  segLines <- do.call(rbind, lapply(fittedLines, function(x) {
    df_seg <- cbind(x[-nrow(x),],x[-1,])
    colnames(df_seg) <- c("x0", "y0", "x1", "y1")
    return(df_seg)
  }))
  segLines <- as.data.frame(segLines)

  return(segLines)
}
