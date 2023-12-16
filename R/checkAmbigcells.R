#' Check the shape of trajectory curves
#'
#' @param obj A basic trajectory evaluation object from \code{prepTraj} function.
#' @param outlierdetect The method used to detect the outlier: "resistant", "adjbox", "asymmetric", "neutral". The default method is "neutral".
#' @param fig A figure showing the ambiguous cells within the embedding.
#'
#' @importFrom FNN get.knn
#' @importFrom scales alpha
#' @importFrom plyr round_any
#' @import RColorBrewer
#' @import graphics
#'
#' @return A list about the results of structure similarity check between high- and low- dimenaional data.
#' \itemize{
#'   \item AmbiguousCells - Vector of ambiguous cells (those considered hard to estimate accurately).
#'   \item NoACells - Number of ambiguous cells.
#'   \item Ambpct - Percentage of ambiguous cells.
#' }
#'
#' @export

UshapeDetector <- function(obj, outlierdetect='neutral', fig=F) {

  dimred <- obj$Embedding
  pse <- obj$pse
  fitLine <- obj$fitLine

  nlineage <- ncol(pse)
  HR_cells_vec <- c()

  if((nrow(dimred)/nlineage)/100<10) {
    num <- 10
  } else {
    num <- plyr::round_any((nrow(dimred)/nlineage)/100, 10, f = round)
  }

  pse_ave <- rowMeans(pse, na.rm = T)

  knn_obj <- FNN::get.knn(dimred, k=num)
  knn_index <- knn_obj$nn.index
  knn_dist <- knn_obj$nn.dist

  outdect <- Escort::boxB(x = as.vector(knn_dist), k = 1.5, method = 'adjbox')
  upperB <- outdect$fences[2]
  filter_ind <- which(knn_dist > upperB, arr.ind = TRUE)
  if(sum(knn_dist > upperB)>0) {
    for (i in 1:nrow(filter_ind)) {
      knn_dist[filter_ind[i, 1], filter_ind[i, 2]] <- NA
      knn_index[filter_ind[i, 1], filter_ind[i, 2]] <- NA
    }
  }

  # check per lineage:
  for (i in 1:nlineage) {
    lineage_pse <- pse[,i]
    lineage_pse[is.na(lineage_pse)] <- pse_ave[is.na(lineage_pse)]
    names(lineage_pse) <- rownames(pse)
    lineage_pse <- lineage_pse[match(rownames(dimred), names(lineage_pse))]

    knn_pt <- matrix(apply(knn_index, 1, function(x) lineage_pse[x]), ncol = num, byrow = T)
    knn_sd <- apply(knn_pt, 1, function(x) sd(x, na.rm = T))

    outdect0 <- Escort::boxB(x = knn_sd, k = 1.5, method="resistant")
    upperB0 <- outdect0$fences[2]

    outdect1 <- Escort::boxB(x = knn_sd, k = 1.5, method="adjbox")
    upperB1 <- outdect1$fences[2]

    outdect2 <- Escort::boxB(x = knn_sd, k = 1.5, method="asymmetric")
    upperB2 <- outdect2$fences[2]

    upperB12 <- median(knn_sd[knn_sd<=max(upperB1, upperB2) & knn_sd>=min(upperB1, upperB2)], na.rm = T)

    if(outlierdetect=="resistant") upperB <- upperB0
    if(outlierdetect=="adjbox") upperB <- upperB1
    if(outlierdetect=="asymmetric") upperB <- upperB2
    if(outlierdetect=="neutral") upperB <- upperB12

    if (any(knn_sd>upperB, na.rm = T)) HR <- which(knn_sd>upperB) else HR <- NULL
    rm_cells <- as.numeric(unique(filter_ind[,1]))

    HR <- HR[! HR %in% rm_cells]
    HR_cells <- rownames(dimred)[HR]
    HR_cells_vec <- c(HR_cells_vec, HR_cells)
  }
  allHR_cells <- unique(HR_cells_vec)
  allHR.no <- length(unique(HR_cells_vec))

  if(fig) {
    colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
    pse$Ave <- rowMeans(pse, na.rm = T)
    plotcol <- colors[cut(pse$Ave, breaks=100)]
    plot(dimred, col = scales::alpha(plotcol, 0.7), pch=16)
    segments(x0 = fitLine$x0,
             y0 = fitLine$y0,
             x1 = fitLine$x1,
             y1 = fitLine$y1, lwd = 3)
    points(dimred[allHR_cells,1], dimred[allHR_cells,2], col = "black", pch=16)
  }

  allHR_rate <- allHR.no/nrow(dimred)
  return(list(AmbiguousCells=allHR_cells, NoACells=allHR.no, Ambpct=allHR_rate))
}
