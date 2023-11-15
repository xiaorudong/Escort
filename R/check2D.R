

#' Structure Similarity between Orignal data and Embedding
#'
#' @param Cluters object from \code{DCClusterscheck} fucntion.
#' @param norm_counts a normalized count data matrix: row:genes, column:cells
#' @param dimred A data frame. a 2D embedding from a DR method and the row name is the cell name.
#'
#' @importFrom FNN get.knn
#' @importFrom scales alpha
#'
#' @return A list about the results of structure similarity check between high- and low- dimenaional data.
#' \itemize{
#'   \item SimilarityLevel - Same group level
#'   \item GoodRate - High same group level percentage
#'   \item all - Same group level for each cell
#'   \item cls - clusters
#' }
#'
#' @export

Similaritycheck <- function(Cluters, norm_counts, dimred) {

  dimred <- as.data.frame(dimred)
  dimred <- dimred[colnames(norm_counts),]
  cut_avg <- Cluters$Clusters
  cut_avg <- cut_avg[colnames(norm_counts)]
  K <- Cluters$K

  knn_index <- FNN::get.knn(dimred, k=3)$nn.index
  knn_group <- matrix(apply(knn_index, 1, function(x) cut_avg[x]), ncol = 3, byrow = T)
  knn_df <- cbind(cut_avg, knn_group)
  knn_overlap <- apply(knn_df, 1, function(x) sum(x[2:4]!=x[1]))
  t_knn <- table(knn_overlap)

  good_rate <- sum(t_knn[as.numeric(names(t_knn))<=1])/ncol(norm_counts)

  plotcol <- as.numeric(as.factor(cut_avg))

  plot(dimred, col = scales::alpha(plotcol,0.7), pch=16)

  return(list(SimilarityLevel=t_knn, GoodRate=good_rate, all=knn_overlap, cls=cut_avg))
}



#' Goodness of Fit
#'
#' @param dimred A data frame. a 2D embedding from a DR method
#' @param alpha alpha value to calculate the alpha-convex hull. If left NULL, the alpha will be calculated based on the  first derivative.
#'
#' @import alphahull
#' @importFrom shotGroups getMinCircle
#' @importFrom rstatix mahalanobis_distance
#'
#' @return A list about the results of the GOF
#' \itemize{
#'   \item cellArea - The area of the alpha-convex hull of a sample of points.
#'   \item alpha - alpha value used to calculate the alpha-convex hull
#'   \item alphahull_obj - 	Object of class "ahull"
#'   \item occupiedRate - the proportion of alpha-convex hull area to the minimum circle that could enclose all cells in the embedding
#' }
#' @export

GOFeval <- function(dimred, alpha=NULL) {

  dimred <- as.data.frame(dimred)

  # remove the outliers:
  outlier_md <- rstatix::mahalanobis_distance(dimred)
  outliers <- rownames(outlier_md)[outlier_md$is.outlier]

  clean_dimred <- dimred
  if(length(outliers)>0) clean_dimred <- dimred[!rownames(dimred) %in% outliers,]

  minCircle <- shotGroups::getMinCircle(as.matrix(clean_dimred))
  circle_rad <- minCircle$rad
  circlearea <- pi*minCircle$rad^2

  if(is.null(alpha)) {

    upper_bound <- round(circle_rad+0.01, 2)

    if(circle_rad<1) {
      seq2 <- seq(0, upper_bound, by=0.05)
    } else {
      seq2 <- seq(0, upper_bound, length.out=15)
    }
    area <- c()
    for (a in seq2) {
      outline_dimred <- alphahull::ahull(clean_dimred, alpha=a)
      area <- c(area, alphahull::areaahull(outline_dimred))
    }
    # norm_area <- (area-min(area))/(max(area)-min(area))

    data <- data.frame(x=seq2, y=area)
    data <- data[order(data$x, decreasing = T),]
    # plot(data[,1:2])
    data$slope = c(NA, diff(data$y)/diff(data$x))
    data$slope_chg = c(NA, round(diff(data$slope),5))
    data$change = ifelse(data$slope_chg != 0, "change","")

    big_slope <- min(head(order(data$slope_chg), 2))
    alpha=data$x[big_slope]

    outline_dimred <- alphahull::ahull(clean_dimred, alpha=alpha)
    area <- alphahull::areaahull(outline_dimred)
    plot(outline_dimred)
    points(dimred[rownames(dimred) %in% outliers,1], dimred[rownames(dimred) %in% outliers,2])
  } else {
    outline_dimred <- alphahull::ahull(clean_dimred, alpha=alpha)
    area <- alphahull::areaahull(outline_dimred)
    plot(outline_dimred)
    points(dimred[rownames(dimred) %in% outliers,1], dimred[rownames(dimred) %in% outliers,2])
  }

  op_rate <- round(area/circlearea, 3)

  return(list(cellArea=area, alpha=alpha,
              alphahull_obj = outline_dimred,
              occupiedRate=op_rate))
}


