#' paste the boxB() function from univOutl package
#' @import robustbase
#' @importFrom dplyr src summarize
#' @importFrom shinyjs html
#' 
#' @export

boxB <- function (x, k = 1.5, method = "asymmetric", weights = NULL,
          id = NULL, exclude = NA, logt = FALSE)
{
  if (is.null(id))
    id <- 1:length(x)
  tst <- x %in% exclude
  if (sum(tst) == 0) {
    to.check <- integer(0)
    yy <- x
    lab <- id
    ww <- weights
  }
  else {
    to.check <- id[tst]
    yy <- x[!tst]
    lab <- id[!tst]
    ww <- weights[!tst]
  }
  if (logt) {
    yy <- log(yy + 1)
    warning("Please note that log(x+1) transforation is considered")
  }
  if (is.null(weights))
    qq <- quantile(x = yy, probs = c(0.25, 0.5, 0.75))
  else qq <- Hmisc::wtd.quantile(x = yy, weights = ww, probs = c(0.25,
                                                                 0.5, 0.75))
  if (method == "resistant" | method == "boxplot") {
    vql <- vqu <- qq[3] - qq[1]
    if (abs(vql) < 1e-06)
      warning("IQR is 0")
    low.b <- qq[1] - k * vql
    up.b <- qq[3] + k * vqu
  }
  if (method == "asymmetric") {
    vql <- qq[2] - qq[1]
    if (abs(vql) < 1e-06)
      warning("(Q2-Q1) = 0")
    vqu <- qq[3] - qq[2]
    if (abs(vqu) < 1e-06)
      warning("(Q3-Q2) = 0")
    low.b <- qq[1] - 2 * k * vql
    up.b <- qq[3] + 2 * k * vqu
  }
  if (method == "adjbox") {
    # warning("With method='adjbox' the argument k is set equal to 1.5")
    ck <- inherits(try(robustbase::mc(yy), silent = TRUE),
                   "try-error")
    if (ck)
      medc <- robustbase::mc(yy, doScale = TRUE)
    else medc <- robustbase::mc(yy, doScale = FALSE)
    # message("The MedCouple skewness measure is: ", round(medc,
    #                                                      4))
    if (is.null(weights)) {
      if (ck)
        aa <- robustbase::adjboxStats(x = yy, doScale = TRUE)
      else aa <- robustbase::adjboxStats(x = yy, doScale = FALSE)
      low.b <- aa$fence[1]
      up.b <- aa$fence[2]
    }
    else {
      if (medc >= 0) {
        low.b <- qq[1] - 1.5 * exp(-4 * medc) * (qq[3] -
                                                   qq[1])
        up.b <- qq[3] + 1.5 * exp(3 * medc) * (qq[3] -
                                                 qq[1])
      }
      else {
        low.b <- qq[1] - 1.5 * exp(-3 * medc) * (qq[3] -
                                                   qq[1])
        up.b <- qq[3] + 1.5 * exp(4 * medc) * (qq[3] -
                                                 qq[1])
      }
    }
  }
  names(low.b) <- "low"
  names(up.b) <- "up"
  outl <- (yy < low.b) | (yy > up.b)
  lower <- (yy < low.b)
  upper <- (yy > up.b)
  # if (sum(outl) == 0)
  #   message("No outliers found")
  # else {
  #   message("No. of outliers in left tail: ", sum(yy < low.b))
  #   message("No. of outliers in right tail: ", sum(yy >
  #                                                    up.b), "\n")
  # }
  fences <- c(low.b, up.b)
  names(fences) <- c("lower", "upper")
  if (sum(outl) == 0) {
    fine <- list(quartiles = qq, fences = fences, excluded = to.check,
                 outliers = integer(0))
  }
  else {
    fine <- list(quartiles = qq, fences = fences, excluded = to.check,
                 outliers = lab[outl], lowOutl = lab[lower], upOutl = lab[upper])
  }
  fine
}
