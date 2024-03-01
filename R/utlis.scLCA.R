#' Useful functions
#'
#' This file contains a few useful functions for scLCA. This is from https://bitbucket.org/scLCA/single_cell_lca/src/master/
#'
# matrix power operator: computes M^power (M must be diagonalizable)
"%^%" <- function(M, power)
  with(eigen(M), vectors %*% (values^power * solve(vectors)))
# cosine distance
cosDist <- function(x) {
  val = (1 - x %*% t(x) / (sqrt(rowSums(x^2) %*% t(rowSums(x^2)))))
  val[val < 0] = 0
  return (val)
}
# cosine distance special
cosDist.part <- function(x, idx1, idx2) {
  idx1 = (1:nrow(x))[idx1]
  idx2 = (1:nrow(x))[idx2]
  u = x[idx1,,drop=F]
  v = x[idx2,,drop=F]
  if(length(idx1) == 1){ u =matrix(u,nrow=1)}
  if(length(idx2) == 1){ v= matrix(v,nrow=1)}
  val = (1 - u %*% t(v) / sqrt((rowSums(u^2)) %*% (t((rowSums(v^2))))))
  val[val < 0] = 0
  return (val)
}
# Tracy Widom test
TWtest <- function(vec) {
  v1 <- sort(vec, decreasing = TRUE)
  sv1 <- v1^2
  m <- length(vec)
  m1 = seq(from=m, to=1)
  pvalue <- c()
  x <- c()
  for (i in 1:m) {
    sumsqrt <- (sum(v1[i:m]))^2
    sum <- sum(sv1[i:m])
    n = ((m1[i]+2)*sumsqrt)/(m1[i]*sum - sumsqrt)
    L = (m1[i]*v1[i])/sum(v1[i:m])
    sigma <- ((1/sqrt(n-1)+1/sqrt(m1[i]))^(1/3))*((sqrt(n-1)+sqrt(m1[i]))/n)
    miu <- ((sqrt(n-1)+sqrt(m1[i]))^2)/n
    x[i] = (L-miu)/sigma
    pvalue[i] <- RMTstat::ptw(x[i], lower.tail = FALSE)
  }
  TWresult <- data.frame(x, pvalue)
  return(TWresult)
}
# sum of squares
sumsq <- function(x) sum(x^2)
# max N elements
maxN <- function(x, N=2){
  len <- length(x)
  sort(x,partial=len-N+1)[len-N+1]
}
# top3
top3 <- function(x, decreasing=T, num = 3){
  sum(sort(x, decreasing=decreasing)[1:num])
}