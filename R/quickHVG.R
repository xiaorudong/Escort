#' Quick HVG (simplified from scran)
#'
#' @param normcounts A normalized count data matrix: row:genes, column:cells
#' @importFrom limma weightedLowess weighted.median
#'
#' @export


quick_model_gene_var <- function(normcounts) {
  
  means <- rowMeans(normcounts)
  RowVar <- function(x) {
    rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
  }
  vars <- RowVar(normcounts)
  names(means) <- names(vars) <- rownames(normcounts)
  
  
  # collected <- scran:::.decompose_log_exprs(x.stats$means, x.stats$vars, fit.stats$means, fit.stats$vars, 
                                    # x.stats$ncells)
  # Filtering out zero-variance and low-abundance genes.
  is.okay <- !is.na(vars) & vars > 1e-8 & means >= 0.1 
  v <- vars[is.okay]
  m <- means[is.okay]
  
  out <- density(m, adjust=1, from=min(m), to=max(m))
  w <- 1/approx(out$x, out$y, xout=m)$y 
  w <- w/mean(w)
  
  to.fit <- log(v)
  left.edge <- min(m)
  PARAMFUN <- function(x) { pmin(1, x/left.edge) } 
  
  nls.args <- list()
  # nls.args <- scran:::.setup_nls_args(nls.args, start.args=list(vars=v, means=m))
  
  start.args=list(vars=v, means=m)
  
  nls.call <- do.call(call, c("nls", nls.args))
  nls.call <- match.call(nls, nls.call)
  nls.args <- as.list(nls.call)[-1]
  
  control <- nls.control(warnOnly=TRUE, maxiter=500)
  
  o <- order(m)
  n <- length(v)
  # Estimating the gradient from the left.
  left.n <- min(100, n*0.1)
  keep <- head(o, max(1, left.n))
  y <- v[keep]
  x <- m[keep]
  grad <- coef(lm(y~0+x))
  
  # Two-dimensional grid search is the most reliable way of estimating the remaining parameters.
  b.grid.pts <- 2^seq(from=-5, to=5, length.out=10)
  n.grid.pts <- 2^seq(from=0, to=10, length.out=10)
  hits <- expand.grid(B=b.grid.pts, n=n.grid.pts)
  
  grid.ss <- mapply(B=hits$B, n=hits$n, FUN=function(B, n) {
    resid <- v - (grad*B*m)/(m^n + B)
    sum(resid^2)
  })
  
  chosen <- which.min(grid.ss)
  N <- hits$n[chosen]
  B <- hits$B[chosen]
  A <- B * grad
  
  raw.start <- list(n=N, b=B, a=A)
    
  start <- list(A=log(raw.start$a),
                B=log(raw.start$b),
                N=log(pmax(1e-8, raw.start$n-1))) # reflects enforced positivity in formula.
  
  altogether <- c(nls.args, list(control=control, start=start))
  keep <- !duplicated(names(altogether)) # nls.args are favoured.
  nls.args <- altogether[keep]
  
  
  nls.args$formula <- v ~ (exp(A)*m)/(m^(1+exp(N)) + exp(B))
  nls.args$weights <- w
  nls.args$control$warnOnly <- FALSE
  init.fit <- try(do.call(nls, nls.args), silent=TRUE) 
  if (is(init.fit, "try-error")) {
      Aest <- exp(nls.args$start$A)
      Best <- exp(nls.args$start$B)
      Nest <- exp(nls.args$start$N)+1
      PARAMFUN <- function(x) { Aest * x / (x^Nest + Best) }
      to.fit <- to.fit - log(PARAMFUN(m))
    } else {
      to.fit <- to.fit - log(fitted(init.fit))
      PARAMFUN <- function(x) { predict(init.fit, data.frame(m=x)) }
    }

  lfit <- limma::weightedLowess(m, to.fit, weights=w)
  LOESSFUN <- approxfun(m, lfit$fitted, rule=2)
  
  UNSCALEDFUN <- function(x) { 
    exp(LOESSFUN(x)) * PARAMFUN(x)
  }

  # Adjusting for any scale shift due to fitting to the log-values.
  leftovers <- v/UNSCALEDFUN(m)
  med <- limma::weighted.median(leftovers, w, na.rm=TRUE)
  
  OUT <- function(x) { 
    output <- UNSCALEDFUN(x) * med
    names(output) <- names(x)
    output
  }
  std.dev <- unname(limma::weighted.median(abs(leftovers/med - 1), w, na.rm=TRUE)) * 1.4826 
  
  fit <- list(trend=OUT, std.dev=std.dev)
  
  output <- data.frame(mean=means, total=vars, tech=fit$trend(means))
  output$bio <- output$total - output$tech
  output$p.value <- pnorm(output$bio/output$tech, sd=fit$std.dev, lower.tail=FALSE)
  output$FDR <- p.adjust(output$p.value, method="BH")

  output <- output[order(output$bio, decreasing = TRUE),]
  
  return(output)
}
