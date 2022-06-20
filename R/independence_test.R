#' Independence Test
#'
#' This function carries out hypothesis test for independence between disjoint groups of covariates.
#' @param X a matrix or dataframe, the observation data.
#' @param part a vector, how covariates are divided in order.
#' @examples
#' # test whether the three variables of sweat dataset are independent.
#' independence_test(sweat,c(1,1,1))
#' @export
independence_test <- function(X,part){
  X <- as.matrix(X)
  p <- ncol(X)
  n <- nrow(X)
  A <- cov(X)*(n-1)
  ind <- c(0,cumsum(part))
  miditem <- 0
  k <- length(part)
  for (i in 1:k) {
    start <- ind[i]+1
    end <- ind[i+1]
    cache <- A[start:end,start:end]
    if(!is.matrix(cache)){
      miditem <- miditem + log(cache)
    }else{
      miditem <- miditem + log(det(cache))
    }
  }
  llambda <- n/2*(log(det(A))-miditem)
  b <- n-3/2-(p^3-sum(part^3))/(3*(p^2-sum(part^2)))
  f <- 1/2*(p*(p+1)-sum(part*(part+1)))
  statistic <- -b*2/n*llambda
  pvalue <- 1-pchisq(statistic, df = f)
  cat("Independence test.\n")
  cat("chisq-statistic = ", statistic, "df = ", f, "p-value = ", pvalue, "\n")
  cat("alternative hypothesis: the divided covariate groups are dependent","\n")
  invisible(list(statistic = statistic, pvalue = pvalue))
}
