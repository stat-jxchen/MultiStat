#' Hotelling Distribution
#'
#' Calculate the value of the distribution function of the Hotelling distribution.
#' @param q vector of quantiles.
#' @param p the number of order.
#' @param n the degree of freedom.
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x) otherwise, P(X > x).
#' @examples pHotelling(2,2,2)
#' @export
pHotelling <- function(q, p, n, lower.tail = TRUE){
  f <- (n-p+1)/(n*p)*q
  return(pf(f, df1 = p, df2 = n-p+1, lower.tail = lower.tail))
}

