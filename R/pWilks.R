#' Wilks Distribution
#'
#' Calculate the value or approximate value of the distribution function of the Wilks distribution.
#' @param q vector of quantiles.
#' @param p the number of order.
#' @param n1 the first degree of freedom.
#' @param n2 the second degree of freedom.
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x) otherwise, P(X > x).
#' @examples 1-pWilks(0.6621,4,57,2)
#' @export
pWilks <- function(q, p, n1, n2, lower.tail = TRUE){
  if(n2 == 1){
    f <- (n1+1-p)/p*(1-q)/q
    return(pf(f, df1 = p, df2 = n1+1-p, lower.tail = lower.tail))
  }
  if(n2 == 2){
    f <- (n1+1-p)/p*(1-sqrt(q))/sqrt(q)
    return(pf(f, df1 = 2*p, df2 = 2*(n1+1-p), lower.tail = lower.tail))
  }
  if(p == 1){
    f <- n1/n2*(1-q)/q
    return(pf(f, df1 = n2, df2 = n1, lower.tail = lower.tail))
  }
  if(p == 2){
    f <- (n1-1)/n2*(1-sqrt(q))/sqrt(q)
    return(pf(f, df1 = 2*n2, df2 = 2*(n1-1), lower.tail = lower.tail))
  }
  chisq <- ((p-n2+1)/2-n1)*log(q)
  return(pchisq(chisq, df = p*n2, lower.tail = lower.tail))
}
