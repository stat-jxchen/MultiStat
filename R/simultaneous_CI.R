#' Simultaneous Confidence Interval
#'
#' This function estimates simultaneous confidence interval for the mean of a
#' multivariate normal distribution with Bonferroni method and Scheffe method.
#' @param X a matrix or dataframe, the observation data.
#' @param method string, "B" for Bonferroni method and "S" for Scheffe method.
#' @param alpha number, the confidence level.
#' @export
simultaneous_CI <- function(X, method = "B",alpha = 0.05){
  Xbar <- apply(X,2,mean)
  n <- nrow(X)
  p <- ncol(X)
  S <- cov(X)
  re <- matrix(0,nrow = p, ncol = 3)
  re[,1] <- Xbar
  if(method == "B"){
    margin <- -qt(alpha/(2*p),df = n-1)*sqrt(diag(S))/n
  }else if(method == "S"){
    margin <- sqrt((n-1)*p/(n-p)*qf(alpha,df1 = p, df2 = n-p))*sqrt(diag(S))/n
  }
  re[,2] <- Xbar-margin
  re[,3] <- Xbar+margin
  colnames(re) <- c("mean","lower","upper")
  return(re)
}
