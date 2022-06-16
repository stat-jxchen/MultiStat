#' Distance Discrimination
#'
#' This function do distance discrimination with Mahalanobis distance.
#'
#' @param train_data a list containing the training data of each class.
#' If it is \code{NULL}, then \code{mu} and \code{Sig} must be given.
#' @param x a vector, a single sample data to be predicted.
#' @param mu a list containing true mean vector of each class. If \code{train_data}
#' has been given, it will be ignored.
#' @param Sig a list containing true covariance matrix of each class. If set
#' \code{Sig_equal} as \code{TRUE}, it just need to be the common covariance matrix.
#' If \code{train_data} has been given, it will be ignored.
#' @param Sig_equal logical, which indicates whether the covariance matrix of each
#' class is equal.
#' @export
dist_discri <- function(train_data = NULL, x, mu, Sig, Sig_equal = FALSE){
  # known mu, sig
  if(is.null(train_data)){
    if(missing(mu)|missing(Sig)){
      stop("Mean and covariance matrix must be specified when there
           is no training set !")
    }else{
      group_num <- length(mu)
      if(Sig_equal){
        Sig <- rep(list(Sig),group_num)
      }
    }
  }else{
      # unknown mu,sig
      group_num <- length(train_data)
      mu <- train_data %>% map(colMeans)
      Sig <- train_data %>% map(cov)
      if(Sig_equal){
        nj <- train_data %>% map_dbl(nrow)
        n <- sum(nj)
        Sig_com <- map2(Sig,(nj-1),`*`) %>% Reduce(`+`,.) %>% `/`(n-group_num)
        Sig <- rep(list(Sig_com),group_num)
      }
    }
  # Calculate the Mahalanobis distance from the test point to each center
  re <- map2_dbl(mu,Sig,mahalanobis,x = x) %>% which.min()
  return(re)
}
