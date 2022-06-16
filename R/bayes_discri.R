#' Bayes Discrimination
#'
#' This function do Bayes discrimination.
#'
#' @param x a vector, a single sample data to be predicted.
#' @param prior a vector, prior probability of each class.
#' @param loss a matrix, the ijth element denotes the loss of misjudging
#' class j as class i. If the loss is equal, it could be set to a number.
#' @param pdf a function, the common form of pdf of each class.
#' @param para_list a list, each element is an additional single parameter of the pdf.
#' @export
bayes_discri <- function(x, prior, loss, pdf, para_list){
  # calculate the density of each class at x according to the pdf
  nclass <- length(prior)
  f <- numeric(nclass)
  for (i in 1:nclass) {
    f[i] <- pdf(x,para_list[i])
  }
  # when loss functions are the same
  if(dim(as.matrix(loss)) == c(1,1)){
    ones <- matrix(1,nclass,nclass)
    loss <- ones-diag(1,nclass,nclass)
  }
  mean_loss <- numeric(nclass)
  for (i in 1:nclass) {
    mean_loss[i] <- prior*loss[i,]*f
  }
  re <- which.min(mean_loss)
  return(re)
}
