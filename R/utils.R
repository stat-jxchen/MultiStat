#' Drop population-distinguishing variables
#' @param mydata a tibble or matrix or dataframe, which contains a variable to show the
#' population.
#' @param label a character, the name of the variable to show the population.
#' @export
#' @noRd
drop_label <- function(mydata, label){
  mydata <- select(mydata,-c(label))
  return(mydata)
}

#' Calculate the trace of a matrix
#' @param x a matrix
#' @export
#' @noRd
tr <- function(x){
  return(sum(diag(x)))
}

#' Calculate the inverse of matrx
#' @param x a matrix
#' @export
#' @noRd
mysolve <- function(x){
  if(rcond(x)<=.Machine$double.eps){
    if(det(x) != 0){
      warning("Condition number of x is too small,
            will use the eigenvalue to find the inverse,
            the result is not necessarily reliable")
      eig_re <- eigen(x)
      x_inv <- eig_re$vectors%*%diag(1/(eig_re$values),dim(x))%*%t(eig_re$vectors)
      return(x_inv)
    }else{
      stop("x is not invertible")
    }
  }else{
    return(solve(x))
  }
}
