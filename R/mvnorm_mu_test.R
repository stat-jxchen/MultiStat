#' Single or Bivariate Population Mean Test
#'
#' This function carries out hypothesis test for the means of one or two multivariate normal distributions.
#'
#' @param X a matrix or a data frame which denotes a sample.
#' @param Y a matrix or a data frame which denotes the second sample.
#' @param mu0 a vector, the value of mu under null hypothesis.
#' @param Sig0 a matrix or \code{unequal}(only when \code{Y} is not \code{NULL}),
#' when it is a matrix, it will be used in test and means covariance matrices are
#' equal and known, when it is \code{unequal}, it means the covariance matrix of two
#' population are not equal, and when it is \code{NULL} means covariance matrices
#' are equal but unknown.
#' @note When \code{Y} is \code{NULL}, it will carry out one-sample test. If \code{Y}
#' is not \code{NULL}, \code{mu0} can be \code{NULL}.
#' @examples
#' # One-Sample Mean Test with Unknown Covariance Matrix
#' mu0 <- c(4, 50, 10)
#' mvnorm_mu_test(X = sweat, mu0 = mu0)
#' # One-Sample Mean Test with known Covariance Matrix
#' sig0 <- cov(sweat)
#' mvnorm_mu_test(sweat,mu0 = mu0, Sig0 = sig0)
#' @export
mvnorm_mu_test <- function(X, Y = NULL, mu0 = NULL, Sig0 = NULL){
  n <- nrow(X)
  p <- ncol(X)
  Xbar <- apply(X, 2, mean)
  if(is.null(Y)){
    if(is.null(mu0)){
      stop("One-sample test requires the mean given the null hypothesis!")
    }
    S <- cov(X)
    if(is.null(Sig0)){
      # One-Sample Mean Test with Unknown Covariance Matrix
      t2 <- n*(Xbar-mu0)%*%solve(S,Xbar-mu0)
      statistic <- (n-p)/((n-1)*p)*t2
      pvalue <- 1-pf(statistic, df1 = p, df2 = n-p)
      cat("One-sample mean test when covariance matrix is unknown.\n")
      cat("F-statistic = ", statistic, "df1 = ", p, "df2 = ", n-p,
          "p-value = ", pvalue, "\n")
      cat("alternative hypothesis: true mean is not equal to mu0","\n")
    }else{
      # One-Sample Mean Test with Known Covariance Matrix
      statistic <- n*(Xbar-mu0)%*%solve(Sig0, Xbar-mu0)
      pvalue <- 1-pchisq(statistic, df = p)
      cat("One-sample mean test when covariance matrix is known.\n")
      cat("chisq-statistic = ", statistic, "df = ", p, "p-value = ", pvalue, "\n")
      cat("alternative hypothesis: true mean is not equal to mu0","\n")
    }
    invisible(list(statistic = statistic, pvalue = pvalue,
                   estimate = list(Xbar = Xbar, S = S)))
  }else{
    if(ncol(Y) != p){
      stop("The number of columns of X and Y must be equal!")
    }
    S1 <- cov(X)
    S2 <- cov(Y)
    Ybar <- apply(Y, 2, mean)
    m <- nrow(Y)
    if(is.null(Sig0)){
      # Two-Sample Means Test When Covariance Matrices Are Equal But Unknown
      A1 <- S1*(n-1)
      A2 <- S2*(m-1)
      t2 <- n*m*(n+m-2)/(n+m)*(Xbar-Ybar)%*%solve(A1+A2, Xbar-Ybar)
      statistic <- (n+m-2-p+1)/((n+m-2)*p)*t2
      pvalue <- 1-pf(statistic, df1 = p, df2 = n+m-p-1)
      cat("Two-sample mean test when covariance matrix is equal but unknown.\n")
      cat("F-statistic = ", statistic, "df1 = ", p, "df = ",n+m-p-1,
          "p-value = ", pvalue, "\n")
      cat("alternative hypothesis: true mean of two population is not equal","\n")
      invisible(list(statistic = statistic, pvalue = pvalue,
                     estimate = list(Xbar = Xbar, Ybar = Ybar, S1 = S1, S2 = S2)))
    }else if(is.matrix(Sig0)){
      # Two-Sample Means Tests When Covariance Matrices Are Equal and Known
      statistic <- (n*m/(n+m))*(Xbar-Ybar)%*%solve(Sig0, Xbar-Ybar)
      pvalue <- 1-pchisq(statistic, df = p)
      cat("Two-sample mean test when covariance matrix is equal and known.\n")
      cat("chisq-statistic = ", statistic, "df = ", p, "p-value = ", pvalue, "\n")
      cat("alternative hypothesis: true mean of two population is not equal","\n")
      invisible(list(statistic = statistic, pvalue = pvalue,
                     estimate = list(Xbar = Xbar, Ybar = Ybar, S1 = S1, S2 = S2)))
    }else if(Sig0 == "unequal"){
      if(n == m){
        # Two-Sample Means Tests with Unequal and Unknown Covariance Matrix but Equal Sample Sizes
        Z <- X - Y
        Sz <- cov(Z)
        Zbar <- apply(Z, 2, mean)
        t2 <- n*Zbar%*%solve(Sz,Zbar)
        statistic <- (n-p)/((n-1)*p)*t2
        pvalue <- 1-pf(statistic, df1 = p, df2 = n-p)
        cat("Two-sample mean test when covariance matrix is unequal
            and unknown with equal size.\n")
        cat("F-statistic = ", statistic, "df1 = ", p, "df2 = ", n-p,
            "p-value = ", pvalue, "\n")
        cat("alternative hypothesis: true mean of two population is not equal","\n")
        invisible(list(statistic = statistic, pvalue = pvalue,
                       estimate = list(Xbar = Xbar, Ybar = Ybar, S1 = S1, S2 = S2),
                       construct = list(Z = Z, Zbar = Zbar, Sz = Sz)))
      }else{
        # Two-Sample Means Tests with Unequal and Unknown Covariance Matrix and Unequal Sample Sizes
        if(n<m){
          mid_amout <- 1/sqrt(n*m)*apply(Y[1:n,], 2, sum)
          Z <- X - sqrt(n/m)*Y[1:n,] + matrix(mid_amout - Ybar, n, p, byrow = TRUE)
          Sz <- cov(Z)
          Zbar <- apply(Z, 2, mean)
          t2 <- n*Zbar%*%solve(Sz,Zbar)
          statistic <- (n-p)/((n-1)*p)*t2
          pvalue <- 1-pf(statistic, df1 = p, df2 = n-p)
        }else{
          mid_amout <- 1/sqrt(n*m)*apply(X[1:m,], 2, sum)
          Z <- Y - sqrt(m/n)*X[1:m,] + matrix(mid_amout - Xbar, m, p, byrow = TRUE)
          Sz <- cov(Z)
          Zbar <- apply(Z, 2, mean)
          t2 <- m*Zbar%*%solve(Sz,Zbar)
          statistic <- (m-p)/((m-1)*p)*t2
          pvalue <- 1-pf(statistic, df1 = p, df2 = m-p)
        }

        cat("Two-sample mean test when covariance matrix is unequal
            and unknown with unequal size.\n")
        cat("F-statistic = ", statistic, "df1 = ", p, "df2 = ", min(m,n)-p,
            "p-value = ", pvalue, "\n")
        cat("alternative hypothesis: true mean of two population is not equal","\n")
        invisible(list(statistic = statistic, pvalue = pvalue,
                       estimate = list(Xbar = Xbar, Ybar = Ybar, S1 = S1, S2 = S2),
                       construct = list(Z = Z, Zbar = Zbar, Sz = Sz)))
      }
    }
  }
}


#' Multiple Population Mean Test
#'
#' This function carries out hypothesis test for the means of two or more multivariate normal distributions.
#' @param mydata a matrix/data frame/tibble which contains all the data with a variable to
#' distinguish populations.
#' @param label a character denotes the name of the distinguish variable.
#' @export
mvnorm_multi_mu_test <- function(mydata, label){
  group_data <- mydata %>%
    group_by(.data[[label]]) %>%
    group_split() %>%
    map2(c(label),drop_label)
  mydata <- drop_label(mydata, label)
  ni <- group_data %>% map_dbl(nrow)
  Ai <- group_data %>% map(cov) %>% map2(ni-1,function(x,y) x*y)
  A <- Reduce('+', Ai)
  n1 <- sum(ni-1)
  n2 <- length(ni)-1
  p <- ncol(group_data[[1]])
  n <- sum(ni)
  # T = A + B, but it can be computed directly
  t <- mydata %>% cov() %>% `*`(n-1)
  # No need to count, but I want to record the algorithm
  # Xbar <- mydata %>% drop_label("ind") %>% map_dbl(mean)
  # Bi <- group_data %>% map(function(x) map_dbl(x,mean)) %>%
  #   map2(ni,function(x,y) (y+1)*(x-Xbar)%*%t(x-Xbar))
  # B <- Reduce('+', Bi)
  wilks <- det(A)/det(t)
  pvalue <- pWilks(wilks,p = p,n1 = n1, n2 = n2, lower.tail = FALSE)
  cat("Multi-sample mean test when covariance matrix is equal but unknown.\n")
  cat("Wilks-statistic = ", wilks, "p = ", p, "n1 = ", n1, "n2 = ", n2,
      "p-value = ", pvalue, "\n")
  cat("alternative hypothesis: there exists a pair of unequal mean","\n")
  invisible(list(wilks = wilks, pvalue = pvalue))
}
