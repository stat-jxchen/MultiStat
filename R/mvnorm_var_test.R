#' Single Population Covariance Matrix Test
#'
#' This function carries out hypothesis test for the covariance matrix of a multivariate normal distributions.
#'
#' @param X a matrix or a data frame which denotes a sample.
#' @param Sig0_type a character denotes the type of covariance matrix of null
#' hypothesis.
#' @param Sig0 a matrix denotes the covariance matrix under null hypothesis.
#' @examples
#' # for general case: Sigma = Sigma0
#' mu <- c(10,20,30)
#' sig <- diag(c(2,5,9))
#' sim_time <- 1000
#' pvalue <- numeric(sim_time)
#' for (i in seq_len(sim_time)) {
#' X <- MASS::mvrnorm(n = 200, mu = mu, Sigma = sig)
#' re <- mvnorm_var_test(X, Sig0_type = "general", Sig0 = sig)
#' pvalue[i] <- as.vector(re$pvalue)
#' }
#' sum(pvalue<0.05)/sim_time
#' # for sphere case: Sigma = sigma^2*Sigma0
#' sigma2 <- 1.5
#' sim_time <- 1000
#' pvalue <- numeric(sim_time)
#' for (i in seq_len(sim_time)) {
#' X <- MASS::mvrnorm(n = 200, mu = mu, Sigma = sig*sigma2)
#' re <- mvnorm_var_test(X, Sig0_type = "sphere", Sig0 = sig)
#' pvalue[i] <- as.vector(re$pvalue)
#' }
#' sum(pvalue<0.05)/sim_time
#' @export
mvnorm_var_test <- function(X, Sig0_type, Sig0){
  n <- nrow(X)
  p <- ncol(X)
  A <- cov(X)*(n-1)
  if(Sig0_type == "general"){
    miditem <- A%*%mysolve(Sig0)
    llambda <- (-tr(miditem)/2)+(n/2)*log(det(miditem))+(n*p/2)*(1-log(n))
    statistic <- -2*llambda
    pvalue <- 1-pchisq(statistic,df = p*(p+1)/2)
    cat("One-sample covariance matrix test when covariance matrix is known.\n")
    cat("Likelihood ratio statistic = ", exp(llambda), "\n")
    cat("chisq-statistic = ", statistic, "df = ", p*(p+1)/2, "p-value = ", pvalue, "\n")
    cat("alternative hypothesis: true covariance matrix is not equal to Sig0","\n")
    invisible(list(statistic = statistic, pvalue = pvalue))
  }else if(Sig0_type == "sphere"){
    miditem <- mysolve(Sig0)%*%A
    sigma2_est <- 1/(n*p)*tr(miditem)
    llambda <- (n/2)*log(det(miditem))-(n*p/2)*(log(tr(miditem))-log(p))
    statistic <- -2/n*(n-1-(2*p^2+p+2)/(6*p))*llambda
    pvalue <- 1-pchisq(statistic,df = p*(p+1)/2-1)
    cat("One-sample covariance matrix test when covariance matrix is known but with an
        unkonwn coefficient.\n")
    cat("Likelihood ratio statistic = ", exp(llambda), "\n")
    cat("chisq-statistic = ", statistic, "df = ", p*(p+1)/2-1, "p-value = ", pvalue, "\n")
    cat("alternative hypothesis: true covariance matrix is not equal to sigma2*Sig0","\n")
    invisible(list(statistic = statistic, pvalue = pvalue, sigma2_est = sigma2_est))
  }
}

#' Multiple Population Covariance Matrix Test
#'
#' This function carries out hypothesis test for the covariance matrices of two or more multivariate normal distributions.
#'
#' @param mydata a matrix/data frame/tibble which contains all the data with a variable to
#' distinguish populations.
#' @param label a character denotes the name of the distinguish variable.
#' @examples
#' mvnorm_multi_var_test(health,"ind")
#' @export
mvnorm_multi_var_test <- function(mydata, label){
  group_data <- mydata %>%
    group_by(.data[[label]]) %>%
    group_split() %>%
    map2(c(label),drop_label)
  mydata <- drop_label(mydata, label)
  p <- ncol(mydata)
  n <- nrow(mydata)
  ni <- group_data %>% map_dbl(nrow)
  k <- length(ni)
  Ai <- group_data %>% map(cov) %>% map2(ni-1,function(x,y) x*y)
  A <- Reduce('+', Ai)
  miditem <- map2(ni-1,Ai,function(x,y) x*log(det(y/x))) %>% reduce(sum)
  M <- (n-k)*log(det(A/(n-k)))-miditem
  lambda <- exp(-M/2)
  d <- (2*p^2+3*p-1)/(6*(p+1)*(k-1))*(sum(1/(ni-1))-1/(n-k))
  statistic <- (1-d)*M
  f <- p*(p+1)*(k-1)/2
  pvalue <- 1-pchisq(statistic, df = f)
  cat("Multi-sample covariance matrix test.\n")
  cat("Likelihood ratio statistic = ", lambda, "\n")
  cat("chisq-statistic = ", statistic, "df = ", f, "p-value = ", pvalue, "\n")
  cat("alternative hypothesis: there exists a pair of unequal covariance matrices","\n")
  invisible(list(statistic = statistic, pvalue = pvalue))
}

#' Multiple Simultaneous Test
#'
#' This function carries out simultaneous hypothesis test for the means and covariance matrices of one or two multivariate normal distributions.
#'
#' @param mydata a matrix/data frame/tibble which contains all the data with a variable to
#' distinguish populations.
#' @param label a character denotes the name of the distinguish variable.
#' @examples
#' mvnorm_multi_simult_test(health,"ind")
#' @export
mvnorm_multi_simult_test <- function(mydata, label){
  group_data <- mydata %>%
    group_by(.data[[label]]) %>%
    group_split() %>%
    map2(c(label),drop_label)
  mydata <- drop_label(mydata, label)
  p <- ncol(mydata)
  n <- nrow(mydata)
  t <- cov(mydata)*(n-1)
  ni <- group_data %>% map_dbl(nrow)
  k <- length(ni)
  Ai <- group_data %>% map(cov) %>% map2(ni-1,function(x,y) x*y)
  miditem1 <- map2(ni-1,Ai,function(x,y) x/2*log(det(y))) %>% reduce(sum)
  miditem2 <- sum((ni-1)*p/2*log(ni-1))
  llambda <- miditem1-(n-k)/2*log(det(t))+(n-k)*p/2*log(n-k)-miditem2
  f <- p*(p+3)*(k-1)/2
  b <- (sum(1/(ni-1))-1/(n-k))*(2*p^2+3*p-1)/(6*(p+3)*(k-1))-(p-k+2)/((n-k)*(p+3))
  statistic <- -2*(1-b)*llambda
  pvalue <- 1-pchisq(statistic, df = f)
  cat("Multi-sample simultaneous test.\n")
  cat("Likelihood ratio statistic = ", exp(llambda), "\n")
  cat("chisq-statistic = ", statistic, "df = ", f, "p-value = ", pvalue, "\n")
  cat("alternative hypothesis: there exists a pair of unequal mean or
      unequal covariance matrices","\n")
  invisible(list(statistic = statistic, pvalue = pvalue))
}
