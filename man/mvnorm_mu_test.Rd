% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvnorm_mu_test.R
\name{mvnorm_mu_test}
\alias{mvnorm_mu_test}
\title{Single or Bivariate Population Mean Test}
\usage{
mvnorm_mu_test(X, Y = NULL, mu0 = NULL, Sig0 = NULL)
}
\arguments{
\item{X}{a matrix or a data frame which denotes a sample.}

\item{Y}{a matrix or a data frame which denotes the second sample.}

\item{mu0}{a vector, the value of mu under null hypothesis.}

\item{Sig0}{a matrix or string "unequal (only when \code{Y} is not \code{NULL}),
when it is a matrix, it will be used in test and means covariance matrices are
equal and known, when it is \code{unequal}, it means the covariance matrix of two
population are not equal, and when it is \code{NULL} means covariance matrices
are equal but unknown.}
}
\description{
This function carries out hypothesis test for the means of one or two multivariate normal distributions.
}
\note{
When \code{Y} is \code{NULL}, it will carry out one-sample test. If \code{Y}
is not \code{NULL}, \code{mu0} can be \code{NULL}.
}
\examples{
# One-Sample Mean Test with Unknown Covariance Matrix
mu0 <- c(4, 50, 10)
mvnorm_mu_test(X = sweat, mu0 = mu0)
# One-Sample Mean Test with known Covariance Matrix
sig0 <- cov(sweat)
mvnorm_mu_test(sweat,mu0 = mu0, Sig0 = sig0)
# Two-Sample Means Test When Covariance Matrices Are Equal But Unknown
mvnorm_mu_test(X = EconomyJapan, Y = EconomyUSA)
# Two-Sample Means Test When Covariance Matrices Are Equal and Known
X <- rbind(EconomyJapan, EconomyUSA)
Sig0 <- var(X)
mvnorm_mu_test(X = EconomyJapan, Y = EconomyUSA, Sig0 = Sig0)
# Two-Sample Means Test with Unequal and Unknown Covariance Matrix and Unequal Sample Sizes
X <- subset(milk,V4 == "gasoline")[,-4]
Y <- subset(milk,V4 == "diesel")[,-4]
mvnorm_mu_test(X = X, Y = Y, Sig0 = "unequal")
}
