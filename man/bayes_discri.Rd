% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayes_discri.R
\name{bayes_discri}
\alias{bayes_discri}
\title{Bayes Discrimination}
\usage{
bayes_discri(x, prior, loss, pdf, para_list)
}
\arguments{
\item{x}{a vector, a single sample data to be predicted.}

\item{prior}{a vector, prior probability of each class.}

\item{loss}{a matrix, the ijth element denotes the loss of misjudging
class j as class i. If the loss is equal, it could be set to a number.}

\item{pdf}{a function, the common form of pdf of each class.}

\item{para_list}{a list, each element is an additional single parameter of the pdf.}
}
\description{
This function do Bayes discrimination.
}
\examples{
{
# ToDo
}
}
