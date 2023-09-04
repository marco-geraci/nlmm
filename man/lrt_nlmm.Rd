\name{lrt_nlmm}
\alias{lrt_nlmm}
\alias{print.lrt_nlmm}

\title{
Likelihood Ratio Test for Generalized Laplace Mixed-Effects Models
}
\description{
This function is used to perform a likelihood ratio test for two fitted generalized Laplace mixed-effects models.
}
\usage{
lrt_nlmm(object0, object1)
\method{print}{lrt_nlmm}(x, digits = max(3, getOption("digits") - 3), ...)
}

\arguments{
  \item{object0}{
object of class \code{nlmm} with estimates of the constrained model.
}
  \item{object1}{
object of class \code{nlmm} with estimates of the unconstrained model.
}
  \item{x}{
a \code{lrt_nlmm} object.
}
  \item{digits}{
a non-null value for digits specifies the minimum number of significant digits to be printed in values.
}

  \item{\dots}{
not used.
}
}

\note{

The function \code{lrt_nlmm} is a wrapper for routines developed by Alessio Farcomeni.

}

\value{
An object of class \code{lrt_nlmm} is a list containing the following components:

\item{statistic}{the value of the test statistic}
\item{p.value}{the p-value of the test}
\item{df}{either the degrees of freedom of a chi-squared test or the weights of a chi-bar-squared test}
\item{V}{the matrix based on which the weights for the chi-bar-squared test are calculated}
\item{alpha}{values of the shape parameter in the constrained object}
\item{alpha.index}{index of the constrained shape parameter}
\item{chibar}{logical flag. If \code{TRUE}, the test statistic is a chi-bar}

}
\references{
Geraci, M. and Farcomeni A. (2020). A family of linear mixed-effects models using the generalized Laplace distribution. Statistical Methods in Medical Research, 29(9), 2665-2682.
}
\author{
Marco Geraci and Alessio Farcomeni
}

\seealso{
\code{\link{nlmm}}
}


\keyword{Likelihood ratio test}
