\name{MultivariateGenLaplace}
\alias{MultivariateGenLaplace}
\alias{dmgl}
\alias{rmgl}

\title{
The Multivariate Symmetric Generalized Laplace Distribution
}

\description{
Density, distribution function, quantile function and random generation for the multivariate generalized Laplace distribution.}

\usage{
dmgl(x, mu = rep(0, n), sigma = diag(n), shape = 1, log = FALSE)
rmgl(n, mu, sigma, shape = 1)
}

\arguments{
  \item{x}{
vector of quantiles.
}
  \item{n}{
number of observations.
}
  \item{mu}{
location parameter.
}
  \item{sigma}{
scale parameter -- positive-definite matrix.
}
  \item{shape}{
shape parameter.
}
  \item{log}{
logical; if \code{TRUE}, probabilities are log--transformed.
}

}
\details{

This is the multivariate extension of the (univariate) generalized Laplace distribution (\code{\link{GenLaplace}}).
}


\references{
Kozubowski, T. J., K. Podgorski, and I. Rychlik (2013). Multivariate generalized Laplace distribution and related random fields. Journal of Multivariate Analysis 113, 59-72.
}
\author{
Marco Geraci
}


\seealso{
\code{\link{GenLaplace}}
}

\keyword{Generalized Laplace Distribution}
