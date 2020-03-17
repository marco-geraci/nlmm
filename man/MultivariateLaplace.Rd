\name{MultivariateLaplace}
\alias{MultivariateLaplace}
\alias{dmal}
\alias{rmal}

\title{
The Multivariate Asymmetric Laplace Distribution
}

\description{
Density and random generation for the multivariate asymmetric Laplace distribution.}

\usage{
dmal(x, m = rep(0, nrow(sigma)), sigma, log = FALSE)
rmal(n, m = rep(0, nrow(sigma)), sigma)
}

\arguments{
  \item{x}{
vector of quantiles.
}
  \item{n}{
number of observations.
}
  \item{m}{
asymmetry parameter.
}
  \item{sigma}{
scale parameter -- positive-definite matrix.
}
  \item{log}{
logical; if \code{TRUE}, probabilities are log--transformed.
}
}
\details{

This is the multivariate extension of the (univariate) asymmetric Laplace distribution.
}


\references{
Kotz, S., Kozubowski, T., and Podgorski, K. (2001). The Laplace distribution and generalizations. Boston, MA: Birkhauser.
}
\author{
Marco Geraci
}


\seealso{
\code{\link{Laplace}}, \code{\link{MultivariateGenLaplace}}
}

\keyword{Multivariate Asymmetric Laplace Distribution}
