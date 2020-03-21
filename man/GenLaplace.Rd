\name{GenLaplace}
\alias{GenLaplace}
\alias{dgl}
\alias{pgl}
\alias{qgl}
\alias{rgl}

\title{
The Symmetric Generalized Laplace Distribution
}

\description{
Density, distribution function, quantile function and random generation for the generalized Laplace distribution.}

\usage{
dgl(x, mu = 0, sigma = 1, shape = 1, log = FALSE)
pgl(x, mu = 0, sigma = 1, shape = 1, lower.tail = TRUE, log.p = FALSE)
qgl(p, mu = 0, sigma = 1, shape = 1, lower.tail = TRUE, log.p = FALSE)
rgl(n, mu = 0, sigma = 1, shape = 1)
}

\arguments{
  \item{x}{
vector of quantiles.
}
  \item{p}{
vector of probabilities.
}
  \item{n}{
number of observations.
}
  \item{mu}{
location parameter.
}
  \item{sigma}{
positive scale parameter.
}
  \item{shape}{
shape parameter.
}
  \item{log,log.p}{
logical; if \code{TRUE}, probabilities are log--transformed.
}
  \item{lower.tail}{
logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.
}

}
\details{

The generalized Laplace distribution has density
\deqn{
f(x) =
\frac{1}{\sqrt{\frac{\pi}{2}}\Gamma(\beta)\sigma^{\beta + 1/2}}(\frac{|y - \mu|}{\sqrt{2}})^{\omega}B[\omega](\frac{\sqrt{2}|y - \mu|}{\sigma})
}
where \eqn{\omega = \beta - 1/2}, \eqn{\mu} is the location parameter, \eqn{\sigma} is the scale parameter, and \eqn{\beta} is the shape parameter (note that the parameterization in \code{\link{nlmm}} is \eqn{\alpha = \frac{1}{\beta}}). \eqn{\Gamma()} denotes the Gamma function and \eqn{B[u]()} the modified Bessel function of the third kind with index \eqn{u}.
}

\value{

\code{dgl} gives the density, \code{pgl} gives the distribution function, \code{qgl} gives the quantile function, and \code{rgl} generates random deviates.

}


\references{
Kotz, S., Kozubowski, T., and Podgorski, K. (2001). The Laplace distribution and generalizations. Boston, MA: Birkhauser.
}
\author{
Marco Geraci
}


\seealso{
\code{\link{MultivariateGenLaplace}}
}

\keyword{Generalized Laplace Distribution}
