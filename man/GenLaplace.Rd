\name{GenLaplace}
\alias{GenLaplace}
\alias{dgl}
\alias{pgl}
\alias{qgl}
\alias{rgl}

\title{
The Univariate Symmetric Generalized Laplace Distribution
}

\description{
Density, distribution function, quantile function and random generation for the univariate symmetric generalized Laplace distribution.}

\usage{
dgl(x, sigma = 1, shape = 1, log = FALSE)
pgl(x, sigma = 1, shape = 1, lower.tail = TRUE, log.p = FALSE)
qgl(p, sigma = 1, shape = 1, lower.tail = TRUE, log.p = FALSE)
rgl(n, sigma = 1, shape = 1)
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
The univariate symmetric generalized Laplace distribution (Kotz et al, 2001, p.190) has density
\deqn{
f(x) =
\frac{2}{\sqrt{2\pi}\Gamma(s)\sigma^{s+1/2}}(\frac{|x|}{\sqrt{2}})^{\omega}B_{\omega}(\frac{\sqrt{2}|x|}{\sigma})
}
where \eqn{\sigma} is the scale parameter, \eqn{\omega = s - 1/2}, and \eqn{s} is the shape parameter. \eqn{\Gamma} denotes the Gamma function and \eqn{B_{u}} the modified Bessel function of the third kind with index \eqn{u}. The variance is \eqn{s\sigma^{2}}.

This distribution is the univariate and symmetric case of \code{\link{MultivariateGenLaplace}}.
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
