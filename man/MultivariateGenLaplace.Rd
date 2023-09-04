\name{MultivariateGenLaplace}
\alias{MultivariateGenLaplace}
\alias{dmgl}
\alias{rmgl}

\title{
The Multivariate Asymmetric Generalized Laplace Distribution
}

\description{
Density and random generation for the multivariate asymmetric generalized Laplace distribution.}

\usage{
dmgl(x, mu = rep(0, d), sigma = diag(d), shape = 1, log = FALSE)
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
simmetry parameter.
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
This is the distribution described by Kozubowski et al (2013) and has density
\deqn{
f(x) =
\frac{2\exp(\mu'\Sigma^{-1}x)}{(2\pi)^{d/2}\Gamma(s)|\Sigma|^{1/2}}(\frac{Q(x)}{C(\Sigma,\mu)})^{\omega}B_{\omega}(Q(x)C(\Sigma,\mu))
}
where \eqn{\mu} is the symmetry parameter, \eqn{\Sigma} is the scale parameter, \eqn{Q(x)=\sqrt{x'\Sigma^{-1}x}}, \eqn{C(\Sigma,\mu)=\sqrt{2+\mu'\Sigma^{-1}\mu}}, \eqn{\omega = s - d/2}, \eqn{d} is the dimension of \eqn{x}, and \eqn{s} is the shape parameter (note that the parameterization in \code{\link{nlmm}} is \eqn{\alpha = \frac{1}{s}}). \eqn{\Gamma} denotes the Gamma function and \eqn{B_{u}} the modified Bessel function of the third kind with index \eqn{u}. The parameter \eqn{\mu} is related to the skewness of the distribution (symmetric if \eqn{\mu = 0}). The variance-covariance matrix is \eqn{s(\Sigma + \mu\mu')}. The multivariate asymmetric Laplace is obtained when \eqn{s = 1} (see \code{\link{MultivariateLaplace}}). 

In the symmetric case (\eqn{\mu = 0}), the multivariate GL distribution has two special cases: multivariate normal for \eqn{s \rightarrow \infty} and multivariate symmetric Laplace for \eqn{s = 1}.

The \bold{univariate symmetric} GL distribution is provided via \code{\link{GenLaplace}}, which gives the distribution and quantile functions in addition to the density and random generation functions.
}


\value{

\code{dmgl} gives the GL density of a \eqn{d}-dimensional vector \code{x}. \code{rmgl} generates a sample of size \code{n} of \eqn{d}-dimensional random GL variables.
}

\references{
Geraci, M. and Farcomeni A. (2020). A family of linear mixed-effects models using the generalized Laplace distribution. Statistical Methods in Medical Research, 29(9), 2665-2682.

Kozubowski, T. J., K. Podgorski, and I. Rychlik (2013). Multivariate generalized Laplace distribution and related random fields. Journal of Multivariate Analysis 113, 59-72.
}
\author{
Marco Geraci
}


\seealso{
\code{\link{GenLaplace}}
}

\keyword{Generalized Laplace Distribution}
