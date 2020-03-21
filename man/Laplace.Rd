\name{Laplace}
\alias{Laplace}
\alias{dl}
\alias{rl}

\title{
The Laplace Distribution
}

\description{
Density and random generation for the (symmetric) Laplace distribution.}

\usage{
dl(x, mu = 0, sigma = 1, log = FALSE)
rl(n, mu = 0, sigma = 1)
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
positive scale parameter.
}
  \item{log}{
logical; if \code{TRUE}, probabilities are log--transformed.
}

}
\details{

The Laplace distribution has density
\deqn{
f(x) =
\frac{1}{\sqrt{2}\sigma}e^{-\frac{sqrt(2)}{\sigma} |x - \mu|}
}
where \eqn{\mu} is the location parameter and \eqn{\sigma} is the scale parameter.
}

\value{

\code{dl} gives the density and \code{rl} generates random deviates.

}

\references{
Kotz, S., Kozubowski, T., and Podgorski, K. (2001). The Laplace distribution and generalizations. Boston, MA: Birkhauser.
}
\author{
Marco Geraci
}


\seealso{
\code{\link{MultivariateLaplace}}, \code{\link{GenLaplace}}
}

\keyword{Laplace Distribution}
