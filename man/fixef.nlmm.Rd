\name{fixef.nlmm}
\alias{fixef}
\alias{fixef.nlmm}

\title{
Extract Generalized Mixed-Effects Models Coefficients
}
\description{
\code{fixef} extracts estimated fixed effects from \code{nlmm} objects.
}
\usage{
\method{fixef}{nlmm}(object, ...)
}

\arguments{
  \item{object}{
a \code{nlmm} object.
}
  \item{\dots}{
not used.
}
}

\value{
a vector of estimated fixed effects. 
}

\author{
Marco Geraci
}

\seealso{
\code{\link{nlmm}}
\code{\link{summary.nlmm}}
}

\keyword{coefficients}

