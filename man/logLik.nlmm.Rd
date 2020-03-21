\name{logLik.nlmm}
\alias{logLik.nlmm}

\title{
Extract Log-Likelihood
}
\description{
\code{logLik.nlmm} extracts the log-likelihood of a fitted nlmm.
}
\usage{
\method{logLik}{nlmm}(object, ...)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} "nlmm".

}
  \item{\dots}{
not used.
}
}

\value{

Returns the loglikelihood of the fitted model. This is a number with at one attribute, "df" (degrees of freedom), giving the number of (estimated) parameters in the model.

}


\author{
Marco Geraci
}


\seealso{
\code{\link{nlmm}}
}

