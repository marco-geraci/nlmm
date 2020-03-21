\name{summary.nlmm}
\alias{summary.nlmm}

\title{Summary for an \code{nlmm} Object}

\description{
Summary method for class \code{nlmm}.}

\usage{
\method{summary}{nlmm}(object, alpha = 0.05, ...)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{nlmm}.
	}

  \item{alpha}{
significance level.
	}
  \item{\dots}{
not used.
	}
}

\details{
\code{print.summary.nlmm} formats the coefficients, standard errors, etc. and additionally gives `significance stars'.}

\value{
an object of class \code{summary.nlmm}. The function \code{summary.nlmm} computes and returns a list of summary statistics of the fitted generalized Laplace mixed-effects model given in \code{object}, using the components (list elements) from its argument, plus

\item{tTable}{a matrix with estimates, standard errors, etc.}
}

\author{
Marco Geraci}

\seealso{
\code{\link{print.summary.nlmm}}
\code{\link{nlmm}}
}


\examples{
\dontrun{
data(rats)
fit <- nlmm(y ~ trt*time, random = ~ 1, group = id, data = rats,
control = nlmmControl(multistart = FALSE))
summary(fit)

}
}

\keyword{standard errors}
