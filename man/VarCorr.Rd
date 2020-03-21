\name{VarCorr.nlmm}
\alias{VarCorr}
\alias{VarCorr.nlmm}

\title{
Extract Variance-Covariance Matrix
}
\description{
This function extracts the variance-covariance matrix of the random effects from a fitted \code{nlmm} object.
}
\usage{
\method{VarCorr}{nlmm}(x, sigma = NULL, ...)
}

\arguments{
  \item{x}{
an object of \code{\link{class}} "nlmm".
}

 \item{sigma}{
not used.
}

 \item{...}{
not used.
}

}

\details{
This function returns the variance or the variance-covariance matrix of the random effects. The generic function \code{VarCorr} is imported from the \code{nlme} package (Pinheiro et al, 2014).
}

\references{
Pinheiro J, Bates D, DebRoy S, Sarkar D and R Core Team (2014). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-117, \url{https://CRAN.R-project.org/package=nlme}.
}

\author{
Marco Geraci
}

\seealso{
\code{\link{nlmm}}
}


\examples{
\dontrun{
data(rats)
fit <- nlmm(y ~ trt*time, random = ~ time, group = id, data = rats, cov = "pdSymm",
control = nlmmControl(multistart = FALSE))

# Symmetric variance-covariance of random intercepts and slopes
VarCorr(fit)
}
}

\keyword{covariance}
\keyword{coefficients}
