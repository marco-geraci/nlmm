\name{vcov.nlmm}
\alias{vcov.nlmm}

\title{
Calculate Variance-Covariance Matrix for a Fitted Generalized Laplace Mixed-Effects Object
}

\description{
Returns the variance-covariance matrix of the all the parameters of a fitted \code{nlmm} object.
}
\usage{
\method{vcov}{nlmm}(object, ...)
}

\arguments{
  \item{object}{
an \code{nlmm} object.
}
  \item{\dots}{
not used.
}
}

\details{
Gives the variance-covariance matrix of the GLME estimator, on the scale of the unconstrained, unrestricted parameters. The size is \eqn{d x d}, \eqn{d = p + r + 2 + 1 + s}, with \eqn{p} fixed coefficients, \eqn{r} non-redundant parameters of the random effects distribution, 2 shape parameters, 1 scale parameter, \eqn{s} parameters of the residual variance function (if specified in the model), in this order.


}

\value{
a matrix.
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

# Number of parameters is d = 6 + 3 + 2 + 1 + 0 = 12
fit <- nlmm(y ~ trt*time, random = ~ time, group = id, data = rats,
cov = "pdSymm", control = nlmmControl(multistart = FALSE))
fit$par
vcov(fit)

# Number of parameters is d = 6 + 1 + 2 + 1 + 4 = 14
fit <- nlmm(y ~ trt*time, random = ~ 1, group = id, data = rats,
control = nlmmControl(multistart = FALSE), weights = varIdent(form = ~ 1|time))
fit$par
vcov(fit)

# Number of parameters is d = 6 + 1 + 0 + 1 + 0 = 8
# Note that the shape parameters are now constrained
fit <- nlmm(y ~ trt*time, random = ~ 1, group = id, data = rats,
control = nlmmControl(alpha.index = 0, multistart = FALSE))
fit$par
vcov(fit)

}
}

\keyword{covariance}
\keyword{standard errors}