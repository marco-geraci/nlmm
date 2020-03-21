\name{ranef.nlmm}
\alias{ranef}
\alias{ranef.nlmm}

\title{
Extract Random Effects
}
\description{
This function computes random effects for a linear quantile mixed model.
}
\usage{
\method{ranef}{nlmm}(object, ...)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{nlmm}.
}
  \item{\dots}{
not used.
}

}

\details{
The prediction of the random effects is done via estimated best linear prediction (Geraci and Farcomeni, 2019). The generic function \code{ranef} is imported from the \code{nlme} package (Pinheiro et al, 2014).
}

\value{
a data frame of predicted random effects.
}

\references{
Geraci, M. and Farcomeni A (2020). A family of linear mixed-effects models using the generalized Laplace distribution. Statistical Methods in Medical Research, \url{https://doi.org/10.1177/0962280220903763}.

Pinheiro J, Bates D, DebRoy S, Sarkar D and R Core Team (2014). nlme: Linear and Nonlinear Mixed Effects Models. R package version 3.1-117, \url{https://CRAN.R-project.org/package=nlme}.
}

\author{
Marco Geraci
}

\seealso{
\code{\link{nlmm}}, \code{\link{fixef.nlmm}}
}


\examples{
\dontrun{
data(rats)
fit <- nlmm(y ~ trt*time, random = ~ time, group = id, data = rats,
control = nlmmControl(multistart = FALSE))

# Predicted random intercepts and slopes
ranef(fit)
}
}

\keyword{random effects}
\keyword{coefficients}
