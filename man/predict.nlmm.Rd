\name{predict.nlmm}
\alias{predict.nlmm}

\title{
Predictions from an \code{nlmm} Object
}

\description{
The predictions at level 0 correspond to predictions based only on the fixed effects estimates. The predictions at level 1 are obtained by adding the best linear predictions of the random effects to the predictions at level 0.}

\usage{
\method{predict}{nlmm}(object, level = 0, ...)
}

\arguments{
  \item{object}{
an \code{nlmm} object.
}
  \item{level}{
an integer vector giving the level of grouping to be used in obtaining the predictions.
}
  \item{\dots}{
not used.
}
}


\value{
a vector of predictions.
}

\references{
Geraci, M. and Farcomeni A (2020). A family of linear mixed-effects models using the generalized Laplace distribution. Statistical Methods in Medical Research, \url{https://doi.org/10.1177/0962280220903763}.
}

\author{
Marco Geraci
}


\seealso{
\code{\link{nlmm}}, \code{\link{ranef.nlmm}}, \code{\link{fixef.nlmm}}
}

\examples{
\dontrun{
data(rats)
fit <- nlmm(y ~ trt*time, random = ~ 1, group = id, data = rats,
control = nlmmControl(multistart = FALSE))

# Individual growth trajectories
predict(fit, level = 1)
}
}


\keyword{prediction}
