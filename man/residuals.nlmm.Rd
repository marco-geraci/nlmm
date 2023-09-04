\name{residuals.nlmm}
\alias{residuals.nlmm}

\title{
Residuals from an \code{nlmm} Object
}

\description{
The residuals at level 0 correspond to population residuals (based only on the fixed effects estimates). The residuals at level 1 are obtained by adding the best linear predictions of the random effects to the predictions at level 0 and the subtracting these from the model response.
}

\usage{
\method{residuals}{nlmm}(object, level = 0, ...)
}


\arguments{
  \item{object}{
an \code{nlmm} object.
}
  \item{level}{
an optional integer vector giving the level of grouping to be used in obtaining the predictions. Level zero corresponds to the population residuals.
}
  \item{\dots}{
not used.
}
}

\value{
a matrix of residuals.
}

\references{
Geraci, M. and Farcomeni A. (2020). A family of linear mixed-effects models using the generalized Laplace distribution. Statistical Methods in Medical Research, 29(9), 2665-2682.
}

\author{
Marco Geraci}


\seealso{
\code{\link{nlmm}}, \code{\link{predict.nlmm}}, \code{\link{fixef.nlmm}}, \code{\link{ranef.nlmm}},
}

\keyword{residuals}
