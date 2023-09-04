\name{generate.data}
\alias{generate.data}

\title{
Simulate Data from Mixed-Effects Models
}
\description{
This function generates data from a 2-level hierarchical design.
}
\usage{
generate.data(R, n, M, sigma_1 = NULL, sigma_2 = NULL,
	shape_1 = NULL, shape_2 = NULL, dist.u, dist.e,
	beta, gamma, fixed = FALSE, seed = round(runif(1,1,1000)))
}

\arguments{
  \item{R}{
number of replications.
}
  \item{n}{
number of observations within cluster.
}

  \item{M}{
number of clusters.
}

  \item{sigma_1}{
scale parameter for the random effects.
}
  \item{sigma_2}{
scale parameter for the errors.
}
  \item{shape_1}{
shape parameter for the random effects.
}
  \item{shape_2}{
shape parameter for the errors.
}
  \item{dist.u}{
distribution of the random effects.
}

  \item{dist.e}{
distribution of the errors.
}

  \item{beta}{
vector of coefficients for fixed effects.
}
  \item{gamma}{
vector of coefficients for heteroscedasticity.
}
  \item{fixed}{
logical flag. See details.
}

  \item{seed}{
seed for random number generation.
}

}

\details{

This function generates data as in the simulation study by Geraci and Farcomeni (2020). The data-generating model is
\deqn{
y[ij] = \beta[0] + \beta[1]x[ij] + \beta[2]z[ij] + u[i] + v[i]x[ij] + (\gamma[0] + \gamma[1]x[ij])e[ij]
}
where \eqn{(u[i],v[i])} follows a distribution with scale \code{sigma_1} and shape \code{shape_1}, and \eqn{e} follows a distribution with scale \code{sigma_2} and shape \code{shape_2}.

The scale parameter \code{sigma_1} must be a 1 by 1 or a 2 by 2 matrix. In the former case, the model will include only random intercepts. In the latter case, then both random intercepts and slopes will be included. Currently, no more than 2 random effects can be specified. The scale parameter \code{sigma_2} must be a matrix \code{n} by \code{n}.

The options for \code{dist.u} and \code{dist.e} are: multivariate normal ("norm") (\code{\link{rmvnorm}}), multivariate symmetric Laplace ("laplace") (\code{\link{rmal}}), multivariate symmetric generalized Laplace ("genlaplace") \code{\link{rmgl}}, and multivariate Student's t ("t") (\code{\link{rmvt}}).

The shape parameter specifies the degrees of freedom for Student's t and chi-squared, and the kurtosis of the generalized Laplace.

The values \eqn{x[ij]} are generated as \eqn{x[ij] = \delta[i] + \zeta[ij]}, where \eqn{\delta[i]} and \eqn{\zeta[ij]} are independent standard normal. If the argument \code{fixed = TRUE}, then \eqn{x[ij] = j}. The values \eqn{z[ij]} are generated from Bernoullis with probability 0.5.

}

\value{
\code{nlmm} returns an object of \code{\link{class}} \code{nlmm}.

The function \code{summary} is used to obtain and print a summary of the results.

An object of class \code{nlmm} is a list containing the following components:

\item{Y}{a matrix \eqn{R x N}, where \eqn{N = n x M}, with responses}
\item{X}{an array \eqn{N x 3 x R} with fixed design matrix}
\item{group}{vector of length \eqn{N} with cluster labels}
\item{u}{an array \eqn{M x 2 x R} with random effects}
\item{e}{a matrix \eqn{R x N} with errors}
}


\references{
Geraci, M. and Farcomeni A. (2020). A family of linear mixed-effects models using the generalized Laplace distribution. Statistical Methods in Medical Research, 29(9), 2665-2682.
}

\author{
Marco Geraci
}


\seealso{
\code{\link{nlmm}}
}


\examples{

# Simulate 10 replications from a homoscedastic normal mixed model.
generate.data(R = 10, n = 3, M = 5, sigma_1 = diag(2), sigma_2 = diag(3),
shape_1 = NULL, shape_2 = NULL, dist.u = "norm", dist.e = "norm",
beta = c(1,2,1), gamma = c(1,0))

# Simulate 10 replications from a generalized Laplace. Note: the shape
# parameter that is passed to rmgl corresponds to the reciprocal of the
# parameter alpha in Geraci and Farcomeni (2020)
generate.data(R = 10, n = 3, M = 5, sigma_1 = diag(2), sigma_2 = diag(3),
shape_1 = 1/0.5, shape_2 = 1/0.5, dist.u = "genlaplace", dist.e = "genlaplace",
beta = c(1,2,1), gamma = c(1,0))
 
}

\keyword{Simulation}
