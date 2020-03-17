\name{nlmmControl}
\alias{nlmmControl}

\title{
Control parameters for nlmm estimation
}
\description{
A list of parameters for controlling the fitting process.
}
\usage{
nlmmControl(method = "Nelder-Mead", nK = 8, multistart = TRUE,
	grid = c(0.001, 0.5, 0.999), alpha = c(0.5, 0.5), alpha.index = 9,
	lme = TRUE, lmeMethod = "REML", lmeOpt = "nlminb", verbose = FALSE)
}

\arguments{
  \item{method}{
character vector that specifies the optimization algorithm to fit a generalized Laplace mixed-effects model. The default is \code{"Nelder-Mead"}.
}
  \item{nK}{
number of knots for each of the two quadratures.
}
  \item{multistart}{
logical flag. If \code{TRUE} (default), the algorithm is run with multiple starting values for the parameter alpha. See also \code{grid}.
}
  \item{grid}{
a vector of values for multi-start optimization. It can be used in conjunction with constrained estimation.
}
  \item{alpha}{
vector of length 2 with starting values between 0 and 1 for the parameter alpha (ignored if multistart is \code{TRUE}) or values at which alpha is constrained if \code{alpha.index} is one of 0, 1, or 2. The first element is for the shape parameter of the random effects, the second for the error term. See Geraci and Farcomeni (2020).
}
  \item{alpha.index}{
the estimation with parameter alpha is unconstrained if \code{alpha.index} is equal to 9 (default). If equal to 0, both shape parameters (random effects and error term) are constrained during estimation and set equal to \code{alpha}. If equal to 1, the first shape parameter (random effects) is constrained during estimation and set equal to the first element of \code{alpha}. If equal to 2, the second shape parameter (error term) is constrained during estimation and set equal to the second element of \code{alpha}.
}
  \item{lme}{
logical flag. Should \code{\link{lme}} be used to get starting values? If \code{FALSE}, \code{\link{lm}} is used instead.
}
  \item{lmeMethod}{
fitting method for \code{\link{lme}} when obtaining starting values. If "REML" the model is fit by maximizing the restricted log-likelihood. If "ML" the log-likelihood is maximized.
}
  \item{lmeOpt}{
optimization algorithm for \code{\link{lme}}. Either either \code{\link{nlminb}} (the default) or \code{\link{optim}}. This is passed to argument \code{opt} in \code{\link{lmeControl}}.
}
  \item{verbose}{
logical flag. If \code{TRUE}, information about the fitting process is printed out.
}
}

\details{
The estimation algorithm for fitting generalized Laplace mixed-effects (GLME) models is described in Geraci and Farcomeni (2020). For unconstrained estimation, it is recommended to leave the default arguments in \code{nlmmControl} unchanged.

The integrated log-likelihood is maximized with either \code{\link{optim}}, in which case \code{method} has to be one of \code{optim}'s options ("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), or \code{\link{nlminb}}, in which case one must use \code{method = "nlminb"}.

Since the parameter alpha is bidimensional, care should be taken when increasing the number of quadrature knots \code{nK} since the total number of quadrature points is given by \eqn{2^{nK}}{2^nK}. For the same reason, care should be taken when providing the \code{grid} values for multi-start optimization since the total number of starting points will be \eqn{s^{2}}{s^2}, where \code{s = length(grid)}.

If \code{alpha.index} is 1 (or 2), the first (or second) element of the alpha parameter is constrained during estimation and set equal to the corresponding value of \code{alpha}. The element of the alpha parameter that is \emph{unconstrained} is initialized with the corresponding element of \code{alpha} (if \code{multistart} is \code{FALSE}) or with values in \code{grid} (if \code{multistart} is \code{TRUE}).

If \code{alpha.index} is 0, both elements of the alpha parameter are fixed and set equal to \code{alpha}. In this case, the argument \code{multistart} is ignored. If \code{alpha} is \code{c(0,0)}, the corresponding model is Normal-Normal and \code{\link{lme}} is used for fitting (only via maximum likelihood). Note that in this case, \code{lmeOpt} can still be used.
}

\value{
a list of control parameters.
}

\references{
Geraci, M. and Farcomeni A (2020). A family of linear mixed-effects models using the generalized Laplace distribution. Statistical Methods in Medical Research, \url{https://doi.org/10.1177/0962280220903763}.
}

\author{
Marco Geraci
}

\seealso{
\code{\link{nlmm}}
}


\keyword{fitting}
\keyword{control}
