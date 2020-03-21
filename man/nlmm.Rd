\name{nlmm}
\alias{nlmm}

\title{
Fitting Generalized Laplace Mixed-Effects Models
}
\description{
\code{nlmm} is used to fit mixed-effects models based on the generalized Laplace distribution.
}
\usage{
nlmm(fixed, random, group, covariance = "pdDiag", data = sys.frame(sys.parent()),
subset, weights = NULL, na.action = na.fail, control = list(), contrasts = NULL,
fit = TRUE)
}

\arguments{
  \item{fixed}{
an object of class \code{\link{formula}} for fixed effects: a symbolic description of the model to be fitted.
}
  \item{random}{
a one-sided formula of the form \code{~x1 + x2 + ... + xn} for random effects: a symbolic description of the model to be fitted.
}
  \item{group}{
grouping factor.
}
  \item{covariance}{
variance--covariance matrix of the random effects. Default is \code{pdDiag} (see details)
}
  \item{data}{
an optional data frame containing the variables named in \code{fixed}, \code{random}, \code{group}, and \code{weights}. By default the variables are taken from the environment from which \code{nlmm} is called.
}
  \item{subset}{
an optional vector specifying a subset of observations to be used in the fitting process.
}
  \item{weights}{
an optional \code{\link{varFunc}} object or one-sided formula describing the within-group heteroscedasticity structure. If given as a formula, it is used as the argument to \code{\link{varFixed}}, corresponding to fixed variance weights. See the documentation on \code{\link{varClasses}} in \code{nlme} package for a description of the available \code{\link{varFunc}} classes. Defaults to \code{NULL}, corresponding to homoscedastic within-group errors.
}
  \item{na.action}{
	a function that indicates what should happen when the
   data contain \code{NA}s.  The default action (\code{na.fail}) causes
   \code{nlmm} to print an error message and terminate if there are any
   incomplete observations.
}
  \item{control}{
list of control parameters of the fitting process. See \code{\link{nlmmControl}}.
}
  \item{contrasts}{
not yet implemented.
}
  \item{fit}{
logical flag. If \code{FALSE} the function returns a list of objects used for fitting.
}
}

\details{

The function fits a generalized Laplace mixed-effects model conditional on the covariates, as specified by the \code{formula} argument, and on random effects, as specified by the \code{random} argument. The predictor is assumed to be linear. The function maximizes the (log)likelihood of the generalized Laplace regression as proposed by Geraci and Farcomeni (2020). The likelihood is numerically integrated via Gaussian quadrature techniques. The optimization algorithm can be either \code{\link{optim}} (Nelder-Mead by default) or \code{\link{nlminb}}. See \code{\link{nlmmControl}} for more details.

By default, the function fits a mixed-effects model where both random effects and error term follow a generalized Laplace distribution (\code{\link{GenLaplace}}). This is a family of distributions that includes the normal and the Laplace distributions as special cases. Constrained fitting can be controlled via the arguments \code{alpha.index} and \code{alpha} in \code{\link{nlmmControl}}. For example, if \code{alpha.index = 0}, the model is Normal-Normal if \code{alpha = c(0,0)}, Normal-Laplace if \code{alpha = c(0,1)}, Laplace-Normal if \code{alpha = c(1,0)}, and Laplace-Laplace if \code{alpha = c(1,1)}. But any value of alpha between 0 (normal distribution) and 1 (Laplace distribution) is allowed.

Different standard types of positive--definite matrices for the random effects can be specified: \code{pdIdent} multiple of an identity; \code{pdCompSymm} compound symmetry structure (constant diagonal and constant off--diagonal elements); \code{pdDiag} diagonal; \code{pdSymm} general positive--definite matrix, with no additional structure.

Within-group heteroscedasticity can be modeled via the \code{weights} argument using \code{varClasses} in the \code{nlme} packages.
}

\value{
\code{nlmm} returns an object of \code{\link{class}} \code{nlmm}.

The function \code{summary} is used to obtain and print a summary of the results.

An object of class \code{nlmm} is a list containing the following components:

\item{theta}{a vector containing (in this order) fixed regression coefficients, parameters of the variance--covariance matrix of the random effects, shape parameter, scale parameter, and (optional) within-group variance function parameters. All parameters are unconstrained. See \code{\link{VarCorr.nlmm}} to extract the variance--covariance of the random effects from an "nlmm" object.}
\item{theta_x,theta_z}{partition of \code{theta}: fixed regression coefficients (\code{theta_x}) and unique variance--covariance unconstrained parameters (\code{theta_z}).}
\item{tau}{unconstrained shape parameter.}
\item{alpha}{constrained shape parameter.}
\item{phi}{unconstrained scale parameter.}
\item{sigma}{constrained scale parameter.}
\item{vf}{(fitted) within-group variance function of class \code{\link{varFunc}}. S3 methods (\code{\link{summary.varFunc}}, \code{\link{varWeights.varFunc}}, \code{\link{coef.varFunc}}) can be applied.}
\item{value}{negative log--likelihood.}
\item{call}{the matched call.}
\item{nn}{column names of \code{mmf}.}
\item{mm}{column names of \code{mmr}.}
\item{nobs}{the number of observations.}
\item{dim_theta}{the number of columns in \code{mmf} and \code{mmr}.}
\item{dim_theta_z}{the length of \code{theta_z}.}
\item{mmf}{the model matrix -- fixed effects.}
\item{mmr}{the model matrix -- random effects.}
\item{y}{the model response.}
\item{revOrder}{original order of observations (now ordered according to \code{group}).}
\item{group}{the grouping factor.}
\item{ngroups}{the number of groups.}
\item{InitialPar}{starting values for theta, included the fitted lme or lm object from where starting values have been taken.}
\item{control}{list of control parameters used for optimization (see \code{\link{nlmmControl}}).}
\item{cov_name}{class of variance-covariance matrix for the random effects.}
\item{mfArgs}{arguments for \code{\link{model.frame}} to return the full data frame.}
\item{sc}{model's distribution. "Generalized Laplace" if unconstrained estimation, or one of four special case for specific values of alpha ("Normal-Normal", "Normal-Laplace", "Laplace-Normal", "Laplace-Laplace").}

}

\references{
Geraci M (2017). Mixed-effects models using the normal and the Laplace distributions: A 2 x 2 convolution scheme for applied research. arXiv:1712.07216v1 [stat.ME]. URL: \url{https://arxiv.org/abs/1712.07216v1}.

Geraci, M. and Farcomeni A (2020). A family of linear mixed-effects models using the generalized Laplace distribution. Statistical Methods in Medical Research, \url{https://doi.org/10.1177/0962280220903763}.
}

\author{
Marco Geraci
}


\seealso{
\code{\link{summary.nlmm}, \link{fixef.nlmm}, \link{ranef.nlmm}, \link{VarCorr.nlmm}, \link{predict.nlmm}, \link{residuals.nlmm}}, \code{\link{nlmmControl}}
}

\examples{

data(rats)

nlmm(y ~ trt*time, random = ~ 1, group = id, data = rats,
control = nlmmControl(multistart = FALSE))
 
}
