################################################################################
# This is copyrighted material by Marco Geraci, University of South Carolina.
# All Rights Reserved.
################################################################################


#########################################################
### Fitting
#########################################################

nlmmControl <- function(method = "Nelder-Mead", nK = 8, multistart = TRUE, grid = c(0.001, 0.5, 0.999), alpha = c(0.5, 0.5), alpha.index = 9, lme = TRUE, lmeMethod = "REML", lmeOpt = "nlminb"){

if(length(alpha) != 2) stop("Provide starting values for alpha")
if(!alpha.index %in% c(0,1,2,9)) stop("alpha.index is one of c(0,1,2,9)")
if(any(alpha < 0 | alpha > 1)) stop("values for alpha must be between 0 and 1")
if(any(grid < 0 | grid > 1)) stop("values for alpha.grid must be between 0 and 1")

# 0 = both alpha fixed
# 1 = first alpha fixed
# 2 = second alpha fixed
# 9 = both alpha free

list(method = method, nK = nK, multistart = multistart, grid = grid, alpha = alpha, alpha.index = alpha.index, lme = lme, lmeMethod = lmeMethod, lmeOpt = lmeOpt)

}

loglik_nlmm <- function(theta, y, x, z, group, nK, P, Q, S, M, N, cov_name, vf){

dim_theta_w <- length(coef(vf))
if(dim_theta_w > 0){
	theta_w <- rev(rev(theta)[1:dim_theta_w])
	coef(vf) <- theta_w
}
w <- varWeights(vf)
y <- y*w
x <- sweep(x, 1, w, "*")
z <- sweep(z, 1, w, "*")

Y <- split(y, group)
Y <- lapply(Y, function(x) matrix(x, nrow = 1))
ni <- as.numeric(table(group))

theta_x <- theta[1:P]
theta_z <- theta[(P + 1):(P + S)]
tau <- theta[(P + S + 1) : (P + S + 2)]
alpha <- invlogit(tau) # inverse logit
sigma <- exp(theta[P + S + 3]) # inverse log
Sigma1 <- as.matrix(pdMat(value = theta_z, pdClass = cov_name, nam = 1:Q))*sigma^2
Sigma2 <- mapply(function(x, sigma) diag(sigma^2, x, x), ni, MoreArgs = list(sigma = sigma), SIMPLIFY = FALSE)

# quadrature
q1 <- gauss.quad.prob(nK, "gamma", alpha = min(1/alpha[1], 1e10), beta = 1)
q2 <- gauss.quad.prob(nK, "gamma", alpha = min(1/alpha[2], 1e10), beta = 1)
QUAD <- list(nodes = cbind(q1$nodes, q2$nodes), weights = cbind(q1$weights, q2$weights))

val <- C_ll(knots = QUAD$nodes, weights = QUAD$weights, beta = theta_x, Sigma1 = Sigma1, Sigma2 = Sigma2, Y = Y, x = x, z = z, M = M, N = N, ni = ni, P = P, Q = Q, K = nK)

# differential
val <- val + sum(log(w))

# negative integrated log-likelihood
return(-val)

}

loglik_alpha_nlmm <- function(theta, y, x, z, group, nK, P, Q, S, M, N, cov_name, vf, tau, index){

dim_theta_w <- length(coef(vf))
if(dim_theta_w > 0){
	theta_w <- rev(rev(theta)[1:dim_theta_w])
	coef(vf) <- theta_w
}

w <- varWeights(vf)
y <- y*w
x <- sweep(x, 1, w, "*")
z <- sweep(z, 1, w, "*")

Y <- split(y, group)
Y <- lapply(Y, function(x) matrix(x, nrow = 1))
ni <- as.numeric(table(group))

theta_x <- theta[1:P]
theta_z <- theta[(P + 1):(P + S)]

if(length(tau) == 1){
	if(index == 1){
		tau <- c(tau, theta[(P + S + 1)])
	}
	if(index == 2){
		tau <- c(theta[(P + S + 1)], tau)
	}
	sigma <- exp(theta[P + S + 2]) # inverse log
}

if(length(tau) == 2){
	sigma <- exp(theta[P + S + 1]) # inverse log
}

alpha <- invlogit(tau) # inverse logit

Sigma1 <- as.matrix(pdMat(value = theta_z, pdClass = cov_name, nam = 1:Q))*sigma^2
Sigma2 <- mapply(function(x, sigma) diag(sigma^2, x, x), ni, MoreArgs = list(sigma = sigma), SIMPLIFY = FALSE)

# quadrature
q1 <- gauss.quad.prob(nK, "gamma", alpha = min(1/alpha[1], 1e10), beta = 1)
q2 <- gauss.quad.prob(nK, "gamma", alpha = min(1/alpha[2], 1e10), beta = 1)
QUAD <- list(nodes = cbind(q1$nodes, q2$nodes), weights = cbind(q1$weights, q2$weights))

val <- C_ll(knots = QUAD$nodes, weights = QUAD$weights, beta = theta_x, Sigma1 = Sigma1, Sigma2 = Sigma2, Y = Y, x = x, z = z, M = M, N = N, ni = ni, P = P, Q = Q, K = nK)

# differential
val <- val + sum(log(w))

# negative integrated log-likelihood
return(-val)

}

nlmm <- function (fixed, random, group, covariance = "pdDiag", data = sys.frame(sys.parent()), subset, weights = NULL, na.action = na.fail, control = list(), contrasts = NULL, fit = TRUE){
Call <- match.call()
if (!is.data.frame(data)) 
	stop("`data' must be a data frame")
if (!inherits(fixed, "formula") || length(fixed) != 3) {
	stop("\nFixed-effects model must be a formula of the form \"y ~ x\"")
}
if (!inherits(random, "formula") || length(random) != 2) {
	stop("\nRandom-effects model must be a formula of the form \" ~ x\"")
}
groupFormula <- asOneSidedFormula(Call[["group"]])
group <- groupFormula[[2]]
mfArgs <- list(formula = asOneFormula(random, fixed, group, formula(varFunc(weights))), data = data, na.action = na.action)
if (!missing(subset)) {
	mfArgs[["subset"]] <- asOneSidedFormula(Call[["subset"]])[[2]]
}
mfArgs$drop.unused.levels <- TRUE
dataMix <- do.call("model.frame", mfArgs)
origOrder <- row.names(dataMix)
for (i in names(contrasts)) contrasts(dataMix[[i]]) = contrasts[[i]]
grp <- model.frame(groupFormula, dataMix)
ord <- order(unlist(grp, use.names = FALSE))
grp <- grp[ord, , drop = TRUE]
dataMix <- dataMix[ord, , drop = FALSE]
revOrder <- match(origOrder, row.names(dataMix))
ngroups <- length(unique(grp))
y <- eval(fixed[[2]], dataMix)
mmr <- model.frame(random, dataMix)
mmr <- model.matrix(random, data = mmr)
contr <- attr(mmr, "contr")
mmf <- model.frame(fixed, dataMix)
Terms <- attr(mmf, "terms")
auxContr <- lapply(mmf, function(el) if (inherits(el, "factor") && 
	length(levels(el)) > 1) 
	contrasts(el))
contr <- c(contr, auxContr[is.na(match(names(auxContr), names(contr)))])
contr <- contr[!unlist(lapply(contr, is.null))]
mmf <- model.matrix(fixed, data = mmf)
cov_name <- covariance
dim_theta <- integer(2)
dim_theta[1] <- ncol(mmf)
dim_theta[2] <- ncol(mmr)
dim_theta_z <- theta.z.dim(type = cov_name, n = dim_theta[2])
# weights
if(!is.null(weights)){
	vf <- varFunc(weights)
	vf <- Initialize(vf, data = dataMix)
	dim_theta_w <- length(coef(vf)) # 0 if coef is numeric(0)
} else {
	vf <- Initialize(varIdent(~1), data = dataMix)
	dim_theta_w <- 0
}
# optimization control parameters
if (is.null(names(control))) 
	control <- nlmmControl()
else {
	control_default <- nlmmControl()
	control_names <- intersect(names(control), names(control_default))
	control_default[control_names] <- control[control_names]
	control <- control_default
}
# determine if special case
sc <- "Generalized Laplace"
if(control$alpha.index == 0){
	if(all(control$alpha == 0)) sc <- "Normal-Normal"
	if(all(control$alpha == 1)) sc <- "Laplace-Laplace"
	if(control$alpha[1] == 0 & control$alpha[2] == 1) sc <- "Normal-Laplace"
	if(control$alpha[1] == 1 & control$alpha[2] == 0) sc <- "Laplace-Normal"

	#if(sc == "Normal-Normal") cat("Both alphas are fixed to 0. Fitting a standard linear mixed model with 'nlmm'. Consider using 'lme' instead", "\n")
}
# initialize
if(control$lme){
	reStruct <- list(group = nlme::pdMat(random, pdClass = cov_name))
	names(reStruct) <- as.character(group)
	lmeFit <- nlme::lme(fixed = fixed, random = reStruct, weights = weights, data = dataMix, method = control$lmeMethod, control = lmeControl(opt = control$lmeOpt))
	theta_x <- as.numeric(lmeFit$coefficients$fixed)
	theta_z <- as.numeric(coef(lmeFit[[1]]$reStruct)) # log-Cholesky scale
	theta_w <- as.numeric(coef(lmeFit[[1]]$varStruct, unconstrained = TRUE)) # numeric(0) if coef is NULL
	phi_0 <- log(lmeFit$sigma) # log scale
} else {
	lmFit <- lm(y ~ mmf - 1)
	theta_x <- lmFit$coef
	dim_theta_z <- theta.z.dim(type = cov_name, n = dim_theta[2])
	theta_z <- rep(0, dim_theta_z) # log-Cholesky scale
	theta_w <- rep(0, dim_theta_w) # numeric(0) if dim_theta_w is 0
	phi_0 <- log(mean(lmFit$residuals^2))/2 # log scale
}

alpha_0 <- control$alpha
tau <- logit(alpha_0, omega = 0.001) # logit scale
tau_0 <- if(control$alpha.index != 9) tau[-control$alpha.index] else tau
theta_0 <- c(theta_x, theta_z, tau_0, phi_0, theta_w)
FIT_ARGS <- list(theta = theta_0, y = y, x = mmf, z = mmr, group = grp, nK = control$nK, P = dim_theta[1], Q = dim_theta[2], S = dim_theta_z, M = ngroups, N = length(y), cov_name = cov_name, vf = vf)

if(!fit){
	return(FIT_ARGS)
}

if(sc == "Normal-Normal"){
	cat("Both alphas are fixed to 0. Fitting a standard linear mixed model with 'lme(..., method = 'ML')'", "\n")
	reStruct <- list(group = nlme::pdMat(random, pdClass = cov_name))
	names(reStruct) <- as.character(group)
	lmeFit <- nlme::lme(fixed = fixed, random = reStruct, weights = weights, data = dataMix, method = "ML", control = lmeControl(opt = control$lmeOpt))
	ans <- lme2nlmm(x = lmeFit, Call = Call, mmf = mmf, mmr = mmr, y = y, revOrder = revOrder, vf = vf, contr = contr, grp = grp, control = control, cov_name = cov_name, mfArgs = mfArgs)
	return(ans)

}

if(control$multistart & control$alpha.index == 0){
	control$multistart <- FALSE
	cat("Both alphas are fixed. Ignoring multistart", "\n")
}

if(control$multistart){

	FLAG <- control$alpha.index %in% c(1, 2)

	if(FLAG){
		alpha.grid <- as.matrix(control$grid)
	} else {
		alpha.grid <- as.matrix(expand.grid(control$grid, control$grid))
	}
	tau.grid <- logit(alpha.grid, omega = 0.001) # logit scale
	K <- nrow(tau.grid)
	tmp <- list()
	pb <- txtProgressBar(title = "Multistart for nlmm", label = "0% done", min = 0, max = 100)
	for(k in 1:K){
		info <- sprintf("%d%% done", round((k/K)*100)) 
		setTxtProgressBar(pb, k/K*100, label = info)
		
		FIT_ARGS$theta <- c(theta_x, theta_z, tau.grid[k,], phi_0, theta_w)
		if(FLAG){
			FIT_ARGS$tau <- tau[control$alpha.index]
			FIT_ARGS$index <- control$alpha.index
			tmp[[k]] <- do.call(optim, args = c(list(fn = loglik_alpha_nlmm, par = FIT_ARGS$theta, method = control$method, control = list(trace = 0)), FIT_ARGS[-c(match(c("theta"), names(FIT_ARGS)))]))
		} else {
			tmp[[k]] <- do.call(optim, args = c(list(fn = loglik_nlmm, par = FIT_ARGS$theta, method = control$method, control = list(trace = 0)), FIT_ARGS[-c(match(c("theta"), names(FIT_ARGS)))]))
		}
	}
	close(pb)
	sel <- which.min(sapply(tmp, function(x) x$value))[1]
	fit <- tmp[[sel]]
	theta_0 <- c(theta_x, theta_z, tau.grid[sel,], phi_0, theta_w)
} else {

	if(control$alpha.index == 9){
		fit <- do.call(optim, args = c(list(fn = loglik_nlmm, par = FIT_ARGS$theta, method = control$method, control = list(trace = 0)), FIT_ARGS[-c(match(c("theta"), names(FIT_ARGS)))]))
	} else {
		sel <- if(control$alpha.index == 0) 1:2 else control$alpha.index
		FIT_ARGS$tau <- tau[sel]
		FIT_ARGS$index <- control$alpha.index
		fit <- do.call(optim, args = c(list(fn = loglik_alpha_nlmm, par = FIT_ARGS$theta, method = control$method, control = list(trace = 0)), FIT_ARGS[-c(match(c("theta"), names(FIT_ARGS)))]))
	}
}

nn <- colnames(mmf)
mm <- colnames(mmr)

fit$theta_x <- fit$par[1:dim_theta[1]]
names(fit$theta_x) <- nn
fit$theta_z <- fit$par[(dim_theta[1] + 1):(dim_theta[1] + dim_theta_z)]

if(control$alpha.index == 0){
	fit$alpha <- control$alpha
	fit$tau <- logit(fit$alpha, omega = 0.001)
	fit$phi <- fit$par[dim_theta[1] + dim_theta_z + 1]
	fit$sigma <- exp(fit$phi)
}

if(control$alpha.index == 1){
	fit$alpha <- c(control$alpha[1], invlogit(fit$par[dim_theta[1] + dim_theta_z + 1]))
	fit$tau <- logit(fit$alpha, omega = 0.001)
	fit$phi <- fit$par[dim_theta[1] + dim_theta_z + 2]
	fit$sigma <- exp(fit$phi)
}

if(control$alpha.index == 2){
	fit$alpha <- c(invlogit(fit$par[dim_theta[1] + dim_theta_z + 1]), control$alpha[2])
	fit$tau <- logit(fit$alpha, omega = 0.001)
	fit$phi <- fit$par[dim_theta[1] + dim_theta_z + 2]
	fit$sigma <- exp(fit$phi)
}

if(control$alpha.index == 9){
	fit$alpha <- invlogit(fit$par[(dim_theta[1] + dim_theta_z + 1) : (dim_theta[1] + dim_theta_z + 2)])
	fit$tau <- logit(fit$alpha, omega = 0.001)
	fit$phi <- fit$par[dim_theta[1] + dim_theta_z + 3]
	fit$sigma <- exp(fit$phi)
}

if(dim_theta_w > 0){
	theta_w <- rev(rev(fit$par)[1:dim_theta_w])
}

coef(vf) <- theta_w

fit$call <- Call
fit$nn <- nn
fit$mm <- mm
fit$nobs <- length(y)
fit$dim_theta <- dim_theta
fit$dim_theta_z <- dim_theta_z
fit$edf <- fit$dim_theta[1] + fit$dim_theta_z
fit$rdf <- fit$nobs - fit$edf
fit$df <- dim_theta[1] + dim_theta_z + 1
fit$mmf <- mmf
fit$mmr <- mmr
fit$y <- y
fit$revOrder <- revOrder
fit$vf <- vf
fit$contrasts <- contr
fit$group <- grp
fit$ngroups <- ngroups
fit$InitialPar <- list(fit = if(control$lme) lmeFit else lmFit, theta = theta_0)
fit$control <- control
fit$cov_name <- cov_name
fit$mfArgs <- mfArgs
fit$sc <- sc
class(fit) <- "nlmm"
fit
}

#########################################################
### LRT
#########################################################

wchibarsq <- function(V){

	safeseq <- function(from = 1L, to = 1L, by = 1L,...){
		disc <- by*(from-to)
		if(disc > 0){
			vector(class(disc), 0L)
		} else seq(from = from, to = to, by = by, ...)
	}
	
stopifnot(is.matrix(V) && diff(dim(V)) == 0L && nrow(V) >  0L)
n <- nrow(V)
P <- function(idx){
	pmvnorm(rep(0, n - i), rep(Inf, n - i), sigma = solve(V[-idx, -idx, drop = FALSE])) * pmvnorm(rep(0, i), rep(Inf, i), sigma = V[idx, idx, drop = FALSE] - V[idx, -idx, drop = FALSE] %*% solve(V[-idx, -idx, drop = FALSE], V[-idx, idx, drop = FALSE]))
}
ans <- numeric(n + 1L)
ans[1] <- pmvnorm(rep(0, n), rep(Inf, n), sigma = solve(V))[[1]]
ans[n + 1L] <- pmvnorm(rep(0, n), rep(Inf, n), sigma = V)[[1]]
for (i in safeseq(1L, n - 1L, by = 1L)) ans[i + 1] = sum(combn(n, i, P))
ans
}

pchibarsq <- function (q, w, lower.tail = TRUE, log.p = FALSE){
n <- length(w)
    ans <- pchisq(q, 0, lower.tail = FALSE) * w[1L] + pchisq(q, n - 1, lower.tail = FALSE) * w[n]
if(n > 2) {
	for (i in seq_len(n - 2)) {
		ans = ans + pchisq(q, i, lower.tail = FALSE) * w[i + 1L]
	}
}
ans[q <= 0] <- 1
ans <- if(isTRUE(lower.tail)) {1 - ans} else ans
if(isTRUE(log.p)) ans <- log(ans)
ans
}

Vchibarsq <- function(object){

alpha.index <- object$control$alpha.index
if(!alpha.index %in% c(0, 1, 2)) stop("This fitted model is not constrained. 'alpha.index' must be either 0, 1, or 2")

FIT_ARGS <- list(theta = object$InitialPar$theta, y = object$y, x = object$mmf, z = object$mmr, group = object$group, nK = object$control$nK, P = object$dim_theta[1], Q = object$dim_theta[2], S = object$dim_theta_z, M = object$ngroups, N = length(object$y), cov_name = object$cov_name, vf = object$vf)
# Modify FIT_ARGS$theta as it was an unconstrained fit
dim_theta_w <- length(coef(object$vf))
if(dim_theta_w > 0){
	theta_w <- coef(object$vf, unconstrained = TRUE)
} else {
	theta_w <- NULL
}
FIT_ARGS$theta <- c(object$theta_x, object$theta_z, object$tau, object$phi, theta_w)
alpha <- object$alpha 

if(alpha.index == 0){

       f <- function(tau_tilde, MoreArgs, alpha){
             P <- MoreArgs$P
             S <- MoreArgs$S
             alpha_tilde <- invlogit(tau_tilde)
             tau <- MoreArgs$theta[(P + S + 1):(P + S + 2)]
             tau[1] <- if(alpha[1] == 1) logit(1 - alpha_tilde[1], omega = 1e-3) else if(alpha[1] == 0) tau_tilde[1] else tau[1]
             tau[2] <- if(alpha[2] == 1) logit(1 - alpha_tilde[2], omega = 1e-3) else if(alpha[2] == 0) tau_tilde[2] else tau[2]
             MoreArgs$theta[(P + S + 1):(P + S + 2)] <- tau
             do.call(loglik_nlmm, args = MoreArgs)
       }

val <- hessian(func = f, x = logit(c(0, 0), omega = 1e-3), method = "Richardson", MoreArgs = FIT_ARGS, alpha = alpha)

}

if(alpha.index == 1){

       f <- function(tau_tilde, MoreArgs, alpha){
             P <- MoreArgs$P
             S <- MoreArgs$S
             alpha_tilde <- invlogit(tau_tilde)
             tau <- if(alpha[1] == 1) logit(1 - alpha_tilde, omega = 1e-3) else tau_tilde
             MoreArgs$theta[(P + S + 1)] <- tau
             do.call(loglik_nlmm, args = MoreArgs)
       }
       
val <- hessian(func = f, x = logit(0, omega = 1e-3), method = "Richardson", MoreArgs = FIT_ARGS, alpha = alpha)

}

if(alpha.index == 2){

       f <- function(tau_tilde, MoreArgs, alpha){
             P <- MoreArgs$P
             S <- MoreArgs$S
             alpha_tilde <- invlogit(tau_tilde)
             tau <- if(alpha[2] == 1) logit(1 - alpha_tilde, omega = 1e-3) else tau_tilde
             MoreArgs$theta[(P + S + 2)] <- tau
             do.call(loglik_nlmm, args = MoreArgs)
       }
       
val <- hessian(func = f, x = logit(0, omega = 1e-3), method = "Richardson", MoreArgs = FIT_ARGS, alpha = alpha)

}


if(!is.positive.definite(val)){
	val <- make.positive.definite(val)
}

return(val)


}

lrt_nlmm <- function(object0, object1){

# check object0 is constrained
if(!object0$control$alpha.index %in% c(0, 1, 2)) stop("'object0' must be a constrained 'nlmm' object")
# check object1 is unconstrained
if(object1$control$alpha.index != 9) stop("'object1' must be an unconstrained 'nlmm' object")
# determine if chi or chibar
alpha <- object0$control$alpha
index <- object0$control$alpha.index
sel <- if(index == 0) 1:2 else index
FLAG <- any(alpha[sel] %in% c(0, 1))

if(FLAG){
V <- Vchibarsq(object0)
w <- wchibarsq(V)
statistic <- -2*(object1$val - object0$val)
pval <- pchibarsq(q = statistic, w = w, lower.tail = FALSE, log.p = FALSE)
} else {
V <- NULL
w <- if(index == 0) 2 else 1
statistic <- -2*(object1$val - object0$val)
pval <- pchisq(q = statistic, df = w, lower.tail = FALSE, log.p = FALSE)
}

if(statistic < 0){
	statistic <- 0
	pval <- 1
	warning("Negative LRT: possible local optimum")
}

ans <- list(statistic = statistic, p.value = pval, df = w, V = V, alpha = alpha, alpha.index = index, chibar = FLAG)
class(ans) <- "lrt_nlmm"

ans
}

#########################################################
### Methods
#########################################################

fixef.nlmm <- function(object){

return(object$theta_x)

}

logLik.nlmm <- function(object, ...){

	return(-object$value)

}

predict.nlmm <- function(object, level = 0, ...){

group <- object$group
M <- object$ngroups
q <- object$dim_theta[2]

FXD <- object$mmf %*% matrix(object$theta_x)

if(level == 1) {
	RE <- ranef(object)
	mmr.l <- split(object$mmr, group)
	RE.l <- split(RE, unique(group))
	RND <- NULL
	for (i in 1:M) {
		RND <- rbind(RND, matrix(as.numeric(mmr.l[[i]]), ncol = q) %*% matrix(as.numeric(RE.l[[i]]), nrow = q))
	}
}

if(level == 0) {
	ans <- FXD[object$revOrder, ]
}

if (level == 1) {
	ans <- FXD + RND
	ans <- ans[object$revOrder, ]
}

return(ans)

}

ranef.nlmm <- function(object, ...){

group <- object$group
ni <- table(group)
theta_z <- object$theta_z
alpha <- invlogit(object$tau)
names(alpha) <- c("Random effects", "Error")
q <- object$dim_theta[2]
 
Sigma1 <- VarCorr(object)
w <- varWeights(object$vf)
vv <- 1/w^2
vv <- split(vv, group)

if(object$sc == "Normal-Normal"){
	Sigma2 <- lapply(vv, function(x, sigma) diag(x = x, length(x), length(x))*sigma^2, sigma = object$sigma)
} else {
	Sigma2 <- lapply(vv, function(x, sigma, alpha) diag(x = x, length(x), length(x))*sigma^2/alpha, sigma = object$sigma, alpha = alpha["Error"])
}

RES <- split(object$y - object$mmf %*% matrix(object$theta_x), group)
mmrList <- split(object$mmr, group)

BLPu <- vector("list", object$ngroups)

for(i in 1:object$ngroups){
	Zi <- matrix(mmrList[[i]], nrow = ni[i])
	Psi <- Zi %*% Sigma1 %*% t(Zi) + Sigma2[[i]]
	BLPu[[i]] <- Sigma1 %*% t(Zi) %*% solve(Psi) %*% matrix(RES[[i]])
}
ans <- data.frame(matrix(unlist(BLPu), ncol = q, byrow = TRUE))
rownames(ans) <- unique(group)
colnames(ans) <- object$mm
return(ans)
}

residuals.nlmm <- function(object, level = 0, ...){

	object$y[object$revOrder] - predict(object, level = level)

}

summary.nlmm <- function(object, alpha = 0.05, ...){


	g <- function(x){
		exp(-x)/(1+ exp(-x))^2
	}


V <- vcov(object)

dim_theta <- object$dim_theta
dim_theta_z <- object$dim_theta_z
p <- dim_theta[1]
index <- object$control$alpha.index

theta_x <- object$theta_x
SE_theta_x <- sqrt(diag(V)[1:p])
lower_theta_x <- theta_x - SE_theta_x * qnorm(1 - alpha/2, 0, 1)
upper_theta_x <- theta_x + SE_theta_x * qnorm(1 - alpha/2, 0, 1)

tau <- object$tau
if(index == 0){
SE_tau <- c(NA, NA)
}
if(index == 1){
SE_tau <- c(NA, sqrt(diag(V)[(dim_theta[1] + dim_theta_z + 1)]))
}
if(index == 2){
SE_tau <- c(sqrt(diag(V)[(dim_theta[1] + dim_theta_z + 1)]), NA)
}
if(index == 9){
SE_tau <- sqrt(diag(V)[(dim_theta[1] + dim_theta_z + 1) : (dim_theta[1] + dim_theta_z + 2)])
}
SE_alpha <- sqrt(SE_tau^2*g(tau)^2)
lower_tau <- tau - SE_tau * qnorm(1 - alpha/2, 0, 1)
upper_tau <- tau + SE_tau * qnorm(1 - alpha/2, 0, 1)

alpha <- object$alpha
lower_alpha <- invlogit(lower_tau)
upper_alpha <- invlogit(upper_tau)

tTable <- data.frame(c(theta_x, alpha), c(SE_theta_x, SE_alpha), c(lower_theta_x, lower_alpha), c(upper_theta_x, upper_alpha))
names(tTable) <- c("Estimate", "Std.Err", "Lower", "Upper")
rownames(tTable) <- c(object$nn, "Random effects", "Error")
object$tTable <- tTable
class(object) <- "summary.nlmm"
return(object)
}

VarCorr.nlmm <- function(x, sigma = NULL, ...){

theta_z <- x$theta_z
alpha <- invlogit(x$tau)
names(alpha) <- c("Random effects", "Error")

Sigma1 <- as.matrix(pdMat(value = theta_z, pdClass = x$cov_name, nam = x$mm))*x$sigma^2

if(x$sc != "Normal-Normal") Sigma1 <- Sigma1/alpha['Random effects']

return(Sigma1)

}

vcov.nlmm <- function(object, ...){

if(object$sc == "Normal-Normal"){
	ans <- matrix(0, length(object$par), length(object$par))
	ans <- as.matrix(bdiag(object$InitialPar$fit$varFix, object$InitialPar$fit$apVar))
	return(ans)
}

FIT_ARGS <- list(theta = object$InitialPar$theta, y = object$y, x = object$mmf, z = object$mmr, group = object$group, nK = object$control$nK, P = object$dim_theta[1], Q = object$dim_theta[2], S = object$dim_theta_z, M = object$ngroups, N = length(object$y), cov_name = object$cov_name, vf = object$vf)

index <- object$control$alpha.index
FLAG <- index %in% c(0, 1, 2)

if(FLAG){
	sel <- if(index == 0) 1:2 else index
	FIT_ARGS$tau <- logit(object$control$alpha, omega = 1e-3)[sel]
	FIT_ARGS$index <- index
}

if(FLAG){
	f <- function(theta, args){
		args$theta <- theta
		do.call(loglik_alpha_nlmm, args = args)
	}
} else {
	f <- function(theta, args){
		args$theta <- theta
		do.call(loglik_nlmm, args = args)
	}
}

H <- hessian(func = f, x = object$par, method = "Richardson", args = FIT_ARGS[-c(match(c("theta"), names(FIT_ARGS)))])
ans <- MASS::ginv(H)

if(!is.positive.definite(ans)){
	ans <- make.positive.definite(ans)
}

return(ans)

}

print.nlmm <- function (x, digits = max(3, getOption("digits") - 3), ...){

theta_x <- x$theta_x
theta_z <- x$theta_z
alpha <- invlogit(x$tau)
names(alpha) <- names(x$alpha) <- c("Random effects", "Error")

Sigma1 <- as.matrix(pdMat(value = theta_z, pdClass = x$cov_name, nam = x$mm))*x$sigma^2
Sigma2 <- x$sigma^2
if(x$sc != "Normal-Normal"){
	Sigma1 <- Sigma1/alpha['Random effects']
	Sigma2 <- Sigma2/alpha['Error']
}

cat("Call: ")
dput(x$call)
cat("\n")
cat("Generalized Laplace mixed-effects model", "\n")
if(x$sc != "Generalized Laplace") cat("Special case:", x$sc, "\n")
cat("\n")
cat("Alpha:\n")
print.default(format(x$alpha, digits = digits), print.gap = 2, 
	quote = FALSE)
cat("\n")
cat("Fixed effects:\n")
print.default(format(theta_x, digits = digits), print.gap = 2, 
	quote = FALSE)
cat("\n")
cat("Covariance matrix of the random effects:\n")
print.default(format(as.matrix(Sigma1), digits = digits), 
	quote = FALSE)
cat("\n")
cat(paste("Residual variance: ", format(Sigma2, digits = digits)), 
	"\n")
print(summary(x$vf))
cat("\n")
cat(paste("Log-likelihood:", format(-x$value, digits = digits), 
	"\n"))
cat(paste("\nNumber of observations:", x$nobs, "\n"))
cat(paste("Number of groups:", x$ngroups, "\n"))
invisible(x)
}

print.summary.nlmm <- function(x, digits = max(3, getOption("digits") - 3), ...){

theta_x <- x$theta_x
theta_z <- x$theta_z
alpha <- invlogit(x$tau)
names(alpha) <- c("Random effects", "Error")

Sigma1 <- as.matrix(pdMat(value = theta_z, pdClass = x$cov_name, nam = x$mm))*x$sigma^2
Sigma2 <- x$sigma^2
if(x$sc != "Normal-Normal"){
	Sigma1 <- Sigma1/alpha['Random effects']
	Sigma2 <- Sigma2/alpha['Error']
}

tTable <- x$tTable
xx <- tTable[1:x$dim_theta[1], , drop = FALSE]
yy <- tTable[-c(1:x$dim_theta[1]), , drop = FALSE]

cat("Call: ")
dput(x$call)
cat("\n")
cat("Generalized Laplace mixed-effects model", "\n")
if(x$sc != "Generalized Laplace") cat("Special case:", x$sc, "\n")
cat("\n")
cat("Alpha:\n")
print(yy, ...)
cat("\n")
cat("Fixed effects:\n")
print(xx, ...)
cat("\n")
cat("Covariance matrix of the random effects:\n")
print.default(format(as.matrix(Sigma1), digits = digits), 
	quote = FALSE)
cat("\n")
cat(paste("Residual variance: ", format(Sigma2, digits = digits)), 
	"\n")
print(summary(x$vf))
cat("\n")
cat(paste("Log-likelihood:", format(-x$value, digits = digits), 
	"\n"))
cat(paste("\nNumber of observations:", x$nobs, "\n"))
cat(paste("Number of groups:", x$ngroups, "\n"))
invisible(x)
}

print.lrt_nlmm <- function(x, digits = max(3, getOption("digits") - 3), ...){

txt <- c("random effects", "error")
type <- if(x$chibar) "Chi-bar-squared" else "Chi-squared"

cat("Likelihood ratio test for constrained", "\n", "generalized Laplace mixed-effects models", "\n")
cat("\n")
if(x$alpha.index == 0){
cat("H0: alpha random effects is equal to", x$alpha[1], "and alpha error is equal to", x$alpha[2], "\n")
} else {
cat("H0: alpha", txt[x$alpha.index], "is equal to", x$alpha[x$alpha.index], "\n")
}
cat("LRT statistic =", round(x$statistic, digits = digits), "\n")
cat("df or weights =", round(x$df, digits = 2), "\n")
cat("p-value (", type, ") = ", round(x$p.value, digits = digits), "\n", sep = "")
}

#########################################################
### Auxiliary functions
#########################################################

# convert from lme to nlmm object

lme2nlmm <- function(x, Call, mmf, mmr, y, revOrder, vf, contr, grp, control, cov_name, mfArgs){

	theta_x <- as.numeric(x$coefficients$fixed)
	theta_z <- as.numeric(coef(x[[1]]$reStruct)) # log-Cholesky scale
	theta_w <- as.numeric(coef(x[[1]]$varStruct)) # numeric(0) if coef is NULL
	sigma <- x$sigma
	phi <- log(sigma) # log scale
	alpha <- c(0, 0)
	#tau <- c(Inf, Inf)
	tau <- logit(alpha, omega = 0.001)
	theta <- c(theta_x, theta_z, phi)
	coef(vf) <- theta_w
	nn <- colnames(mmf)
	mm <- colnames(mmr)
	dim_theta <- integer(2)
	dim_theta[1] <- ncol(mmf)
	dim_theta[2] <- ncol(mmr)
	dim_theta_z <- length(theta_z)
	ngroups <- length(unique(grp))

	fit <- list()
	
	fit$value <- -logLik(x, REML = FALSE)
	fit$par <- theta
	fit$theta_x <- theta_x
	names(fit$theta_x) <- nn
	fit$theta_z <- theta_z
	fit$alpha <- alpha
	fit$tau <- tau
	fit$phi <- phi
	fit$sigma <- sigma

	fit$call <- Call
	fit$nn <- nn
	fit$mm <- mm
	fit$nobs <- length(y)
	fit$dim_theta <- dim_theta
	fit$dim_theta_z <- dim_theta_z
	fit$edf <- fit$dim_theta[1] + fit$dim_theta_z
	fit$rdf <- fit$nobs - fit$edf
	fit$df <- dim_theta[1] + dim_theta_z + 1
	fit$mmf <- mmf
	fit$mmr <- mmr
	fit$y <- y
	fit$revOrder <- revOrder
	fit$vf <- vf
	fit$contrasts <- contr
	fit$group <- grp
	fit$ngroups <- ngroups
	fit$InitialPar <- list(fit = x, theta = theta)
	fit$control <- control
	fit$cov_name <- cov_name
	fit$mfArgs <- mfArgs
	fit$sc <- "Normal-Normal"
	class(fit) <- "nlmm"
	fit

}

#########################################################
### Distributions
#########################################################

### Generalized Laplace

# Density of the (symmetric) generalized Laplace

dgl <- function(x, mu = 0, sigma = 1, shape = 1, log = FALSE){

# symmetric generalized Laplace
# Kotz et al p.190
# mu = location
# sigma = scale
# variance = shape*sigma^2

n <- length(x)
if(length(mu) == 1) mu <- rep(mu, n)
if(length(sigma) == 1) sigma <- rep(sigma, n)
if(length(shape) == 1) shape <- rep(shape, n)
if(length(mu) != n) stop("'mu' must be a vector of 1 or of the same length of x")
if(length(sigma) != n) stop("'sigma' must be a vector of 1 or of the same length of x")
if(length(shape) != n) stop("'shape' must be a vector of 1 or of the same length of x")

p <- shape - 1/2

val1 <- sqrt(2)/(sigma^(p + 1)*gamma(shape)*sqrt(pi))
val2 <- (abs(x - mu)/sqrt(2))^p
val3 <- besselK(sqrt(2)*abs(x - mu)/sigma, nu = p, expon.scaled = FALSE)

val <- val1*val2*val3
if(log) val <- log(val)

return(val)
}

# CDF of the (symmetric) generalized Laplace

pgl <- function(x, mu = 0, sigma = 1, shape = 1, lower.tail = TRUE, log.p = FALSE){

# symmetric generalized Laplace
# Kotz et al p.190
# mu = location
# sigma = scale
# variance = shape*sigma^2

n <- length(x)
val <- rep(NA, n)
if(length(mu) == 1) mu <- rep(mu, n)
if(length(sigma) == 1) sigma <- rep(sigma, n)
if(length(shape) == 1) shape <- rep(shape, n)
if(length(mu) != n) stop("'mu' must be a vector of 1 or of the same length of x")
if(length(sigma) != n) stop("'sigma' must be a vector of 1 or of the same length of x")
if(length(shape) != n) stop("'shape' must be a vector of 1 or of the same length of x")
	
if(lower.tail){
	for(i in 1:n){
		val[i] <- integrate(dgl, lower = -Inf, upper = x[i], mu = mu[i], sigma = sigma[i], shape = shape[i])$value
	}
} else {
	for(i in 1:n){
		val[i] <- integrate(dgl, lower = x[i], upper = Inf, mu = mu[i], sigma = sigma[i], shape = shape[i])$value
	}
}

if(log.p) val <- log(val)

return(val)
}

# quantile of the (symmetric) generalized Laplace

qgl <- function(p, mu = 0, sigma = 1, shape = 1, lower.tail = TRUE, log.p = FALSE){

# symmetric generalized Laplace
# Kotz et al p.190
# mu = location
# sigma = scale
# variance = shape*sigma^2

n <- length(p)

if(length(mu) == 1) mu <- rep(mu, n)
if(length(sigma) == 1) sigma <- rep(sigma, n)
if(length(shape) == 1) shape <- rep(shape, n)
if(length(mu) != n) stop("'mu' must be a vector of 1 or of the same length of x")
if(length(sigma) != n) stop("'sigma' must be a vector of 1 or of the same length of x")
if(length(shape) != n) stop("'shape' must be a vector of 1 or of the same length of x")

if(log.p) p <- exp(p)

f <- function(x, p, mu, sigma, shape){
	p - pgl(x, mu = mu, sigma = sigma, shape = shape)
}

V <- shape*sigma^2
val <- rep(NA, n)
for(i in 1:n){
	 tmp <- try(uniroot(f, p = p[i], mu = mu[i], sigma = sigma[i], shape = shape[i], interval = c(mu[i] - 20*sqrt(V[i]), mu[i] + 20*sqrt(V[i]))), silent = TRUE)
	 if(!inherits(tmp, "try-error")) val[i] <- tmp$root
}

val[p == 0.5] <- mu[p == 0.5]

if(!lower.tail) val <- -val

return(val)
}

# Random generation for the (symmetric) generalized Laplace

rgl <- function(n, mu = 0, sigma = 1, shape = 1){

# symmetric generalized Laplace
# Kotz et al p.190
# mu = location
# sigma = scale
# variance = shape*sigma^2

W <- rgamma(n, shape = shape, scale = 1)
N <- rnorm(n, sd = sigma)
val <- mu + sqrt(W)*N

attr(val, "scale") <- W

return(val)

}


### Multivariate generalized Laplace

# Density of the (centered) asymmetric multivariate generalized Laplace

dmgl <- function(x, mu = rep(0, n), sigma = diag(n), shape = 1, log = FALSE){

# asymmetric multivariate generalized Laplace (symmetric if mu = 0)
# Kozubowski et al et (2013, Journal of Multivariate Analysis)
# mu = symmetry
# sigma = scale


Q <- function(x, sigma){
	x <- as.matrix(x)
	val <- sqrt(crossprod(x, sigma) %*% x)
	return(val)
}

C <- function(mu, sigma){
	mu <- as.matrix(mu)
	val <- sqrt(2 + crossprod(mu, sigma) %*% mu)
	return(val)
}

n <- length(x)
mu <- as.matrix(mu)
sigma <- as.matrix(sigma)
x <- as.matrix(x)
p <- shape - n/2
k <- sqrt(det(sigma))
sigma <- solve(sigma)

val1 <- 2*exp(crossprod(mu, sigma) %*% x)/((2*pi)^(n/2)*gamma(shape)*k)
val2 <- (Q(x, sigma)/C(mu, sigma))^p
val3 <- besselK(Q(x, sigma)*C(mu, sigma), nu = p, expon.scaled = FALSE)

val <- val1*val2*val3
if(log) val <- log(val)

attr(val, "terms") <- c(val1, val2, val3)
return(val)

}

# Random generation for the asymmetric multivariate generalized Laplace

rmgl <- function(n, mu, sigma, shape = 1){

# asymmetric multivariate generalized Laplace (symmetric if mu = 0)
# Kozubowski et al et (2013, Journal of Multivariate Analysis)
# mu = symmetry
# sigma = scale

    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
        check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mu) != nrow(sigma)) 
        stop("mu and sigma have non-conforming size")

sigma <- as.matrix(sigma)
W <- rgamma(n, shape = shape, scale = 1)
N <- mvtnorm::rmvnorm(n, sigma = sigma)
mu <- matrix(mu, nrow = 1)

val <- matrix(kronecker(mu, W), ncol = length(mu)) + sweep(N, 1, sqrt(W), "*")

attr(val, "scale") <- W

return(val)
}


################################################################################
# R code to generate data in the simulation study 'A family of linear mixed-effects models using the generalized Laplace distribution' by Geraci and Farcomeni
################################################################################

slash <- function(n, mean, sd){
	val <- rnorm(n, mean = mean, sd = sd)/runif(n)
	return(val)
}

generate.dist <- function(fun, n, q, sigma, shape){

if(q == 1){
	if(fun == "norm"){
		fun <- rnorm
		return(list(fun = fun, n = n, mean = 0, sd = sigma))
	}

	if(fun == "laplace"){
		fun <- rl
		return(list(fun = fun, n = n, mu = 0, sigma = sigma))
	}

	if(fun == "genlaplace"){
		fun <- rmgl
		return(list(fun = fun, n = n, mu = rep(0, 1), sigma = as.matrix(sigma), shape = shape))
	}
	
	if(fun == "chisq"){
		fun <- rchisq
		return(list(fun = fun, n = n, df = sigma))
	}

	if(fun == "t"){
		fun <- rt
		return(list(fun = fun, n = n, df = shape))
	}

	if(fun == "slash"){
		fun <- slash
		return(list(fun = fun, n = n, mean = 0, sd = sigma))
	}


}

if(q > 1){
	if(fun == "norm"){
		fun <- rmvnorm
		return(list(fun = fun, n = n, mean = rep(0, q), sigma = sigma))
	}

	if(fun == "laplace"){
		fun <- rmal
		return(list(fun = fun, n = n, mu = rep(0, q), sigma = sigma))
	}
	
	if(fun == "genlaplace"){
		fun <- rmgl
		return(list(fun = fun, n = n, mu = rep(0, q), sigma = sigma, shape = shape))
	}

	if(fun == "t"){
		fun <- rmvt
		return(list(fun = fun, n = n, sigma = sigma, df = shape))
	}

}

}

generate.design <- function(n, M, fixed = FALSE){

# M groups
# n measurements in each group

N <- n*M
if(fixed){
	x <- rep(1:n, M)
	z <- rbinom(n = N, size = 1, prob = 0.5)
} else {
	delta <- rnorm(M, 0, 1)
	zeta <- rnorm(N, 0, 1)
	x <- rep(delta, each = n) + zeta
	z <- rbinom(n = N, size = 1, prob = 0.5)
}

X <- cbind(1, x, z)
colnames(X) <- c("intercept","x","z")
X
}

generate.data <- function(R, n, M, sigma_1 = NULL, sigma_2 = NULL, shape_1 = NULL, shape_2 = NULL, dist.u, dist.e, beta, gamma, fixed = FALSE, seed = round(runif(1,1,1000))){

# M groups
# n measurements in each group

set.seed(seed)
N <- n*M
id <- rep(1:M, each = n)
beta <- matrix(beta)
gamma <- matrix(gamma)

q_1 <- nrow(as.matrix(sigma_1))
q_2 <- nrow(as.matrix(sigma_2))

if(q_1 > 2) stop("max q = 2 for random effects")

par.u <- generate.dist(fun = dist.u, n = M, q = q_1, sigma = sigma_1, shape = shape_1)
par.e <- generate.dist(fun = dist.e, n = M, q = q_2, sigma = sigma_2, shape = shape_2)

U <- replicate(R, do.call(par.u$fun, args = par.u[-1]))
e <- replicate(R, do.call(par.e$fun, args = par.e[-1]))

if(R == 1){
u <- if(q_1 == 2) rep(U[,1,1], each = n) else rep(U, each = n)
v <- if(q_1 == 2) rep(U[,2,1], each = n) else rep(0,N)
e <- as.vector(t(e[,,1]))
}

if(R > 1){
u <- if(q_1 == 2) apply(U[,1,], 2, function(x, n) rep(x, each = n), n = n) else apply(U, 2, function(x, n) rep(x, each = n), n = n)
v <- if(q_1 == 2) apply(U[,2,], 2, function(x, n) rep(x, each = n), n = n) else rep(0,N)
e <- apply(e, 3, function(x) t(x))

}

D <- replicate(R, generate.design(n = n, M = M, fixed = fixed))
x <- D[,'x',]
if(!is.matrix(x)) x <- matrix(x)

y <- apply(D, 3, function(x,b) x%*%b, b = beta) + u + x*v + apply(x, 2, function(x,g) cbind(1,x)%*%g, g = gamma)*e
ans <- list(Y = y, X = D, group = id, u = U, e = e)
attr(ans, "call") <- match.call()
attr(ans, "seed") <- seed
return(ans)

}


