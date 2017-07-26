#####################################################################################################
# Initialize the factors and FLCs using a SVD
# Inputs: Y, tau, K (see fdlm() for details)
# Returns a list of the main parameters in the model: Beta, d, and splineInfo (spline basis matrices)
#####################################################################################################
fdlm_init = function(Y, tau, K){

  # Rescale observation points to [0,1]
  tau01 = (tau - min(tau))/diff(range(tau))

  # And the dimensions:
  T = nrow(Y); m = ncol(Y)

  # Compute basic quantities for the FLC splines:
  splineInfo = getSplineInfo(tau01, m_avg = floor(mean(rowSums(!is.na(Y)))), orthonormal = TRUE)

  # For initialization: impute
  Y0 = matrix(NA, nr = T, nc = m) # Storage for imputed initialization data matrix
  allMissing.t = (rowSums(!is.na(Y))==0)   # Boolean indicator of times at which no points are observed
  # First: for all times at which we observe a curve, impute the full curve (across tau)
  Y0[!allMissing.t,] = t(apply(Y[!allMissing.t,], 1, function(x) splinefun(tau01, x, method='natural')(tau01)))
  # Second: impute any times for which no curve is observed (i.e., impute across time)
  Y0 = apply(Y0, 2, function(x){splinefun(1:T, x, method='natural')(1:T)})

  # Compute SVD of the (completed) data matrix:
  YB0 = Y0%*%splineInfo$Bmat%*%chol2inv(chol(splineInfo$BtB))
  singVal = svd(YB0)

  # If K is unspecified, select based on cpv
  if(is.null(K)){
    # Cumulative sums of the s^2 proportions (More reasonable when Y has been centered)
    cpv = cumsum(singVal$d^2/sum(singVal$d^2))
    K = max(2, which(cpv >= 0.99)[1])
  }

  # Basis coefficients of FLCs:
  Psi0 = as.matrix(singVal$v[,1:K])

  # Initial FLCs:
  F0 = splineInfo$Bmat%*%Psi0

  # Factors:
  Beta0 = (singVal$u%*%diag(singVal$d))[,1:K]

  # Initialize all curves to have positive sums (especially nice for the intercept)
  negSumk = which(colSums(F0) < 0); Psi0[,negSumk] = -Psi0[,negSumk]; Beta0[,negSumk] = -Beta0[,negSumk]

  list(Beta = Beta0, Psi = Psi0, splineInfo = splineInfo)
}
#####################################################################################################
# Update (or initialize) the SSModel object used for sampling the factors
# Inputs:
# Y.dlm: (T x m0) matrix of response (w/ missing values); m0 is usually either K (fastImpute = TRUE) or m (fastImpute = FALSE)
# Zt: (m0 x K0) observation matrix or (m0 x K0 x T) observation array
# Optional inputs:
# Gt: (K0 x K0) evolution matrix or (K0 x K0 x T) array; if NULL, set as identity (for random walk)
# sigma_et: observation error SD(s): can be
# (T x m0) matrix, assuming time-dependent diagonal obs error variance (columns give diagonal elements)
# T-dimensional vector, assuming time-dependent scalar multiple of the identity
# 1-dimensional vector, assuming time-invariant scalar multiple of the identity
# if NULL, set as identity (for dimension purposes)
# Wt: (K0 x K0) matrix or (K0 x K0 x T) array of evolution error covariance; if NULL, set as identity (for dimension purposes)
# W0: (K0 x K0) matrix of initial evolution error covariance; if NULL, set as identity*10^4
# kfas_model: SSModel object from KFAS package; if NULL, construct model (might be slower!)
#####################################################################################################
#' @import KFAS
update_kfas_model = function(Y.dlm, Zt, sigma_et = NULL, Gt = NULL, Wt = NULL, W0 = NULL, kfas_model = NULL){

  if (!requireNamespace("KFAS", quietly = TRUE)) {
    stop("KFAS needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # Compute these locally
  T = nrow(Y.dlm); m0 = ncol(Y.dlm) # m0 is usually either K (fastImpute = TRUE) or m (fastImpute = FALSE)
  K0 = ncol(Zt) # This is the dimension of the state vector

  # Evolution matrix: if NULL (unspecified), set as identity (for random walk model);
  # if (K0 x K0) matrix, transform to (K0 x K0 x T) array
  if(is.null(Gt)) {Gt = array(diag(K0), c(K0,K0,T))} else{if(length(dim(Gt)) == 2) Gt = array(Gt, c(K0, K0, T))}

  # Evolution error covariance matrix: if NULL (unspecified), set as identity (for dimension purposes only);
  # if (K0 x K0) matrix, transform to (K0 x K0 x T) array
  if(is.null(Wt)) {Wt = array(diag(K0), c(K0,K0,T))} else{if(length(dim(Wt)) == 2) Wt = array(Wt, c(K0, K0, T))}

  # Initial variance matrix (not an array in kfas!)
  if(is.null(W0)) W0 = diag(10^4, K0)

  # Observation error variance, which can be time-dependent (if NULL, just set it as the identity)
  Ht = array(diag(m0), c(m0, m0, T))
  if(!is.null(sigma_et)){
    if(length(sigma_et) == 1) sigma_et = rep(sigma_et, T);  #if(length(sigma_et) == 1) Ht = array(diag(sigma_et^2, m0), c(m0, m0, T))
    if(length(sigma_et) == T) sigma_et = tcrossprod(sigma_et, rep(1,m0))
    for(j in 1:m0) Ht[j,j,] = sigma_et[,j]^2 # for(j in 1:m0) Ht[j,j,] = sigma_et^2
  }

  # We can either initialize the SSModel object, kfas_model (when NULL), or update the parameters
  if(is.null(kfas_model)){
    kfas_model = SSModel(Y.dlm ~ -1 + (SSMcustom(Z = Zt, T = Gt, Q = Wt, a1 = matrix(0, nr = K0), P1 = W0, n = T, index = 1:m0)), H = Ht)
  } else {kfas_model$y = Y.dlm; kfas_model$Z = Zt; kfas_model$T = Gt; kfas_model$Q = Wt; kfas_model$P1 = W0; kfas_model$H = Ht}


  # Check for errors
  if(!is.SSModel(kfas_model)) stop("Error: Model has incorrect dimensions")

  # Return the SSModel object
  kfas_model
}
#####################################################################################################
# Imputes the missing observation
# Inputs:
# Yna: (T x m) matrix of observations, including missing values
# mu: (T x m) matrix of conditional expectations (no missing values!)
# sigma_et: T-dimensional vector of observation error SD(s)
# Bmat: (m x J) matrix of spline basis coefficients
#####################################################################################################
fdlm_impute = function(Yna, mu, sigma_et, Bmat){

  # These are the missing indicese:
  na.inds = which(is.na(Yna), arr.ind = TRUE)

  # No missing data: no need to loop, or sample
  if(length(na.inds)==0) return(list(Y = Yna, BtY = tcrossprod(t(Bmat), Y)))

  # Impute by sampling:
  Y = Yna; Y[na.inds] = mu[na.inds] + sigma_et[na.inds[,1]]*rnorm(nrow(na.inds))

  BtY = tcrossprod(t(Bmat), Y) #BtY = t(Y%*%Bmat)

  list(Y = Y, BtY = BtY)
}
#####################################################################################################
# getSplineInfo() initializes (and transforms) the spline basis
# Inputs:
# tau01: all observation points, scaled to [0,1]
# m_avg (=NULL): average number of obs points; if NULL, set m_avg = length(tau01)
# orthonormal: orthonormalize the basis (TRUE/FALSE)
# Notes:
# Uses quantile-based placement of knots for a cubic spline basis
# Enfoces a penalty on the integrated squared second derivative
# Computes the matrix of integrals for the orthonormality constraints
# Transforms the basis, penalty, and matrix of integrals so that:
# d_k is decomposed into linear (first 2 elements) and nonlinear components
# the resulting prior for d_k is diagonal, which helps with posterior mixing, and proper
# Follows Wand and Ormerod (2008)
# Note: for orthonormalized, this is no longer true
#####################################################################################################
getSplineInfo = function(tau01, m_avg = NULL, orthonormal = TRUE){

  # This is the linear component, which appears in both cases
  X<-cbind(1, tau01)

  m = length(tau01);

  # Average number of observation points
  if(is.null(m_avg)) m_avg = m

  # Low-rank thin plate spline

  # Number of knots: if m > 25, use fewer
  if(m > 25){
    num.knots = max(20, min(ceiling(m_avg/4), 150))
  } else num.knots = max(3, ceiling(m_avg/2))

  knots<-quantile(unique(tau01), seq(0,1,length=(num.knots+2))[-c(1,(num.knots+2))])

  # SVD-type reparam (see Ciprian's paper)
  Z_K<-(abs(outer(tau01,knots,"-")))^3; OMEGA_all<-(abs(outer(knots,knots,"-")))^3
  svd.OMEGA_all<-svd(OMEGA_all)
  sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*%(t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d)))

  # The nonlinear component:
  Z<-t(solve(sqrt.OMEGA_all,t(Z_K)))

  # Now combine the linear and nonlinear pieces:
  Bmat = cbind(X, Z);

  # The penalty matrix:
  Omega = diag(c(rep(10^-8, 2), rep(1, (ncol(Bmat) - 2))))

  if(orthonormal){
    # QR decomposition:
    qrd = qr(Bmat, complete = TRUE);  R.t = t(qr.R(qrd));
    # Update hte basis and the penalty matrix:
    Bmat = qr.Q(qrd); Omega = forwardsolve(R.t, t(forwardsolve(R.t, Omega, upper.tri = FALSE)), upper.tri = FALSE)

    BtB = diag(1, ncol(Bmat))
  } else BtB = crossprod(Bmat)

  # Return the matrix, the penalty, and the cross product (of the basis)
  list(Bmat = Bmat, Omega = Omega, BtB = BtB)
}
#####################################################################################################
#' Compute Simultaneous Credible Bands
#'
#' Compute (1-alpha)% credible BANDS for a function based on MCMC samples using Crainiceanu et al. (2007)
#'
#' @param sampFuns \code{Nsims x m} matrix of \code{Nsims} MCMC samples and \code{m} points along the curve
#' @param alpha confidence level
#'
#' @return \code{m x 2} matrix of credible bands; the first column is the lower band, the second is the upper band
#'
#' @note The input needs not be curves: the simultaneous credible "bands" may be computed
#' for vectors. The resulting credible intervals will provide joint coverage at the (1-alpha)%
#' level across all components of the vector.
#'
#' @export
credBands = function(sampFuns, alpha = .05){

  N = nrow(sampFuns); m = ncol(sampFuns)

  # Compute pointwise mean and SD of f(x):
  Efx = colMeans(sampFuns); SDfx = apply(sampFuns, 2, sd)

  # Compute standardized absolute deviation:
  Standfx = abs(sampFuns - tcrossprod(rep(1, N), Efx))/tcrossprod(rep(1, N), SDfx)

  # And the maximum:
  Maxfx = apply(Standfx, 1, max)

  # Compute the (1-alpha) sample quantile:
  Malpha = quantile(Maxfx, 1-alpha)

  # Finally, store the bands in a (m x 2) matrix of (lower, upper)
  cbind(Efx - Malpha*SDfx, Efx + Malpha*SDfx)
}
#####################################################################################################
#' Estimate the remaining time in the MCMC based on previous samples
#' @param nsi Current iteration
#' @param timer0 Initial timer value, returned from \code{proc.time()[3]}
#' @param nsims Total number of simulations
#' @param nrep Print the estimated time remaining every \code{nrep} iterations
#' @return Table of summary statistics using the function \code{summary}
computeTimeRemaining = function(nsi, timer0, nsims, nrep=100){

  # Only print occasionally:
  if(nsi%%nrep == 0 || nsi==20) {
    # Current time:
    timer = proc.time()[3]

    # Simulations per second:
    simsPerSec = nsi/(timer - timer0)

    # Seconds remaining, based on extrapolation:
    secRemaining = (nsims - nsi -1)/simsPerSec

    # Print the results:
    if(secRemaining > 3600) {
      print(paste(round(secRemaining/3600, 1), "hours remaining"))
    } else {
      if(secRemaining > 60) {
        print(paste(round(secRemaining/60, 2), "minutes remaining"))
      } else print(paste(round(secRemaining, 2), "seconds remaining"))
    }
  }
}
#----------------------------------------------------------------------------
#' Summarize of effective sample size
#'
#' Compute the summary statistics for the effective sample size (ESS) across
#' posterior samples for possibly many variables
#'
#' @param postX An array of arbitrary dimension \code{(nsims x ... x ...)}, where \code{nsims} is the number of posterior samples
#' @return Table of summary statistics using the function \code{summary()}.
#'
#' @examples
#' # ESS for iid simulations:
#' rand_iid = rnorm(n = 10^4)
#' getEffSize(rand_iid)
#'
#' # ESS for several AR(1) simulations with coefficients 0.1, 0.2,...,0.9:
#' rand_ar1 = sapply(seq(0.1, 0.9, by = 0.1), function(x) arima.sim(n = 10^4, list(ar = x)))
#' getEffSize(rand_ar1)
#'
#' @import coda
#' @export
getEffSize = function(postX) {
  if(is.null(dim(postX))) return(effectiveSize(postX))
  summary(effectiveSize(as.mcmc(array(postX, c(dim(postX)[1], prod(dim(postX)[-1]))))))
}
#----------------------------------------------------------------------------
#' Compute the ergodic (running) mean.
#' @param x vector for which to compute the running mean
#' @return A vector \code{y} with each element defined by \code{y[i] = mean(x[1:i])}
#' @examples
#' # Compare:
#' ergMean(1:10)
#' mean(1:10)
#'
#'# Running mean for iid N(5, 1) samples:
#' x = rnorm(n = 10^4, mean = 5, sd = 1)
#' plot(ergMean(x))
#' abline(h=5)
#' @export
ergMean = function(x) {cumsum(x)/(1:length(x))}
#----------------------------------------------------------------------------
#' Plot the factors
#'
#' Plot posterior mean of the factors together with the simultaneous and pointwise
#' 95\% credible bands.
#'
#' @param post_beta the \code{Nsims x T x K} array of \code{Nsims} draws from the posterior
#' distribution of the \code{T x K} matrix of factors, \code{beta}
#' @param dates \code{T x 1} vector of dates or labels corresponding to the time points
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics abline lines  par plot polygon
#' @import viridis
#' @import coda
#' @export
plot_factors = function(post_beta, dates = NULL){
  K = dim(post_beta)[3] # Number of factors
  if(is.null(dates)) dates = seq(0, 1, length.out = dim(post_beta)[2])

  dev.new(); par(mai = c(.8,.9,.4,.4), bg = 'gray90');
  plot(dates, post_beta[1,,1], ylim = range(post_beta), xlab = 'Dates', ylab = '', main = paste('Dynamic Factors', sep=''), type='n', cex.lab = 2, cex.axis=2,cex.main=2)
  abline(h = 0, lty=3, lwd=2);
  for(k in K:1){
    cb = credBands(as.mcmc(post_beta[,,k])); ci = HPDinterval(as.mcmc(post_beta[,,k]));
    polygon(c(dates, rev(dates)), c(cb[,2], rev(cb[,1])), col='grey50', border=NA);
    polygon(c(dates, rev(dates)), c(ci[,2], rev(ci[,1])), col='grey', border=NA);
    lines(dates,colMeans(post_beta[,,k]), lwd=8, col=viridis(K, end = .8)[k])
  }
}
#----------------------------------------------------------------------------
#' Plot the factor loading curves
#'
#' Plot posterior mean of the factor loading curves together with the simultaneous
#' and pointwise 95\% credible bands.
#'
#' @param post_fk the \code{Nsims x m x K} array of \code{Nsims} draws from the posterior
#' distribution of the \code{m x K} matrix of FLCs, \code{fk}
#' @param tau \code{m x 1} vector of observation points
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics abline lines  par plot polygon
#' @import viridis
#' @import coda
#' @export
plot_flc = function(post_fk, tau = NULL){
  K = dim(post_fk)[3] # Number of factors
  if(is.null(tau)) tau = seq(0, 1, length.out = dim(post_fk)[2])

  dev.new(); par(mai = c(.9,.9,.4,.4), bg = 'gray90');
  plot(tau, post_fk[1,,1], ylim = range(post_fk), xlab = expression(tau), ylab = '', main = 'Factor Loading Curves', type='n', cex.lab = 2, cex.axis=2,cex.main=2)
  abline(h = 0, lty=3, lwd=2);
  for(k in K:1){
    cb = credBands(as.mcmc(post_fk[,,k])); ci = HPDinterval(as.mcmc(post_fk[,,k]));
    polygon(c(tau, rev(tau)), c(cb[,2], rev(cb[,1])), col='grey50', border=NA);
    polygon(c(tau, rev(tau)), c(ci[,2], rev(ci[,1])), col='grey', border=NA);
    lines(tau,colMeans(post_fk[,,k]), lwd=8, col=viridis(K, end = .8)[k])
  }
}
# Just add these for general use:
#' @importFrom stats quantile rgamma rnorm sd splinefun var
NULL
