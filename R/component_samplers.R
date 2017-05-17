#' Sample the dynamic factors
#'
#' Sample the dynamic factors (latent state space variables) using the simulation smoothing
#' algorithm of Koopman and Durbin (2001), implemented in the KFAS package.

#' @param Y the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param sigma_et \code{T}- or \code{1}-dimensional vector of observation error standard deviation(s)
#' @param Wt \code{K x K} matrix or \code{K x K x T} array of evolution error covariances
#' @param Fmat \code{m x K} matrix of FLCs; only needed for \code{useFastImpute = FALSE}
#' @param YF \code{T x K} matrix of data \code{Y} projected onto FLCs, \code{Y\%*\%Fmat}, which takes the place
#' of \code{Y}; only needed for \code{useFastImpute = TRUE}
#' @param Gt \code{K x K} evolution matrix; if NULL, set as identity (for random walk)
#' @param kfas_model \code{SSModel} object from KFAS package; if NULL, construct model w/in the sampler (might be slower!)
#' @param useFastImpute logical; when TRUE, use imputation/projection scheme for the dynamic factors; otherwise use full state space model for factors (slower)
#' @return The \code{T x K} matrix of dynamic factors, \code{Beta}.
#'
#' @note The sampler has two options: \code{useFastImpute = TRUE}, in which the response is
#' \code{YF = Y\%*\%Fmat} (\code{T x K}) and the observation matrix is the identity (\code{K x K});
#' and \code{useFastImpute = FALSE}, in which the response is \code{Y}  (\code{T x m})
#' and the observation matrix is \code{Fmat} (\code{m x K}). Recall that typically \code{K < < m},
#' so \code{useFastImpute = TRUE} is often much faster.
#'
#' @examples
#' # Read in the yield curve data:
#' data("US_Yields")
#'
#' # Restrict to dates since 2006:
#' Y = Y[which(dates > as.Date("2006-01-01")),];
#'
#' # This is a simple (yet meaningless) example:
#' K = 3
#' Beta = fdlm_factor(Y, sigma_et = 1, Wt = diag(K), Fmat = diag(ncol(Y))[,1:K])
#####################################################################################################
#' @import KFAS
#' @export
fdlm_factor = function(Y, sigma_et, Wt,  Fmat = NULL, YF = NULL, Gt = NULL, kfas_model = NULL, useFastImpute = FALSE){

  # Define these locally
  T = nrow(Y); m = ncol(Y)

  # Update the SSModel object given the new parameters
  if(useFastImpute){ # Fast imputation case: project Y onto F, then sample (requires no missing data in this case!)
    if(is.null(YF)) stop('YF must be specified for fastImpute factor sampler')
    K = ncol(YF) # Store locally
    kfas_model = update.kfas_model(Y.dlm = YF, Zt = array(diag(K), c(K,K,1)), sigma_et = sigma_et, Gt = Gt, Wt = Wt, W0 = NULL, kfas_model = kfas_model)
  } else {# Standard DLM case: Y is the response, Fmat is the observation matrix
    if(is.null(Fmat)) stop('Fmat must be specified for non-fastImpute factor sampler')
    K = ncol(Fmat) # Store locally
    kfas_model = update.kfas_model(Y.dlm = Y, Zt = array(Fmat, c(m,K,1)), sigma_et = sigma_et, Gt = Gt, Wt = Wt, W0 = NULL, kfas_model = kfas_model)
  }

  # Run the sampler
  simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1]
}
#####################################################################################################
#' Evolution error variance sampler
#'
#' Sample the (non-dynamic) evolution error variance matrix.
#'
#' @param resBeta \code{(T-1) x K} matrix of residuals from the evolution equation
#' @param Rinv \code{K x K}; only needed for \code{useDiagonal = FALSE}
#' @param rh0 scalar prior scale parameter; only needed for \code{useDiagonal = FALSE}
#' @param useDiagonal logical; when TRUE, use a diagonal covariance matrix; otherwise, use a full matrix

#' @return The \code{K x K} evolution error covariance matrix, \code{Wt}.
#'
#' @note The sampler has two options: \code{useDiagonal = TRUE}, in which the evolution
#' error variance matrix is diagional with independent Gamma(0.001,0.001) priors on the
#' \code{K} diagonal component precisions; and \code{useDiagonal = FALSE}, in which the
#' \code{K x K} evolution error variance matrix is full with an inverse Wishart prior
#' with prior precision \code{Rinv} and scale parameter \code{rh0}. The default is
#' \code{Rinv = diag(1, K)} and \code{rh0 = K}.
#'
#' @examples
#' # Simulate an example:
#' n = 100; K = 3
#' resBeta = matrix(rnorm(n*K), nr = n, nc = K)
#'
#' # Sample the variance:
#' sample_Wt(resBeta)
#####################################################################################################
#' @export
sample_Wt = function(resBeta, Rinv = diag(1, K), rh0 = K, useDiagonal=FALSE){
  K = ncol(resBeta)
  if(useDiagonal){
    # Assumes independent Gamma(0.001, 0.001) priors for each component
    diag(apply(resBeta, 2, function(x) 1/rgamma(n=1, shape = 0.001 + (length(x)-1)/2, rate = sum(x^2)/2 + 0.001)))
  } else chol2inv(chol(MCMCpack::rwish(rh0 + nrow(resBeta), chol2inv(chol(Rinv*rh0 + crossprod(resBeta))))))
}
#####################################################################################################
#' Factor loading curve smoothing parameter sampler
#'
#' Sample the smoothing parameters for each factor loading curve.
#'
#' @param lambda \code{K}-dimensional vector of smoothing parameters (prior precisions)
#' from previous MCMC iteration
#' @param Psi \code{J x K} matrix of basis coefficients, where \code{J} is the number of
#' basis functions and \code{K} is the number of factors
#' @param Omega \code{J x J} penalty matrix; if NULL, assume it is diag(0, 0, 1,...,1)
#' @param uniformPrior logical; when TRUE, use a uniform prior on prior standard deviations,
#' \code{1/sqrt{lambda[k]}}; otherwise use independent Gamma(0.001, 0.001) prior for each \code{lambda[k]}
#' @param orderLambda logical; when TRUE, enforce the ordering constraint \code{lambda[1] > ... > lambda[K]}
#' for identifiability

#' @return The \code{K}-dimensional vector of samoothing parameters, \code{lambda}.
#####################################################################################################
#' @export
sample_lambda = function(lambda, Psi, Omega = NULL, uniformPrior = TRUE, orderLambdas = TRUE){
  J = nrow(Psi); K = ncol(Psi)

  if(uniformPrior){shape0 = (J + 1)/2} else shape0 = (J - 2)/2 + 0.001; # for Gamma(0.001, 0.001) prior

  for(k in 1:K){
    if(is.null(Omega)){rate0 = crossprod(Psi[-(1:2),k])/2} else rate0 = crossprod(Psi[,k], Omega)%*%Psi[,k]/2
    if(!uniformPrior) rate0 = rate0 + 0.001  # for Gamma(0.001, 0.001) prior

    # Lower and upper bounds, w/ ordering constraints (if specified):
    if(orderLambdas){
      lam.l = 10^-8; lam.u = Inf; if(k != 1) lam.u = lambda[k-1];  # if(k != K) lam.l = lambda[k+1];
      lambda[k] = truncdist::rtrunc(1, 'gamma', a=lam.l, b=lam.u, shape=shape0, rate=rate0) # more stable, possibly faster
    } else lambda[k] = rgamma(1, shape = shape0, rate = rate0)
  }
  lambda
}
