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
#' @param W0 \code{K x K} matrix of initial evolution error covariances; if NULL, set to diag(10^-4, K)
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
fdlm_factor = function(Y, sigma_et, Wt,  Fmat = NULL, YF = NULL, Gt = NULL, W0 = NULL, kfas_model = NULL, useFastImpute = FALSE){

  # Define these locally
  T = nrow(Y); m = ncol(Y)

  # Update the SSModel object given the new parameters
  if(useFastImpute){ # Fast imputation case: project Y onto F, then sample (requires no missing data in this case!)
    if(is.null(YF)) stop('YF must be specified for fastImpute factor sampler')
    K = ncol(YF) # Store locally
    kfas_model = update_kfas_model(Y.dlm = YF, Zt = array(diag(K), c(K,K,1)), sigma_et = sigma_et, Gt = Gt, Wt = Wt, W0 = NULL, kfas_model = kfas_model)
  } else {# Standard DLM case: Y is the response, Fmat is the observation matrix
    if(is.null(Fmat)) stop('Fmat must be specified for non-fastImpute factor sampler')
    K = ncol(Fmat) # Store locally
    kfas_model = update_kfas_model(Y.dlm = Y, Zt = array(Fmat, c(m,K,1)), sigma_et = sigma_et, Gt = Gt, Wt = Wt, W0 = NULL, kfas_model = kfas_model)
  }

  # Run the sampler
  as.matrix(simulateSSM(kfas_model, "states", nsim = 1, antithetics=FALSE, filtered=FALSE)[,,1])
}
#' Factor Loading Curve Sampling Algorithm
#'
#' Sample the factor loading curve basis coefficients subject to an orthonormality constraint.
#' Additional linear constraints may be included as well.
#'
#' @param BtY \code{J x T} matrix \code{B.t()*Y} for basis matrix B
#' @param Beta \code{T x K} matrix of factors
#' @param Psi \code{J x K} matrix of previous factor loading curve coefficients
#' @param BtB \code{J x J} matrix of \code{B.t()*B}
#' @param Omega \code{J x J} prior precision/penalty matrix
#' @param BtCon (optional) \code{J x Jc} matrix of additional constraints, pre-multiplied by B.t()
#' @param lambda \code{K}-dimensional vector of prior precisions
#' @param sigmat2 \code{T}-dimensional vector of time-dependent observation error variances
#' @return Psi \code{J x K} matrix of factor loading curve coefficients
#'
#' @note This is a wrapper for Rcpp functions for the special cases of
#' \code{K = 1} and whether or not additional (linear) constraints are included,
#' i.e., whether or not \code{BtCon} is non-\code{NULL}.
#' @export
fdlm_flc = function(BtY, Beta, Psi, BtB, Omega, BtCon = NULL, lambda, sigmat2){

  # Obtain the dimensions, in order of appearance:
  J = nrow(BtY); T = ncol(BtY); K = ncol(Beta);

  # Allow for scalar variance input
  if(length(sigmat2) == 1) sigmat2 = rep(sigmat2, T)

  # Check dimensions:
  if( (nrow(Beta) != T) ||
      (nrow(Psi) != J) || (ncol(Psi) != K) ||
      (nrow(BtB) != J) || (ncol(BtB) != J) ||
      (nrow(Omega) != J) || (ncol(Omega) != J) ||
      (length(lambda) != K) ||
      (length(sigmat2) != T)
  ) stop("Mismatched dimensions in FLC sampler")

  # No additional constraints (besides orthonormality of FLCs themselves)
  if(is.null(BtCon)){

    if(K == 1){
      # Special case: (FLC) orthogonality not necessary
      Psi = sampleFLC_1(BtY = BtY, Beta = Beta, Psi = Psi, BtB = BtB, Omega = Omega, lambda = lambda, sigmat2 = sigmat2)
    } else {
      Psi = sampleFLC(BtY = BtY, Beta = Beta, Psi = Psi, BtB = BtB, Omega = Omega, lambda = lambda, sigmat2 = sigmat2)
    }


  } else {
    # Additional constraints: orthogonal to BtCon


    # Special case: (FLC) orthogonality not necessary
    if(K == 1){
      Psi = sampleFLC_cons_1(BtY = BtY, Beta = Beta, Psi = Psi, BtB = BtB, Omega = Omega, BtCon = BtCon, lambda = lambda, sigmat2 = sigmat2)
    } else {
      Psi = sampleFLC_cons(BtY = BtY, Beta = Beta, Psi = Psi, BtB = BtB, Omega = Omega, BtCon = BtCon, lambda = lambda, sigmat2 = sigmat2)
    }

  }
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
  K = ncol(as.matrix(resBeta))
  if(useDiagonal){
    # Assumes independent Gamma(0.001, 0.001) priors for each component
    if(K > 1){
      diag(apply(resBeta, 2, function(x) 1/rgamma(n=1, shape = 0.001 + (length(x)-1)/2, rate = sum(x^2)/2 + 0.001)))
    } else  diag(1/rgamma(n=1, shape = 0.001 + (length(resBeta)-1)/2, rate = sum(resBeta^2)/2 + 0.001), 1)
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
#' @param d dimension of \code{tau}; default is 1
#' @param uniformPrior logical; when TRUE, use a uniform prior on prior standard deviations,
#' \code{1/sqrt{lambda[k]}}; otherwise use independent Gamma(0.001, 0.001) prior for each \code{lambda[k]}
#' @param orderLambdas logical; when TRUE, enforce the ordering constraint \code{lambda[1] > ... > lambda[K]}
#' for identifiability

#' @return The \code{K}-dimensional vector of samoothing parameters, \code{lambda}.
#####################################################################################################
#' @export
sample_lambda = function(lambda, Psi, Omega = NULL, d = 1, uniformPrior = TRUE, orderLambdas = TRUE){
  J = nrow(Psi); K = ncol(Psi)

  if(uniformPrior){shape0 = (J - d + 1 + 1)/2} else shape0 = (J - d - 1)/2 + 0.001; # for Gamma(0.001, 0.001) prior
  #if(uniformPrior){shape0 = (J + 1)/2} else shape0 = (J - 2)/2 + 0.001; # for Gamma(0.001, 0.001) prior

  for(k in 1:K){
    if(is.null(Omega)){rate0 = crossprod(Psi[-(1:(d+1)),k])/2} else rate0 = crossprod(Psi[,k], Omega)%*%Psi[,k]/2
    #if(is.null(Omega)){rate0 = crossprod(Psi[-(1:2),k])/2} else rate0 = crossprod(Psi[,k], Omega)%*%Psi[,k]/2

    if(!uniformPrior) rate0 = rate0 + 0.001  # for Gamma(0.001, 0.001) prior

    # Lower and upper bounds, w/ ordering constraints (if specified):
    if(orderLambdas){
      lam.l = 10^-8; lam.u = Inf; if(k != 1) lam.u = lambda[k-1];  # if(k != K) lam.l = lambda[k+1];
      lambda[k] = truncdist::rtrunc(1, 'gamma', a=lam.l, b=lam.u, shape=shape0, rate=rate0) # more stable, possibly faster
    } else lambda[k] = rgamma(1, shape = shape0, rate = rate0)
  }
  lambda
}
#' Sample the autoregressive coefficients in an AR(1) Model
#'
#' Compue one draw of the autoregressive coefficients \code{phi} in an AR(1) model.
#' The sampler also applies to a multivariate case with independent components.
#'
#' Sample the AR(1) coefficients \code{phi_j} using the model
#'
#' \code{y_tj = mu_j + phi_j(y_{t-1,j} - mu_j) + e_tj},
#'
#' with \code{e_tj ~ N(0, sigma[j]^2)}
#'
#' @param yt the \code{T x p} matrix of centered multivariate time series
#' (i.e., the time series minus the unconditional means, \code{mu})
#' @param phi_j the \code{p x 1} vector of previous AR(1) coefficients
#' @param sigma_tj the \code{(T-1) x p} matrix or \code{p x 1} vector of error standard deviations
#' @param prior_phi the parameters of the prior for the AR(1) coefficients \code{phi_j};
#' either \code{NULL} for uniform on [-0.99,0.99] or a 2-dimensional vector of (shape1, shape2) for a Beta prior
#' on \code{[(phi_j + 1)/2]}
#'
#' @return \code{p x 1} vector of sampled AR(1) coefficient(s)
#'
#' @note For the standard AR(1) case, \code{p = 1}. However, the function applies more
#' generally for sampling \code{p > 1} independent AR(1) processes (jointly).
#'
#' @import truncdist
#' @export
sampleARphi = function(yt, phi_j, sigma_tj, prior_phi = NULL){

  # Just in case:
  yt = as.matrix(yt)

  # Store dimensions locally:
  T = nrow(yt); p = ncol(yt)

  if(length(sigma_tj) == p) sigma_tj = matrix(rep(sigma_tj, each = T-1), nrow = T-1)

  # Loop over the j=1:p
  for(j in 1:p){

    # Compute "regression" terms for dhs_phi_j:
    y_ar = yt[-1,j]/sigma_tj[,j] # Standardized "response"
    x_ar = yt[-T,j]/sigma_tj[,j] # Standardized "predictor"

    # Using Beta distribution:
    if(!is.null(prior_phi)){

      # Check to make sure the prior params make sense
      if(length(prior_phi) != 2) stop('prior_phi must be a numeric vector of length 2')

      phi01 = (phi_j[j] + 1)/2 # ~ Beta(prior_phi[1], prior_phi[2])

      # Slice sampler when using Beta prior:
      phi01 = uni.slice(phi01, g = function(x){
        -0.5*sum((y_ar - (2*x - 1)*x_ar)^2) +
          dbeta(x, shape1 = prior_phi[1], shape2 = prior_phi[2], log = TRUE)
      }, lower = 0.005, upper = 0.995)[1]

      phi_j[j] = 2*phi01 - 1

    } else {
      # For phi_j ~ Unif(-0.99, 0.99), the posterior is truncated normal
      phi_j[j] = rtrunc(n = 1, spec = 'norm',
                        a = -0.99, b = 0.99,
                        mean = sum(y_ar*x_ar)/sum(x_ar^2),
                        sd = 1/sqrt(sum(x_ar^2)))
    }
  }
  return(phi_j)
}
#' Sample the unconditional mean in an AR(1) Model
#'
#' Compue one draw of the unconditional mean \code{mu} in an AR(1) model assuming a
#' Gaussian prior (with mean zero).
#'
#' Sample the unconditional mean \code{mu} using the model
#'
#' \code{y_tj = mu_j + phi_j(y_{t-1,j} - mu_j) + e_tj},
#'
#' with \code{e_tj ~ N(0, sigma[j]^2)} and prior \code{mu ~ N(0, 1/priorPrec[j])}
#'
#' @param yt the \code{T x p} matrix of multivariate time series
#' @param phi_j the \code{p x 1} vector of AR(1) coefficients
#' @param sigma_tj the \code{(T-1) x p} matrix or \code{p x 1} vector of error standard deviations
#' @param priorPrec the \code{p x 1} vector of prior precisions;
#' if \code{NULL}, use \code{rep(10^-6, p)}
#'
#' @return The \code{p x 1} matrix of unconditional means.
#' @export
sampleARmu = function(yt, phi_j, sigma_tj, priorPrec = NULL){

  # Just in case:
  yt = as.matrix(yt)

  # Store dimensions locally:
  T = nrow(yt); p = ncol(yt)

  if(length(sigma_tj) == p) sigma_tj = matrix(rep(sigma_tj, each = T-1), nrow = T-1)

  # Prior precision:
  if(is.null(priorPrec)) priorPrec = rep(10^-6, p)

  # Now, form the "y" and "x" terms in the (auto)regression
  y_mu = (yt[-1,] - matrix(rep(phi_j, each = T-1), nrow = T-1)*yt[-T,])/sigma_tj
  x_mu = matrix(rep(1 - phi_j, each = T-1), nrow = T-1)/sigma_tj

  # Posterior SD and mean:
  postSD = 1/sqrt(colSums(x_mu^2) + priorPrec)
  postMean = (colSums(x_mu*y_mu))*postSD^2
  mu = rnorm(n = p, mean = postMean, sd = postSD)

  return(mu)
}
#' Sample the unconditional mean in a VAR
#'
#' Compue one draw of the unconditional mean \code{mu} in a VAR assuming a
#' Gaussian prior (with mean zero).
#'
#' Sample the unconditional mean \code{mu} using the model
#'
#' \code{y_t = mu + G(y_{t-1} - mu) + e_t},
#'
#' with \code{e_t ~ N(0, Sigma)} and prior \code{mu ~ N(0, solve(priorPrec))}
#'
#' @param yt the \code{T x p} matrix of multivariate time series
#' @param G the \code{p x p} VAR coefficient matrix
#' @param Sigma the \code{p x p} error variance matrix
#' @param priorPrec the \code{p x p} prior precision matrix;
#' if \code{NULL}, use \code{diag(10^-6, p)}
#'
#' @return The \code{p x 1} matrix of unconditional means.
#' @export
sampleVARmu = function(yt, G, Sigma, priorPrec = NULL){

  # Store dimensions locally:
  T = nrow(yt); p = ncol(yt)

  # Prior precision:
  if(is.null(priorPrec)) priorPrec = diag(10^-6, p)

  # Recurring terms:
  I_p = diag(p)
  I_Gt_SigInv = crossprod(I_p - G, chol2inv(chol(Sigma)))

  # Cholesky of Quadratic term:
  chQ_mu = chol(priorPrec + (T-1)*I_Gt_SigInv%*%(I_p - G))

  # Linear term:
  l_mu = as.matrix(I_Gt_SigInv%*%colSums((yt[-1,] - t(tcrossprod(G, as.matrix(yt[-T,]))))))

  # And return:
  as.matrix(backsolve(chQ_mu,forwardsolve(t(chQ_mu), l_mu) + rnorm(p)))
}
#' Sample the VAR coefficient matrix
#'
#' Compue one draw of the VAR coefficient matrix \code{G} assuming a
#' Gaussian prior (with mean zero). The sampler also assumes the VAR
#' is lag-1.
#'
#' @param ytc the \code{T x p} matrix of (centered) multivariate time series
#' @param Sigma the \code{p x p} error variance matrix
#' @param priorPrec the \code{p x p} prior precision matrix;
#' if \code{NULL}, use \code{diag(10^-6, p)}
#' @param stationary logical; if \code{TRUE}, resample until stationary
#'
#' @return The \code{p x p} VAR coefficient matrix.
#' @export
sampleVAR = function(ytc, Sigma, priorPrec = NULL, stationary = TRUE){

  # Store dimensions locally:
  T = nrow(ytc); p = ncol(ytc)

  # Prior precision:
  if(is.null(priorPrec)) priorPrec = diag(10^-6, p^2)

  Sigma_inv = chol2inv(chol(Sigma))

  # Choleksy of quadratic term:
  chQg = chol(priorPrec + kronecker(crossprod(ytc[-T,]), Sigma_inv))
  # Linear term:
  lg = matrix(Sigma_inv%*%crossprod(ytc[-1,], ytc[-T,]))

  # And sample the VAR coefficient matrix:
  G = matrix(backsolve(chQg,forwardsolve(t(chQg), lg) + rnorm(p^2)), nrow = p, byrow = FALSE)

  # To enforce stationarity, resample until we obtain a stationary matrix
  # (or until the we exceed 1000 tries...)
  if(stationary){
    counter = 0 # So we don't do this forever...
    while(!all(abs(eigen(G, only.values=TRUE)$values) < 1) && counter < 1000){
      # Nonstationary, so resample
      G = matrix(backsolve(chQg,forwardsolve(t(chQg), lg) + rnorm(p^2)), nrow = p, byrow = FALSE)
      # And update counter
      counter = counter + 1
    }
  }

  return(G)
}
#--------------------------------------------------------------
#' Multiplicative Gamma Process (MGP) Sampler
#'
#' Sample the global parameters, delta.h, with tau.h = cumprod(delta.h), from
#' the MGP prior for factor models
#'
#' @param theta.jh the \code{p x K} matrix with entries theta.jh ~ N(0, tau.h^-1), j=1:p, h=1:K
#' @param delta.h the \code{K}-dimensional vector of previous delta.h values, tau.h[h] = prod(delta.h[1:h])
#' @param a1 the prior parameter for factor 1: delta.h[1] ~ Gamma(a1, 1)
#' @param a2 the prior parameter for factors 2:K: delta.h[h] ~ Gamma(a2, 1) for h = 2:K
#' @return \code{delta.h}, the \code{K}-dimensional vector of multplicative components,
#' where tau.h[h] = prod(delta.h[1:h])
#'
#' @note The default \code{a1 = 2} and \code{a2 = 3} appears to offer the best performance
#' in Durante (2017).
#' @export
sampleMGP = function(theta.jh, delta.h, a1 = 2, a2 = 3){

  # Just in case:
  theta.jh = as.matrix(theta.jh)

  # Store the dimensions locally
  p = nrow(theta.jh); K = ncol(theta.jh)

  # Sum over the (squared) replicates:
  sum.theta.l = colSums(theta.jh^2)

  # h = 1 case is separate:
  tau.not.1 = cumprod(delta.h)/delta.h[1]
  delta.h[1] = rgamma(n = 1, shape = a1 + p*K/2,
                      rate = 1 + 1/2*sum(tau.not.1*sum.theta.l))
  # h > 1:
  if(K > 1){for(h in 2:K){
    tau.not.h = cumprod(delta.h)/delta.h[h]
    delta.h[h] = rgamma(n = 1, shape = a2 + p/2*(K - h + 1),
                        rate = 1 + 1/2*sum(tau.not.h[h:K]*sum.theta.l[h:K]))
  }}
  delta.h #list(tau.h = cumprod(delta.h), delta.h = delta.h)
}
