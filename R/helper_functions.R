#####################################################################################################
# Initialize the factors and FLCs using a SVD
# Inputs: Y, tau, K (see fdlm() for details)
# Returns a list of the main parameters in the model:
  # Beta, d, splineInfo (spline basis matrices), and the imputed data matrix Y0
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
  Beta0 = as.matrix((singVal$u%*%diag(singVal$d))[,1:K])

  # Initialize all curves to have positive sums (especially nice for the intercept)
  negSumk = which(colSums(F0) < 0); Psi0[,negSumk] = -Psi0[,negSumk]; Beta0[,negSumk] = -Beta0[,negSumk]

  list(Beta = Beta0, Psi = Psi0, splineInfo = splineInfo, Y0 = Y0)
}
#####################################################################################################
# Initialize the factors and FLCs using a SVD
  # Inputs: Y, tau, K (see fdlm() for details)
# Returns a list of the main parameters in the model:
# Beta, d, splineInfo (spline basis matrices), and the imputed data matrix Y0
#####################################################################################################
fdlm_init_d = function(Y, tau, K){

  # Convert to matrix, if necessary:
  tau = as.matrix(tau)

  # Rescale observation points to [0,1]
  tau01 = apply(tau, 2, function(x) (x - min(x))/(max(x) - min(x)))

  # And the dimensions:
  T = nrow(Y); m = ncol(Y); d = ncol(tau)

  # Compute basic quantities for the FLC splines:
  if(d == 1){
    splineInfo = getSplineInfo(as.numeric(tau01), m_avg = floor(mean(rowSums(!is.na(Y)))), orthonormal = TRUE)
  } else splineInfo = getSplineInfo_d(tau01, num.knots = NULL, orthonormalize = TRUE)

  # For initialization: impute (this is a bit crude)
  Y0 = matrix(NA, nr = T, nc = m) # Storage for imputed initialization data matrix
  allMissing.t = (rowSums(!is.na(Y))==0)   # Boolean indicator of times at which no points are observed
  # First: for all times at which we observe a curve, impute the full curve (across tau)
  #Y0[!allMissing.t,] = t(apply(Y[!allMissing.t,], 1, function(x) splinefun(tau01, x, method='natural')(tau01)))
  Y0[!allMissing.t,] = t(apply(Y[!allMissing.t,], 1, function(x) splinefun(seq(0, 1, length.out = m), x, method='natural')(seq(0, 1, length.out = m))))
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
  Beta0 = as.matrix((singVal$u%*%diag(singVal$d))[,1:K])

  # Initialize all curves to have positive sums (especially nice for the intercept)
  negSumk = which(colSums(F0) < 0); Psi0[,negSumk] = -Psi0[,negSumk]; Beta0[,negSumk] = -Beta0[,negSumk]

  list(Beta = Beta0, Psi = Psi0, splineInfo = splineInfo, Y0 = Y0)
}
#####################################################################################################
# Initialize the parameters in the Dynamic Nelson-Siegel Model
#
dns_init = function(Y, tau, orthogonalize = TRUE){

  # Rescale observation points to [0,1]
  tau01 = (tau - min(tau))/diff(range(tau))

  # And the dimensions:
  T = nrow(Y); m = ncol(Y)

  # Initial value of lambda_ns: Diebold and Li value, scaled for yearly maturities (D+L use months)
  lambda_ns = .0609*12 # 0.0609 is D+L's choice to maximize the curvature factor at 30 months
  F_ns = f_ns(tau, lambda_ns)

  # Orthogonalize:
  if(orthogonalize) F_ns = qr.Q(qr(F_ns))

  # Check for signs, then flip if necessary:
  if(sum(F_ns[,1]) < 0) F_ns[,1] = -F_ns[,1]

  # For initialization: impute
  Y0 = matrix(NA, nr = T, nc = m) # Storage for imputed initialization data matrix
  allMissing.t = (rowSums(!is.na(Y))==0)   # Boolean indicator of times at which no points are observed
  # First: for all times at which we observe a curve, impute the full curve (across tau)
  Y0[!allMissing.t,] = t(apply(Y[!allMissing.t,], 1, function(x) splinefun(tau01, x, method='natural')(tau01)))
  # Second: impute any times for which no curve is observed (i.e., impute across time)
  Y0 = apply(Y0, 2, function(x){splinefun(1:T, x, method='natural')(1:T)})

  # Initial factors:
  if(orthogonalize){
    Beta_ns = Y0%*%F_ns
  } else Beta_ns = Y0%*%F_ns%*%chol2inv(chol(crossprod(F_ns)))

  list(Beta_ns = Beta_ns, F_ns = F_ns, lambda_ns = lambda_ns)
}
#####################################################################################################
#' Initialize the parametric terms
#'
#' Compute initial values for the factors and nonlinear parameter (if necessary)
#' of the parametric component.
#'
#' @param Y the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau vector of observation points (\code{m}-dimensional)
#' @param f_p a function to compute the parametric component, which must return a \code{m x K_p} matrix
#' for \code{K_p} the number of parametric curves; may include a (scalar) nonlinear parameter argument
#' @param orthogonalize logical; when TRUE, orthogonalize the loading curve matrix
#'
#' @return a list containing
#' \itemize{
#' \item \code{Beta_p} the \code{T x K_p} matrix of parametric factors
#' \item \code{F_p} the \code{m x K_p} matrix of parametric loading curves
#' \item \code{lambda_p} the scalar nonlinear parameter
#' }
#'
#' @details Compute initial values via the following algorithm:
#' \enumerate{
#' \item Impute missing values in \code{Y}
#' \item Initialize \code{lambda_p = 1} and compute \code{F_p}; orthogonalize if specified
#' \item Estimate \code{Beta_p} via least squares
#' \item Estimate an initial standard deviation \code{sigma_0} via conditional MLE
#' \item Estimate \code{lambda_p} via conditional MLE and recompute \code{F_p}; orthogonalize if specified
#' }
#' @import KFAS
par_init = function(Y, tau, f_p, orthogonalize = TRUE){

  # Rescale observation points to [0,1]
  tau01 = (tau - min(tau))/diff(range(tau))

  # And the dimensions:
  T = nrow(Y); m = ncol(Y)

  # Initial value of lambda_p
  lambda_p = 1
  F_p = f_p(tau, lambda_p)

  # Orthogonalize:
  if(orthogonalize) F_p = qr.Q(qr(F_p))

  # For initialization: impute
  Y0 = matrix(NA, nr = T, nc = m) # Storage for imputed initialization data matrix
  allMissing.t = (rowSums(!is.na(Y))==0)   # Boolean indicator of times at which no points are observed
  # First: for all times at which we observe a curve, impute the full curve (across tau)
  Y0[!allMissing.t,] = t(apply(Y[!allMissing.t,], 1, function(x) splinefun(tau01, x, method='natural')(tau01)))
  # Second: impute any times for which no curve is observed (i.e., impute across time)
  Y0 = apply(Y0, 2, function(x){splinefun(1:T, x, method='natural')(1:T)})

  # Initial factors:
  if(orthogonalize){
    Beta_p = Y0%*%F_p
  } else Beta_p = Y0%*%F_p%*%chol2inv(chol(crossprod(F_p)))

  sigma_0 = sd(Y - tcrossprod(Beta_p, F_p), na.rm=TRUE)

  # Compute an optimum:
  lambda_p =  optim(lambda_p, fn = function(x){
    F_p_x = f_p(tau, x); if(orthogonalize) F_p_x = qr.Q(qr(F_p_x))
    sum(0.5*rowSums((tcrossprod(Beta_p, F_p_x) - Y)^2, na.rm = TRUE)/sigma_0^2)
  }, method = "L-BFGS-B", lower = 10^-4, upper = 10^4)$par

  # Update F_p
  F_p = f_p(tau, lambda_p)

  # Orthogonalize:
  if(orthogonalize) F_p = qr.Q(qr(F_p))

  # Check for signs, then flip if necessary:
  if(sum(F_p[,1]) < 0) {F_p[,1] = -F_p[,1]; Beta_p[,1] = -Beta_p[,1]}


  list(Beta_p = Beta_p, F_p = F_p, lambda_p = lambda_p)
}
#####################################################################################################
#' Initialize the parametric terms using MLEs
#'
#' Compute initial values for the factors and nonlinear parameter (if necessary)
#' of the parametric component using the Kalman Filter.
#'
#' @param Y the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau vector of observation points (\code{m}-dimensional)
#' @param f_p a function to compute the parametric component, which must return a \code{m x K_p} matrix
#' for \code{K_p} the number of parametric curves; may include a (scalar) nonlinear parameter argument
#' @param orthogonalize logical; when TRUE, orthogonalize the loading curve matrix
#'
#' @return a list containing
#' \itemize{
#' \item \code{Beta_p} the \code{T x K_p} matrix of parametric factors
#' \item \code{F_p} the \code{m x K_p} matrix of parametric loading curves
#' \item \code{lambda_p} the scalar nonlinear parameter
#' }
#'
#' @details Compute initial values via the following algorithm:
#' \enumerate{
#' \item Impute missing values in \code{Y}
#' \item Initialize \code{lambda_p = 1} and compute \code{F_p}; orthogonalize if specified
#' \item Estimate \code{Beta_p} via least squares
#' \item Estimate an initial standard deviation \code{sigma_0} via conditional MLE
#' \item Estimate \code{lambda_p} via conditional MLE and recompute \code{F_p}; orthogonalize if specified
#' }
#' @import KFAS
par_init_ssm = function(Y, tau, f_p, orthogonalize = TRUE){

  # Compute the dimensions:
  T = nrow(Y); m = ncol(Y)

  # Initial loading matrix values (and dimensions), mostly for constructing ssmodel
  F_p = f_p(tau, 1); K_p = ncol(F_p)
  if(orthogonalize) F_p = qr.Q(qr(F_p))

  # Initial parameters: nonlinear parameter, obs log-variance, evol diagonal log-variance
  pars0 = c(1, 2*log(sd(Y, na.rm = TRUE)), rep(0, K_p))

  # Initial state space model (KFAS)
  kfas_model = SSModel(Y ~ -1 +
                         (SSMcustom(Z = F_p,
                                    T = diag(K_p),
                                    Q = diag(K_p),
                                    a1 = matrix(0, nr = K_p),
                                    P1 = diag(10^4, K_p),
                                    n = T, index = 1:m)),
                       H = diag(m))

  # Run the optimizer:
  nmOpt = optim(par = pars0, fn = function(params){
    -logLik(update_model(params,kfas_model, tau, f_p, orthogonalize))
  }, method = 'Nelder-Mead')

  # Loading curves:
  lambda_p = nmOpt$par[1]
  F_p = f_p(tau, lambda_p); if(orthogonalize) F_p = qr.Q(qr(F_p))

  # And the factors:
  kfas_model = update_model(nmOpt$par, kfas_model, tau, f_p, orthogonalize)
  Beta_p = KFS(kfas_model, filtering="none", smoothing="state")$alphahat

  # Return the parameters:
  list(Beta_p = Beta_p, F_p = F_p, lambda_p = lambda_p)
}
#####################################################################################################
#' Update the KFAS model
#'
#' Update the KFAS model using the given parameters, assuming a specific form.
#'
#' @param pars vector of parameters
#' @param model object of class SSModel() which should be updated
#' @param tau vector of observation points (\code{m}-dimensional)
#' @param f_p a function to compute the parametric component, which must return a \code{m x K_p} matrix
#' for \code{K_p} the number of parametric curves; may include a (scalar) nonlinear parameter argument
#' @param orthogonalize logical; when TRUE, orthogonalize the loading curve matrix
#' @return The updated \code{model}
#' @note This function is necessary for state space optimization (initialization)
update_model = function(pars, model, tau, f_p, orthogonalize){

  m = length(tau)

  # Update the loadings matrix:
  F_p = f_p(tau, pars[1]); if(orthogonalize) F_p = qr.Q(qr(F_p))
  model$Z[,,1] = F_p

  K_p = ncol(F_p)

  # Update the observation error variance:
  model$H[1:m,1:m,1] = diag(exp(pars[2]), m)

  # Update the evoluation error variance:
  for(k in 1:K_p) model$Q[k,k,1] = exp(pars[2 + k])
  model
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
#' Construct the spline basis and penalty matrices
#'
#' Given input points in \code{d} dimensions, construct a low-rank thin plate spline basis matrix
#' and penalty matrix.
#'
#' @param tau \code{m x d} matrix of coordinates, where \code{m} is the number of observation points and \code{d} is the dimension
#' @param num.knots the number of knots to include
#' @param orthonormalize logical; if TRUE, orthornomalize the basis matrix
#'
#' @note The knot locations are selected using a space-filling algorithm.
#'
#' @import fields
getSplineInfo_d = function(tau, num.knots = NULL, orthonormalize = TRUE){

  # Just in case:
  tau = as.matrix(tau)

  # Number of observation points
  m = nrow(tau)

  # Dimension:
  d = ncol(tau)

  # Order of derivative in penalty:
  m.deriv = 2

  # This is the linear component
  X = cbind(1, tau)

  # Select the number and location of knots
  if(is.null(num.knots)) num.knots = max(20, min(floor(m/4), 150))

  if(num.knots < m){
    #knots0 = cover.design(tau, num.knots); knots = NULL; for(j in 1:d) knots = cbind(knots, knots0[,j])
    knots =  cover.design(tau, num.knots)$design
  } else knots = tau

  # For the penalty matrix, need to compute distances between obs. points and knots
  dist.mat <- matrix(0, num.knots, num.knots); dist.mat[lower.tri(dist.mat)] <- dist(knots); dist.mat <- dist.mat + t(dist.mat)
  if(d%%2 == 0){
    # Even dim:
    Omega = dist.mat^(2*m.deriv - d)*log(dist.mat)
  } else {
    # Odd dim:
    Omega = dist.mat^(2*m.deriv - d)
  }
  # For numerical stability:
  diag(Omega) = 0

  # Compute the "random effects" matrix
  Zk = matrix(0, nrow=m, ncol=num.knots)
  for (k in 1:num.knots){
    di = sqrt(rowSums((tau - matrix(rep(knots[k,], each = m), nrow=m))^2)) # di = 0; for(j in 1:d) di = di + (tau[,j] - knots[k,j])^2; di = sqrt(di)
    if(d%%2 == 0){# Even dim:
      Zk[,k] = di^(2*m.deriv - d)*log(di)
    } else { # Odd dim:
      Zk[,k] = di^(2*m.deriv - d)
    }
  }
  Zk[is.nan(Zk)] = 0

  # Natural constraints, if necessary:
  if(num.knots > m - 2){Q2 = qr.Q(qr(X), complete=TRUE)[,-(1:2)]; Zk = Zk%*%Q2; Omega = crossprod(Q2, Omega)%*%Q2}

  # SVD of penalty matrix
  # So that the "random effects" have diagonal prior variance
  svd.Omega = svd(Omega)
  sqrt.Omega = t(svd.Omega$v %*%(t(svd.Omega$u)*sqrt(svd.Omega$d)))
  Z = t(solve(sqrt.Omega,t(Zk)))

  # Now combine the linear and nonlinear pieces to obtain the matrix of basis functions evaluated at the obs. points
  Bmat = cbind(X, Z);

  # The penalty matrix:
  Omega = diag(c(rep(0, ncol(X)), rep(1, ncol(Z))))

  if(orthonormalize){
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
#' Compute the Nelson-Siegel Curves
#'
#' Compute the three Nelson-Siegel curves at a given vector of points
#'
#' @param tau the \code{m x 1} vector of observation points (in months)
#' @param lambda the nonlinear parameter
#' @return The \code{m x 3} matrix of Nelson-Siegel curves
#' evaluated at \code{tau}
#'
#' @examples
#' tau = 1:300
#' F_ns = f_ns(tau)
#' plot(tau, F_ns[,1], ylim = range(F_ns), type='l', lwd=4)
#' for(j in 2:3) lines(tau, F_ns[,j], lwd=4,lty=j)
#'
#' @export
f_ns = function(tau, lambda_p = 0.0609){
  f = matrix(1, nr =length(tau), nc = 3)
  f[,2] = (1 - exp(-lambda_p*tau))/(lambda_p*tau)
  f[,3] = ((1 - exp(-lambda_p*tau))/(lambda_p*tau) - exp(-lambda_p*tau))
  return(f)
}
#####################################################################################################
#' Compute Simultaneous Credible Bands
#'
#' Compute (1-alpha)\% credible BANDS for a function based on MCMC samples using Crainiceanu et al. (2007)
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
#' Compute Simultaneous Band Scores (SimBaS)
#'
#' Compute simultaneous band scores (SimBaS) from Meyer et al. (2015, Biometrics).
#' SimBaS uses MC(MC) simulations of a function of interest to compute the minimum  
#' alpha such that the joint credible bands at the alpha level do not include zero.
#' This quantity is computed for each grid point (or observation point) in the domain
#' of the function. 
#'
#' @param sampFuns \code{Nsims x m} matrix of \code{Nsims} MCMC samples and \code{m} points along the curve
#'
#' @return \code{m x 1} vector of simBaS
#'
#' @note The input needs not be curves: the simBaS may be computed
#' for vectors to achieve a multiplicity adjustment.
#' 
#' @note The minimum of the returned value, \code{PsimBaS_t}, 
#' over the domain \code{t} is the Global Bayesian P-Value (GBPV) for testing
#' whether the function is zero everywhere.
#'
#' @export
simBaS = function(sampFuns){
  
  N = nrow(sampFuns); m = ncol(sampFuns)
  
  # Compute pointwise mean and SD of f(x):
  Efx = colMeans(sampFuns); SDfx = apply(sampFuns, 2, sd)
  
  # Compute standardized absolute deviation:
  Standfx = abs(sampFuns - tcrossprod(rep(1, N), Efx))/tcrossprod(rep(1, N), SDfx)
  
  # And the maximum:
  Maxfx = apply(Standfx, 1, max)
  
  # And now compute the SimBaS scores:
  PsimBaS_t = rowMeans(sapply(Maxfx, function(x) abs(Efx)/SDfx <= x))
  
  # Alternatively, using a loop:
  #PsimBaS_t = numeric(T); for(t in 1:m) PsimBaS_t[t] = mean((abs(Efx)/SDfx)[t] <= Maxfx)
  
  PsimBaS_t
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
#' Compute a block diagonal matrix w/ constant blocks
#'
#' The function returns kronecker(diag(nrep), Amat), but is computed more efficiently
#' @param Amat matrix to populate the diagaonal blocks
#' @param nrep number of blocks on the diagonal
#----------------------------------------------------------------------------
blockDiag = function(Amat, nrep){
  nr1 = nrow(Amat); nc1 = ncol(Amat)
  fullMat = matrix(0, nr = nr1*nrep, nc = nc1*nrep)
  rSeq = seq(1, nr1*nrep + nr1, by=nr1) # row sequence
  cSeq = seq(1, nc1*nrep + nc1, by=nc1) # col sequence
  for(i in 1:nrep) fullMat[rSeq[i]:(rSeq[i+1] - 1),  cSeq[i]:(cSeq[i+1] - 1)] = Amat

  fullMat
}
#----------------------------------------------------------------------------
#' Plot a curve given posterior samples
#'
#' Plot the posterior mean, simultaneous and pointwise 95\% credible bands
#' for a curve given draws from the posterior distribution
#' @param post_f \code{Ns x m} matrix of \code{Ns} posterior simulations
#' of the curve at \code{m} points
#' @param tau \code{m x 1} vector of observation points
#' @param alpha confidence level for the bands
#' @param include_joint logical; if TRUE, include joint bands (as well as pointwise)
plot_curve = function(post_f, tau = NULL, alpha = 0.05, include_joint = TRUE){

  Ns = nrow(post_f); m = ncol(post_f)

  if(is.null(tau)) tau = 1:m

  par(mfrow = c(1, 1), mai = c(1, 1, 1, 1))

  # Pointwise intervals:
  dcip = dcib = HPDinterval(as.mcmc(post_f), prob = 1 - alpha);

  # Joint intervals, if necessary:
  if(include_joint) dcib = credBands(post_f, alpha = alpha)

  f_hat = colMeans(post_f)

  plot(tau, f_hat, type = "n", ylim = range(dcib, dcip, na.rm = TRUE),
       xlab = expression(tau), ylab = "", main = "Posterior Mean and Credible Bands",
       cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
  if(include_joint) polygon(c(tau, rev(tau)), c(dcib[, 2], rev(dcib[, 1])), col = "gray50",
                            border = NA)
  polygon(c(tau, rev(tau)), c(dcip[, 2], rev(dcip[, 1])), col = "grey",
          border = NA)
  lines(tau, f_hat, lwd = 8, col = "cyan")
}
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
    lines(dates,colMeans(post_beta[,,k]), lwd=8, col=k)
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
#' @import coda
#' @export
plot_flc = function(post_fk, tau = NULL){
  K = dim(post_fk)[3] # Number of factors
  if(is.null(tau)) tau = seq(0, 1, length.out = dim(post_fk)[2])

  dev.new(); par(mai = c(.9,.9,.4,.4), bg = 'gray90');
  plot(tau, post_fk[1,,1], ylim = range(post_fk), xlab = expression(tau), ylab = '', main = 'Factor Loading Curves', type='n', cex.lab = 2, cex.axis=2,cex.main=2)
  abline(h = 0, lty=3, lwd=2);
  for(k in K:1){
    # Credible intervals:
    ci = HPDinterval(as.mcmc(post_fk[,,k]));
    # Credible bands (w/ error catch):
    cb = try(credBands(as.mcmc(post_fk[,,k])), silent = TRUE)
    if(class(cb) == "try-error") cb = ci
    polygon(c(tau, rev(tau)), c(cb[,2], rev(cb[,1])), col='grey50', border=NA);
    polygon(c(tau, rev(tau)), c(ci[,2], rev(ci[,1])), col='grey', border=NA);
    lines(tau,colMeans(post_fk[,,k]), lwd=8, col=k)
  }
}
#----------------------------------------------------------------------------
#' Univariate Slice Sampler from Neal (2008)
#'
#' Compute a draw from a univariate distribution using the code provided by
#' Radford M. Neal. The documentation below is also reproduced from Neal (2008).
#'
#' @param x0    Initial point
#' @param g     Function returning the log of the probability density (plus constant)
#' @param w     Size of the steps for creating interval (default 1)
#' @param m     Limit on steps (default infinite)
#' @param lower Lower bound on support of the distribution (default -Inf)
#' @param upper Upper bound on support of the distribution (default +Inf)
#' @param gx0   Value of g(x0), if known (default is not known)
#'
#' @return  The point sampled, with its log density attached as an attribute.
#'
#' @note The log density function may return -Inf for points outside the support
#' of the distribution.  If a lower and/or upper bound is specified for the
#' support, the log density function will not be called outside such limits.
uni.slice <- function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL)
{
  # Check the validity of the arguments.

  if (!is.numeric(x0) || length(x0)!=1
      || !is.function(g)
      || !is.numeric(w) || length(w)!=1 || w<=0
      || !is.numeric(m) || !is.infinite(m) && (m<=0 || m>1e9 || floor(m)!=m)
      || !is.numeric(lower) || length(lower)!=1 || x0<lower
      || !is.numeric(upper) || length(upper)!=1 || x0>upper
      || upper<=lower
      || !is.null(gx0) && (!is.numeric(gx0) || length(gx0)!=1))
  {
    stop ("Invalid slice sampling argument")
  }

  # Keep track of the number of calls made to this function.
  #uni.slice.calls <<- uni.slice.calls + 1

  # Find the log density at the initial point, if not already known.

  if (is.null(gx0))
  { #uni.slice.evals <<- uni.slice.evals + 1
    gx0 <- g(x0)
  }

  # Determine the slice level, in log terms.

  logy <- gx0 - rexp(1)

  # Find the initial interval to sample from.

  u <- runif(1,0,w)
  L <- x0 - u
  R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff

  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.

  if (is.infinite(m))  # no limit on number of steps
  {
    repeat
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
    }

    repeat
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
    }
  }

  else if (m>1)  # limit on steps, bigger than one
  {
    J <- floor(runif(1,0,m))
    K <- (m-1) - J

    while (J>0)
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
      J <- J - 1
    }

    while (K>0)
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
      K <- K - 1
    }
  }

  # Shrink interval to lower and upper bounds.

  if (L<lower)
  { L <- lower
  }
  if (R>upper)
  { R <- upper
  }

  # Sample from the interval, shrinking it on each rejection.

  repeat
  {
    x1 <- runif(1,L,R)

    #uni.slice.evals <<- uni.slice.evals + 1
    gx1 <- g(x1)

    if (gx1>=logy) break

    if (x1>x0)
    { R <- x1
    }
    else
    { L <- x1
    }
  }

  # Return the point sampled, with its log density attached as an attribute.

  attr(x1,"log.density") <- gx1
  return (x1)

}

# Just add these for general use:
#' @importFrom stats quantile rgamma rnorm sd splinefun var rexp runif
NULL
