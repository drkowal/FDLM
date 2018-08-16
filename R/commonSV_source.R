#----------------------------------------------------------------------------
#' Sample the (common) stochastic volatility parameters
#'
#' Compute one draw for each of the parameters in the stochastic volatility model
#' where each component of omega_{j,t} has a common dynamic variance.
#'
#' @param omega \code{T x m} matrix of residuals
#' @param svParams list of parameters to be updated (see Value below)
#' @param prior_phi the parameters of the prior for the log-volatilty AR(1) coefficient \code{h_phi};
#' either \code{NULL} for uniform on [-1,1] or a 2-dimensional vector of (shape1, shape2) for a Beta prior
#' on \code{[(h_phi + 1)/2]}
#' @return List of relevant components:
#' \itemize{
#' \item the \code{T x 1} error standard deviations \code{sigma_et},
#' \item the \code{T x 1} log-volatility \code{ht},
#' \item the log-vol unconditional mean \code{h_mu},
#' \item the log-vol AR(1) coefficient \code{h_phi},
#' \item the log-vol innovation standard deviation \code{h_sigma_eta}
#' }
#'
#' @import truncdist
#' @export
sampleCommonSV = function(omega, svParams, prior_phi = c(20,1.5)){

  # Store the SV parameters locally:
  ht = svParams$ht; h_mu = svParams$h_mu; h_phi = svParams$h_phi; h_sigma_eta = svParams$h_sigma_eta;

  # "Local" number of time points
  ht = as.matrix(ht)
  n = nrow(ht); m = ncol(ht)

  # Sample the log-volatilities using AWOL sampler
  ht = sampleCommonLogVols(h_y = omega, h_prev = ht, h_mu = h_mu, h_phi=h_phi, h_sigma_eta = h_sigma_eta)

  # Compute centered log-vols for the samplers below:
  ht_tilde = ht - h_mu

  # Sample AR(1) parameters
  h_phi = sampleAR1(h_yc = ht_tilde, h_phi = h_phi, h_sigma_eta_t = matrix(rep(h_sigma_eta, n-1)), prior_dhs_phi = prior_phi)

  # Sample evolution error: uniform prior on standard deviation
  eta_t = ht_tilde[-1] - h_phi*ht_tilde[-n]       # Residuals
  #h_sigma_eta = 1/sqrt(rgamma(n = 1, shape = 0.01 + length(eta_t)/2, rate = 0.01 + sum(eta_t^2)/2))
  h_sigma_eta = 1/sqrt(truncdist::rtrunc(n = 1,
                                         'gamma',
                                         a = 1/100^2, # Lower interval
                                         b = Inf,   # Upper interval
                                         shape = length(eta_t)/2 - 1/2,
                                         rate =  sum(eta_t^2)/2))

  # Sample the unconditional mean:
    # Prior mean: h_mu ~ N(-10, 100)
  y_mu = (ht[-1] - h_phi*ht[-n])/h_sigma_eta
  x_mu = (1 - h_phi)/h_sigma_eta
  postSD = 1/sqrt((n-1)*x_mu^2 + 1/100)
  postMean = (sum(x_mu*y_mu) + -10/100)*postSD^2
  h_mu = rnorm(n = 1, mean = postMean, sd = postSD)

  # Evolution error SD:
  sigma_et = exp(ht/2)

  # Note: if we have enormous SDs, reduce:
  sigma_et[which(sigma_et > 10^3, arr.ind = TRUE)] = 10^3

  # Return the same list, but with the new values
  list(sigma_et = sigma_et, ht = ht, h_mu = h_mu, h_phi = h_phi, h_sigma_eta = h_sigma_eta)
}
#----------------------------------------------------------------------------
#' Sample the latent log-volatilities, common to m-dimensional time series
#'
#' Compute one draw of the log-volatilities using a discrete mixture of Gaussians
#' approximation to the likelihood (see Omori, Chib, Shephard, and Nakajima, 2007)
#' where the log-vols are assumed to follow an AR(1) model. The model assumes that
#' the volatility are common to an m-dimensional time series.
#'
#' @param h_y the \code{T x m} matrix of data, which follows one join SV model
#' @param h_prev the \code{T x 1} matrix of the previous log-vols
#' @param h_mu the log-vol unconditional mean
#' @param h_phi the log-vol AR(1) coefficient
#' @param h_sigma_eta the log-vol innovation standard deviation
#'
#' @return \code{T x 1} matrix of simulated log-vols
#' @import Matrix
#' @export
sampleCommonLogVols = function(h_y, h_prev, h_mu, h_phi, h_sigma_eta){

  # Compute dimensions:
  h_y = as.matrix(h_y) # Just to be sure (T x m)
  n = nrow(h_y); m = ncol(h_y)

  h_prev = as.matrix(h_prev)

  # Mixture params: mean, variance, and weights
  # Kim, Shephard, Chib (1998) 7-component mixture:
  #m_st  = c(-11.40039, -5.24321, -9.83726, 1.50746,  -0.65098, 0.52478,  -2.35859)
  #v_st2 = c(5.795960,  2.613690, 5.179500, 0.167350, 0.640090, 0.340230, 1.262610)
  #q     = c(0.007300,  0.105560, 0.000020, 0.043950, 0.340010, 0.245660, 0.257500)

  # Omori, Chib, Shephard, Nakajima (2007) 10-component mixture:
  m_st  = c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000)
  v_st2 = c(0.11265, 0.17788, 0.26768, 0.40611,  0.62699,  0.98583,  1.57469,  2.54498,  4.16591,   7.33342)
  q     = c(0.00609, 0.04775, 0.13057, 0.20674,  0.22715,  0.18842,  0.12047,  0.05591,  0.01575,   0.00115)

  # Add an offset: common for all times, but distict for each j=1,...,p
  #yoffset = any(h_y^2 < 10^-16)*max(10^-8, mad(h_y)/10^6)
  yoffset = any(h_y==0)*sd(h_y)/10000

  # This is the response in our DLM, log(y^2)
  ystar = log(h_y^2 + yoffset)

  # Sample the mixture components
  #z = draw.indicators(res = ystar - matrix(rep(h_prev, m), nr = n), nmix = list(m = m_st, v = v_st2, p = q))
  z = sapply(ystar-ystar-matrix(rep(h_prev, m), nr = n), ncind, m_st, sqrt(v_st2), q)

  # Subset mean and variances to the sampled mixture components; (n x p) matrices
  m_st_all = matrix(m_st[z], nr=n); v_st2_all = matrix(v_st2[z], nr=n)

  # Evolution error variance:
  if(length(h_sigma_eta) == 1){
    sigmat2 = rep(h_sigma_eta^2, n);
    # for simplicity, ignore that part:
    #sigmat2[1] = h_sigma_eta^2/(1 - h_phi^2)
  } else sigmat2 = h_sigma_eta^2

  # Quadratic term:
  QHt.Matrix = bandSparse(n, k = c(0,1),
                          diag = list(rowSums(1/v_st2_all) + 1/sigmat2 + c(h_phi^2/sigmat2[-1], 0),
                                      -h_phi/sigmat2[-1]), symm = TRUE)
  chQht_Matrix = Matrix::chol(QHt.Matrix)

  # Linear term:
  linht = matrix(rowSums((ystar - m_st_all - h_mu)/v_st2_all))

  # Sample the log-vols:
  hsamp = h_mu + matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(n)), nr = n)

  # Return the (uncentered) log-vols
  hsamp
}
#----------------------------------------------------------------------------
#' Initialize the (common) SV parameters
#'
#' Compute initial values for common stochastic volatility  parameters
#'
#' @param omega \code{T x m} matrix of residuals
#' @return List of relevant components:
#' \itemize{
#' \item the \code{T x 1} error standard deviations \code{sigma_et},
#' \item the \code{T x 1} log-volatility \code{ht},
#' \item the log-vol unconditional mean \code{h_mu},
#' \item the log-vol AR(1) coefficient \code{h_phi},
#' \item the log-vol innovation standard deviation \code{h_sigma_eta}
#' }
#' @export
initCommonSV = function(omega){

  # "Local" number of time points
  omega = as.matrix(omega)
  n = nrow(omega); m = ncol(omega)

  # Initialize the log-volatilities:
  ht = rowMeans(log(omega^2 + 0.0001))

  # Initialize the AR(1) model to obtain unconditional mean and AR(1) coefficient
  arCoefs = arima(ht, c(1,0,0))$coef
  h_mu = arCoefs[2]; h_phi = arCoefs[1]

  # Initialize the SD of log-vol innovations simply using the expectation:
  h_sigma_eta = sd((ht - h_mu)[-1] - h_phi*(ht - h_mu)[-n] )

  # Evolution error SD:
  sigma_et = exp(ht/2)

  # Return the same list, but with the new values
  list(sigma_et = sigma_et, ht = ht, h_mu = h_mu, h_phi = h_phi, h_sigma_eta = h_sigma_eta)
}
