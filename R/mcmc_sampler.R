#' MCMC Sampling Algorithm for the FDLM
#'
#' Runs the MCMC for the (univariate) FDLM under some default conditions:
#' 1) A random walk model for the dynamic factors;
#' 2) A full \code{K x K} non-dynamic evolution error variance matrix;
#' 3) A non-dynamic scalar multiple of the identity for the observation error variance.
#'
#' @param Y the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau vector of observation points (\code{m}-dimensional)
#' @param K the number of factors; if NULL, use SVD-based proportion of variability explained
#' @param nsims number of MCMC iterations to record
#' @param burnin length of the burnin, which is discarded
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of "beta" (factors),  "fk" (FLCs), "sigma_e" (observation error SD), "Wt" (evolution error variance), "Yhat" (fitted values)
#' @param useFastImpute logical; when TRUE, use imputation/projection scheme for the dynamic factors; otherwise use full state space model for factors (slower)
#' @return A named list of the \code{nsims} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note If \code{Tm} is large, then storing all posterior samples for \code{Yhat}, which is \code{nsims x T x m},  may be inefficient
#'
#' @examples
#' # Read in the yield curve data:
#' data("US_Yields")
#'
#' # Restrict to dates since 2006:
#' Y = Y[which(dates > as.Date("2006-01-01")),];
#'
#' # Run the MCMC:
#' mcmc_output = fdlm(Y, tau, K = 3, nsims = 1000)
#' @export
fdlm = function(Y, tau, K = NULL, nsims = 10000, burnin = 1000, mcmc_params = list("beta", "fk"), useFastImpute = TRUE){

  # Compute the dimensions:
  T = nrow(Y); m = ncol(Y)

  # Rescale observation points to [0,1]
  tau01 = (tau - min(tau))/diff(range(tau))

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init(Y, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo

  # In case K was unspecified, and therefore chosen via SVD initialization:
  if(is.null(K)) K = ncol(Beta)

  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi

  # Initialize the conditional expectation:
  BetaPsit = tcrossprod(Beta, Psi); Btheta = tcrossprod(BetaPsit, splineInfo$Bmat)

  # Initialize the missing obs, as well as the cross products BtY (imputed) and BtY0 (only computed from observed data)
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  Y.samp = fdlm_impute(Yna, Btheta, sigma_et = rep(0,T), Bmat = splineInfo$Bmat)
  Y = Y.samp$Y; BtY = Y.samp$BtY; BtY0 = Y.samp$BtY0

  # Initialize the basis coefficient residuals and the corresponding standard deviation
  includeBasisInnovation = TRUE
  if(includeBasisInnovation){
    omega = t(BtY) - BetaPsit; sigma_w = sd(omega)
    theta = BetaPsit + omega; Btheta = tcrossprod(theta,splineInfo$Bmat)
  } else sigma_w = 0

  # Initialize the (time-dependent) observation error SD:
  sigma_e = sd(Y - Btheta, na.rm=TRUE)
  sigma_et = rep(sigma_e, T)

  # Initialize the FLC smoothing parameters (conditional MLE):
  lambda = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - 2)/crossprod(x, splineInfo$Omega)%*%x)

  if(useFastImpute){

    # The response: Y projected onto FLCs
    YF = crossprod(BtY, Psi)

    # Initialize the SSModel:
    kfas_model = update.kfas_model(Y.dlm = YF, Zt = diag(K))

  } else kfas_model = update.kfas_model(Y.dlm = Y, Zt = Fmat) # Full DLM in non-fast case

  # Initialize the factor evolution error variance:
  Wt = array(var(diff(Beta)), c(K,K,T))

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsims, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsims, m, K))
  if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e = array(NA, c(nsims, 1))
  if(!is.na(match('sigma_w', mcmc_params))) post.sigma_w = array(NA, c(nsims, 1))
  if(!is.na(match('Wt', mcmc_params))) post.Wt = array(NA, c(nsims, K, K))
  if(!is.na(match('Yhat', mcmc_params))) post.Yhat = array(NA, c(nsims, T, m))

  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:(nsims + burnin)){

    # Impute the data, Y:
    if(useFastImpute && any.missing){Y.samp = fdlm_impute(Yna, Btheta, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}

    # Sample the FLCs
    Psi = sampleFLC(BtY, Beta, Psi, BtB = diag(nrow(BtY)), splineInfo$Omega, lambda, sigmat2 = rep(sigma_e^2 + sigma_w^2, T))
    Fmat = splineInfo$Bmat%*%Psi; if(useFastImpute) YF = crossprod(BtY, Psi)

    # Sample the smoothing parameters:
    lambda = sample_lambda(lambda, Psi, Omega = splineInfo$Omega, uniformPrior = TRUE, orderLambdas = TRUE)

    # Sample the factors (note: some of these arguments are unnecessary)
    Beta = fdlm_factor(Y = Y, sigma_et = rep(sqrt(sigma_e^2 + sigma_w^2), T), Wt = Wt,  Fmat = Fmat, YF = YF, Gt = NULL, kfas_model = kfas_model, useFastImpute = useFastImpute)

    # Store this term:
    BetaPsit = tcrossprod(Beta,Psi)

    # Sample the basis coefficients:
    if(includeBasisInnovation){
      chQtheta = sqrt(sigma_e^-2 + sigma_w^-2) # Chol of diag (quadratic term) is just sqrt
      linTheta = t(BtY)/sigma_e^2 + BetaPsit/sigma_w^2 # Linear term from the posterior
      theta = linTheta/chQtheta^2 + 1/chQtheta*rnorm(length(theta))
      Btheta = tcrossprod(theta,splineInfo$Bmat)
      sigma_w = 1/sqrt(rgamma(n = 1, shape = 0.001 + length(theta)/2, rate =0.001 + sum((theta - BetaPsit)^2)/2))
    } else {theta = BetaPsit; sigma_w = 0; Btheta = tcrossprod(theta,splineInfo$Bmat)}

    # Sample the observation error variance (just assume Jeffreys prior)
    sigma_e = 1/sqrt(rgamma(n = 1, shape = sum(!is.na(Y))/2, rate = sum((Y - Btheta)^2, na.rm=TRUE)/2))
    sigma_et = rep(sigma_e, T)

    # Sample the evolution error variance, assuming random walk model:
    Wt[1:K, 1:K,] = sample_Wt(resBeta = diff(Beta), useDiagonal=FALSE)

    # Adjust the ordering:
    if(nsi == 10){adjOrder = order(lambda, decreasing = TRUE); lambda = lambda[adjOrder]; Psi = Psi[,adjOrder]; Beta = Beta[,adjOrder]}

    # Store the MCMC output:
    if(nsi > burnin){
      if(!is.na(match('beta', mcmc_params))) post.beta[nsi - burnin,,] = Beta
      if(!is.na(match('fk', mcmc_params))) post.fk[nsi - burnin,,] = Fmat
      if(!is.na(match('sigma_e', mcmc_params))) post.sigma_e[nsi-burnin,] = sigma_e
      if(!is.na(match('sigma_w', mcmc_params))) post.sigma_w[nsi-burnin,] = sigma_w
      if(!is.na(match('Wt', mcmc_params))) post.Wt[nsi-burnin,,] = Wt[,,1]
      if(!is.na(match('Yhat', mcmc_params))) post.Yhat[nsi-burnin,,] = Btheta # + sigma_e*rnorm(length(Y))
    }
    computeTimeRemaining(nsi, timer0, (nsims + burnin), nrep = 500)
  }

  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('sigma_w', mcmc_params))) mcmc_output$sigma_w = post.sigma_w
  if(!is.na(match('Wt', mcmc_params))) mcmc_output$Wt = post.Wt
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat

  print(paste('Total time: ', round((proc.time()[3] - timer0)/60), 'minutes'))

  return (mcmc_output);
}
