# Note: to update, use "git push -u origin master" (C******7)
#' MCMC Sampling Algorithm for the FDLM
#'
#' Runs the MCMC for the (univariate) FDLM under some default conditions:
#' \enumerate{
#' \item A random walk model for the dynamic factors;
#' \item A full \code{K x K} non-dynamic evolution error variance matrix;
#' \item A non-dynamic scalar multiple of the identity for the observation error variance.
#' }
#'
#' @param Y the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau vector of observation points (\code{m}-dimensional)
#' @param K the number of factors; if NULL, use SVD-based proportion of variability explained
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "beta" (factors)
#' \item "fk" (loading curves)
#' \item "sigma_et" (observation error SD; possibly dynamic)
#' \item "Wt" (evolution error variance)
#' \item "Yhat" (fitted values)
#' }
#' @param h_step integer for h-step forecasting; if NULL, do not compute any forecasts
#' @param useFastImpute logical; when TRUE, use imputation/projection scheme for the dynamic factors; otherwise use full state space model for factors (slower)
#' @param use_obs_SV logical; when TRUE, include a stochastic volatility model
#' for the observation error variance
#' @param includeBasisInnovation logical; when TRUE, include an iid basis coefficient term for residual correlation
#' (i.e., the idiosyncratic error term for a factor model on the full basis matrix)
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note If \code{Tm} is large, then storing all posterior samples for \code{Yhat}, which is \code{nsave x T x m},  may be inefficient
#'
#' @examples
#' # Read in the yield curve data (US or UK):
#' data("US_Yields") # data("UK_Yields")
#'
#' # Restrict to dates since 2012:
#' Y = Y[which(dates > as.Date("2012-01-01")),];
#' dates = dates[which(dates > as.Date("2012-01-01"))] # subset the dates for easy reference
#'
#' # Center and scale for numerical purposes:
#' Y = scale(Y)
#'
#' # Run the MCMC:
#' mcmc_output = fdlm(Y, tau, K = 3,
#'                   nsave = 1000, nburn = 100, nskip = 2,
#'                   mcmc_params = list("beta", "fk", "Yhat", "sigma_et"))
#'
#' # Plot the factors:
#' plot_factors(mcmc_output$beta, dates)
#'
#' # Plot the factor loading curves:
#' plot_flc(mcmc_output$fk, tau)
#'
#' # Some diagnostics: effective sample size(s) (may be slow!)
#' getEffSize(mcmc_output$beta)
#' getEffSize(mcmc_output$fk)
#' # getEffSize(mcmc_output$Yhat)
#'
#' # Check the residuals:
#' Yhat = colMeans(mcmc_output$Yhat)
#' resids = Y - Yhat
#' dev.new(); persp(tau, dates, t(resids))
#'
#' @export
fdlm = function(Y, tau, K = NULL,
                nsave = 1000, nburn = 1000, nskip = 10,
                mcmc_params = list("beta", "fk"),
                h_step = NULL,
                useFastImpute = TRUE,
                use_obs_SV = FALSE,
                includeBasisInnovation = TRUE,
                computeDIC = TRUE){

  # Compute the dimensions:
  T = nrow(Y); m = ncol(Y)

  # Rescale observation points to [0,1]
  tau01 = (tau - min(tau))/diff(range(tau))

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init(Y, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo

  # Also use the imputed data values here for initialization:
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  if(any.missing) Y = inits$Y0
  BtY = tcrossprod(t(splineInfo$Bmat), Y)

  # In case K was unspecified, and therefore chosen via SVD initialization:
  if(is.null(K)) K = ncol(Beta)

  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi

  # Initialize the conditional expectation:
  BetaPsit = tcrossprod(Beta, Psi); Btheta = tcrossprod(BetaPsit, splineInfo$Bmat)

  # Initialize the basis coefficient residuals and the corresponding standard deviation
  if(includeBasisInnovation){
    omega = t(BtY) - BetaPsit; sigma_w = sd(omega)
    theta = BetaPsit + omega; Btheta = tcrossprod(theta,splineInfo$Bmat)
  } else sigma_w = 0

  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Btheta)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Btheta, na.rm=TRUE)
    sigma_et = rep(sigma_e, T)
  }

  # Initialize the FLC smoothing parameters (conditional MLE):
  lambda = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - 2)/crossprod(x, splineInfo$Omega)%*%x)

  if(useFastImpute){

    # The response: Y projected onto FLCs
    YF = crossprod(BtY, Psi)

    # Initialize the SSModel:
    kfas_model = update_kfas_model(Y.dlm = YF, Zt = diag(K))

  } else kfas_model = update_kfas_model(Y.dlm = Y, Zt = Fmat) # Full DLM in non-fast case

  # Initialize the factor evolution error variance:
  Wt = array(var(diff(Beta)), c(K,K,T))

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('sigma_et', mcmc_params)) || computeDIC) post.sigma_et = array(NA, c(nsave, T))
  if(!is.na(match('sigma_w', mcmc_params))) post.sigma_w = array(NA, c(nsave, 1))
  if(!is.na(match('Wt', mcmc_params))) post.Wt = array(NA, c(nsave, K, K))
  if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat = array(NA, c(nsave, T, m))
  if(!is.null(h_step)) post.yfore = array(NA, c(nsave, h_step, m))
  if(computeDIC) post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    # Impute the data, Y:
    if(useFastImpute && any.missing){Y.samp = fdlm_impute(Yna, Btheta, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}

    # Sample the FLCs
    #Psi = sampleFLC(BtY, Beta, Psi, BtB = diag(nrow(BtY)), splineInfo$Omega, lambda, sigmat2 = sigma_et^2 + sigma_w^2)
    Psi = fdlm_flc(BtY = BtY,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = lambda,
                   sigmat2 = sigma_et^2 + sigma_w^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi; if(useFastImpute) YF = crossprod(BtY, Psi)

    # Sample the smoothing parameters:
    lambda = sample_lambda(lambda, Psi, Omega = splineInfo$Omega, uniformPrior = TRUE, orderLambdas = TRUE)

    # Sample the factors (note: some of these arguments are unnecessary)
    Beta = fdlm_factor(Y = Y,
                       sigma_et = sqrt(sigma_et^2 + sigma_w^2),
                       Wt = Wt,
                       Fmat = Fmat,
                       YF = YF,
                       Gt = NULL,
                       W0 = NULL,
                       kfas_model = kfas_model,
                       useFastImpute = useFastImpute)

    # Store this term:
    BetaPsit = tcrossprod(Beta,Psi)

    # Sample the basis coefficients:
    if(includeBasisInnovation){

      # Quad/linear construction a little faster w/o observation SV, but both work:
      if(use_obs_SV){
        Sigma_prec = matrix(rep(sigma_et^-2, ncol(theta)), nr = T)
        chQtheta = sqrt(Sigma_prec + sigma_w^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)*Sigma_prec + BetaPsit/sigma_w^2 # Linear term from the posterior
      } else {
        chQtheta = sqrt(sigma_e^-2 + sigma_w^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)/sigma_e^2 + BetaPsit/sigma_w^2 # Linear term from the posterior
      }
      theta = linTheta/chQtheta^2 + 1/chQtheta*rnorm(length(theta))
      Btheta = tcrossprod(theta,splineInfo$Bmat)
      sigma_w = 1/sqrt(rgamma(n = 1, shape = 0.001 + length(theta)/2, rate =0.001 + sum((theta - BetaPsit)^2)/2))
    } else {theta = BetaPsit; sigma_w = 0; Btheta = tcrossprod(theta,splineInfo$Bmat)}

    # Sample the observation error variance
    if(use_obs_SV){
      svParams = sampleCommonSV(Y - Btheta, svParams)
      sigma_et = svParams$sigma_et
    } else {
      sigma_e = 1/sqrt(rgamma(n = 1, shape = sum(!is.na(Y))/2, rate = sum((Y - Btheta)^2, na.rm=TRUE)/2))
      sigma_et = rep(sigma_e, T)
    }

    # Sample the evolution error variance, assuming random walk model:
    Wt[1:K, 1:K,] = sample_Wt(resBeta = diff(Beta), useDiagonal=FALSE)

    # Adjust the ordering:
    if(nsi == 10 && K > 1){adjOrder = order(lambda, decreasing = TRUE); lambda = lambda[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('sigma_et', mcmc_params)) || computeDIC) post.sigma_et[isave,] = sigma_et
        if(!is.na(match('sigma_w', mcmc_params))) post.sigma_w[isave,] = sigma_w
        if(!is.na(match('Wt', mcmc_params))) post.Wt[isave,,] = Wt[,,1]
        if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat[isave,,] = Btheta # + sigma_e*rnorm(length(Y))
        #if(!is.na(match('yfore', mcmc_params))) post.yfore[isave,] =Fmat%*%kfas_model$T[,,T]%*%Beta[T,]
        if(!is.null(h_step)) {
          GbetaT = kfas_model$T[,,T]%*%Beta[T,]
          post.yfore[isave,1, ] = Fmat%*%GbetaT
          if(h_step > 1){for(h in 2:h_step){
            GbetaT = kfas_model$T[,,T]%*%GbetaT
            post.yfore[isave,h, ] = Fmat%*%GbetaT
          }}
        }
        if(computeDIC) post_loglike[isave] = sum(dnorm(matrix(Yna), mean = matrix(Btheta), sd = rep(sigma_et,m), log = TRUE), na.rm = TRUE)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 500)
  }

  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('sigma_et', mcmc_params))) mcmc_output$sigma_et = post.sigma_et
  if(!is.na(match('sigma_w', mcmc_params))) mcmc_output$sigma_w = post.sigma_w
  if(!is.na(match('Wt', mcmc_params))) mcmc_output$Wt = post.Wt
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  if(!is.null(h_step)) mcmc_output$yfore = post.yfore

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(matrix(Yna),
                            mean = matrix(colMeans(post.Yhat)),
                            sd = rep(colMeans(post.sigma_et), m),
                            log = TRUE), na.rm=TRUE)

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }

  print(paste('Total time: ', round((proc.time()[3] - timer0)/60), 'minutes'))

  return (mcmc_output);
}
#' MCMC Sampling Algorithm for the FDLM
#'
#' Runs the MCMC for the (univariate) FDLM under some default conditions:
#' \enumerate{
#' \item A random walk model for the dynamic factors;
#' \item A full \code{K x K} non-dynamic evolution error variance matrix;
#' \item A non-dynamic scalar multiple of the identity for the observation error variance.
#' }
#'
#' @param Y the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau vector of observation points (\code{m}-dimensional)
#' @param K the number of factors; if NULL, use SVD-based proportion of variability explained
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "beta" (factors)
#' \item "fk" (loading curves)
#' \item "sigma_et" (observation error SD; possibly dynamic)
#' \item "Wt" (evolution error variance)
#' \item "Yhat" (fitted values)
#' }
#' @param h_step integer for h-step forecasting; if NULL, do not compute any forecasts
#' @param useFastImpute logical; when TRUE, use imputation/projection scheme for the dynamic factors; otherwise use full state space model for factors (slower)
#' @param use_obs_SV logical; when TRUE, include a stochastic volatility model
#' for the observation error variance
#' @param includeBasisInnovation logical; when TRUE, include an iid basis coefficient term for residual correlation
#' (i.e., the idiosyncratic error term for a factor model on the full basis matrix)
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#' @note If \code{Tm} is large, then storing all posterior samples for \code{Yhat}, which is \code{nsave x T x m},  may be inefficient
#'
#' @examples
#' # Read in the yield curve data (US or UK):
#' data("US_Yields") # data("UK_Yields")
#'
#' # Restrict to dates since 2012:
#' Y = Y[which(dates > as.Date("2012-01-01")),];
#' dates = dates[which(dates > as.Date("2012-01-01"))] # subset the dates for easy reference
#'
#' # Center and scale for numerical purposes:
#' Y = scale(Y)
#'
#' # Run the MCMC:
#' mcmc_output = fdlm0(Y, tau, K = 3,
#'                   nsave = 1000, nburn = 100, nskip = 2,
#'                   mcmc_params = list("beta", "fk", "Yhat", "sigma_e"))
#'
#' # Plot the factors:
#' plot_factors(mcmc_output$beta, dates)
#'
#' # Plot the factor loading curves:
#' plot_flc(mcmc_output$fk, tau)
#'
#' # Some diagnostics: effective sample size(s) (may be slow!)
#' getEffSize(mcmc_output$beta)
#' getEffSize(mcmc_output$fk)
#' # getEffSize(mcmc_output$Yhat)
#'
#' # Check the residuals:
#' Yhat = colMeans(mcmc_output$Yhat)
#' resids = Y - Yhat
#' dev.new(); persp(tau, dates, t(resids))
#'
#' @export
fdlm0 = function(Y, tau, K = NULL,
                nsave = 1000, nburn = 1000, nskip = 10,
                mcmc_params = list("beta", "fk"),
                h_step = NULL,
                useFastImpute = TRUE,
                use_obs_SV = FALSE,
                includeBasisInnovation = TRUE,
                computeDIC = TRUE){

  # Compute the dimensions:
  T = nrow(Y); m = ncol(Y)

  # Rescale observation points to [0,1]
  tau01 = (tau - min(tau))/diff(range(tau))

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init(Y, tau, K); Beta = inits$Beta; Psi = inits$Psi; splineInfo = inits$splineInfo

  # Also use the imputed data values here for initialization:
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  if(any.missing) Y = inits$Y0
  BtY = tcrossprod(t(splineInfo$Bmat), Y)

  # In case K was unspecified, and therefore chosen via SVD initialization:
  if(is.null(K)) K = ncol(Beta)

  # FLC matrix:
  Fmat = splineInfo$Bmat%*%Psi

  # Initialize the conditional expectation:
  BetaPsit = tcrossprod(Beta, Psi); Btheta = tcrossprod(BetaPsit, splineInfo$Bmat)

  # Initialize the basis coefficient residuals and the corresponding standard deviation
  if(includeBasisInnovation){
    omega = t(BtY) - BetaPsit; sigma_w = sd(omega)
    theta = BetaPsit + omega; Btheta = tcrossprod(theta,splineInfo$Bmat)
  } else sigma_w = 0

  # Initialize the (time-dependent) observation error SD:
  if(use_obs_SV){
    svParams = initCommonSV(Y - Btheta)
    sigma_et = svParams$sigma_et
  } else {
    sigma_e = sd(Y - Btheta, na.rm=TRUE)
    sigma_et = rep(sigma_e, T)
  }

  # Initialize the FLC smoothing parameters (conditional MLE):
  lambda = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - 2)/crossprod(x, splineInfo$Omega)%*%x)

  if(useFastImpute){

    # The response: Y projected onto FLCs
    YF = crossprod(BtY, Psi)

    # Initialize the SSModel:
    kfas_model = update_kfas_model(Y.dlm = YF, Zt = diag(K))

  } else kfas_model = update_kfas_model(Y.dlm = Y, Zt = Fmat) # Full DLM in non-fast case

  # Initialize the factor evolution error variance:
  Wt = array(var(diff(Beta)), c(K,K,T))

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('sigma_et', mcmc_params)) || computeDIC) post.sigma_et = array(NA, c(nsave, T))
  if(!is.na(match('sigma_w', mcmc_params))) post.sigma_w = array(NA, c(nsave, 1))
  if(!is.na(match('Wt', mcmc_params))) post.Wt = array(NA, c(nsave, K, K))
  if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat = array(NA, c(nsave, T, m))
  if(!is.null(h_step)) post.yfore = array(NA, c(nsave, m))
  if(computeDIC) post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    # Impute the data, Y:
    if(useFastImpute && any.missing){Y.samp = fdlm_impute(Yna, Btheta, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}

    # Sample the FLCs
    #Psi = sampleFLC(BtY, Beta, Psi, BtB = diag(nrow(BtY)), splineInfo$Omega, lambda, sigmat2 = rep(sigma_e^2 + sigma_w^2, T))
    Psi = fdlm_flc(BtY = BtY,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = lambda,
                   sigmat2 = sigma_et^2 + sigma_w^2)
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi; if(useFastImpute) YF = crossprod(BtY, Psi)

    # Sample the smoothing parameters:
    lambda = sample_lambda(lambda, Psi, Omega = splineInfo$Omega, uniformPrior = TRUE, orderLambdas = TRUE)

    # Sample the factors (note: some of these arguments are unnecessary)
    Beta = fdlm_factor(Y = Y,
                       sigma_et = sqrt(sigma_et^2 + sigma_w^2),
                       Wt = Wt,
                       Fmat = Fmat,
                       YF = YF,
                       Gt = NULL,
                       W0 = NULL,
                       kfas_model = kfas_model,
                       useFastImpute = useFastImpute)

    # Store this term:
    BetaPsit = tcrossprod(Beta,Psi)

    # Sample the basis coefficients:
    if(includeBasisInnovation){

      # Quad/linear construction a little faster w/o observation SV, but both work:
      if(use_obs_SV){
        Sigma_prec = matrix(rep(sigma_et^-2, ncol(theta)), nr = T)
        chQtheta = sqrt(Sigma_prec + sigma_w^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)*Sigma_prec + BetaPsit/sigma_w^2 # Linear term from the posterior
      } else {
        chQtheta = sqrt(sigma_e^-2 + sigma_w^-2) # Chol of diag (quadratic term) is just sqrt
        linTheta = t(BtY)/sigma_e^2 + BetaPsit/sigma_w^2 # Linear term from the posterior
      }
      theta = linTheta/chQtheta^2 + 1/chQtheta*rnorm(length(theta))
      Btheta = tcrossprod(theta,splineInfo$Bmat)
      sigma_w = 1/sqrt(rgamma(n = 1, shape = 0.001 + length(theta)/2, rate =0.001 + sum((theta - BetaPsit)^2)/2))
    } else {theta = BetaPsit; sigma_w = 0; Btheta = tcrossprod(theta,splineInfo$Bmat)}

    # Sample the observation error variance
    if(use_obs_SV){
      svParams = sampleCommonSV(Y - Btheta, svParams)
      sigma_et = svParams$sigma_et
    } else {
      sigma_e = 1/sqrt(rgamma(n = 1, shape = sum(!is.na(Y))/2, rate = sum((Y - Btheta)^2, na.rm=TRUE)/2))
      sigma_et = rep(sigma_e, T)
    }

    # Sample the evolution error variance, assuming random walk model:
    Wt[1:K, 1:K,] = sample_Wt(resBeta = diff(Beta), useDiagonal=FALSE)

    # Adjust the ordering:
    if(nsi == 10 && K > 1){adjOrder = order(lambda, decreasing = TRUE); lambda = lambda[adjOrder]; Psi = Psi[,adjOrder]; Beta = as.matrix(Beta[,adjOrder])}

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('beta', mcmc_params))) post.beta[isave,,] = Beta
        if(!is.na(match('fk', mcmc_params))) post.fk[isave,,] = Fmat
        if(!is.na(match('sigma_et', mcmc_params)) || computeDIC) post.sigma_et[isave,] = sigma_et
        if(!is.na(match('sigma_w', mcmc_params))) post.sigma_w[isave,] = sigma_w
        if(!is.na(match('Wt', mcmc_params))) post.Wt[isave,,] = Wt[,,1]
        if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat[isave,,] = Btheta # + sigma_e*rnorm(length(Y))
        #if(!is.na(match('yfore', mcmc_params))) post.yfore[isave,] =Fmat%*%kfas_model$T[,,T]%*%Beta[T,]
        if(!is.null(h_step)) {
          GbetaT = kfas_model$T[,,T]%*%Beta[T,]
          post.yfore[isave,1, ] = Fmat%*%GbetaT
          if(h_step > 1){for(h in 2:h_step){
            GbetaT = kfas_model$T[,,T]%*%GbetaT
            post.yfore[isave,h, ] = Fmat%*%GbetaT
          }}
        }
        if(computeDIC) post_loglike[isave] = sum(dnorm(matrix(Yna), mean = matrix(Btheta), sd = rep(sigma_et,m), log = TRUE), na.rm = TRUE)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 500)
  }

  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('sigma_et', mcmc_params))) mcmc_output$sigma_et = post.sigma_et
  if(!is.na(match('sigma_w', mcmc_params))) mcmc_output$sigma_w = post.sigma_w
  if(!is.na(match('Wt', mcmc_params))) mcmc_output$Wt = post.Wt
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  #if(!is.na(match('yfore', mcmc_params))) mcmc_output$yfore = post.yfore
  if(!is.null(h_step)) mcmc_output$yfore = post.yfore

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(matrix(Yna),
                            mean = matrix(colMeans(post.Yhat)),
                            sd = rep(colMeans(post.sigma_et), m),
                            log = TRUE), na.rm=TRUE)

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }

  print(paste('Total time: ', round((proc.time()[3] - timer0)/60), 'minutes'))

  return (mcmc_output);
}
