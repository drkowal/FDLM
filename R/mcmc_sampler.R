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
#' \item "sigma_e" (observation error SD)
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
#' # Run the MCMC:
#' mcmc_output = fdlm(Y, tau, K = 3,
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
fdlm = function(Y, tau, K = NULL,
                nsave = 1000, nburn = 1000, nskip = 10,
                mcmc_params = list("beta", "fk"),
                h_step = NULL,
                useFastImpute = TRUE,
                use_obs_SV = FALSE,
                includeBasisInnovation = TRUE,
                computeDIC = TRUE){

  # Check the model specifications to see if they make sense:
  if(use_obs_SV) stop("SV not yet implemented")

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
  sigma_e = sd(Y - Btheta, na.rm=TRUE)
  sigma_et = rep(sigma_e, T)

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
  if(!is.na(match('sigma_e', mcmc_params)) || computeDIC) post.sigma_e = array(NA, c(nsave, 1))
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
    #Psi = sampleFLC(BtY, Beta, Psi, BtB = diag(nrow(BtY)), splineInfo$Omega, lambda, sigmat2 = rep(sigma_e^2 + sigma_w^2, T))
    Psi = fdlm_flc(BtY = BtY,
                   Beta  = Beta,
                   Psi = Psi,
                   BtB = diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   lambda = lambda,
                   sigmat2 = rep(sigma_e^2 + sigma_w^2, T))
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi; if(useFastImpute) YF = crossprod(BtY, Psi)

    # Sample the smoothing parameters:
    lambda = sample_lambda(lambda, Psi, Omega = splineInfo$Omega, uniformPrior = TRUE, orderLambdas = TRUE)

    # Sample the factors (note: some of these arguments are unnecessary)
    Beta = fdlm_factor(Y = Y,
                       sigma_et = rep(sqrt(sigma_e^2 + sigma_w^2), T),
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
        if(!is.na(match('sigma_e', mcmc_params)) || computeDIC) post.sigma_e[isave,] = sigma_e
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
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('sigma_w', mcmc_params))) mcmc_output$sigma_w = post.sigma_w
  if(!is.na(match('Wt', mcmc_params))) mcmc_output$Wt = post.Wt
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  #if(!is.null(h_step)) mcmc_output$yfore = post.yfore
  if(!is.null(h_step)) mcmc_output$yfore = post.yfore

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(matrix(Yna),
                            mean = matrix(colMeans(post.Yhat)),
                            sd = colMeans(post.sigma_e),
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
#' MCMC Sampling Algorithm for the Semiparametric Functional Dynamic Linear Model
#'
#' Runs the MCMC for the (univariate) Semiparametric FDLM. The options for the dynamic parametric factors are:
#' \enumerate{
#' \item Random walk model with full evolution covariance matrix
#' \item Independent AR(1) models (i.e., diagional covariance matrix)
#' \item Full vector autoregression (VAR) coefficient matrix with full evolution covariance matrix
#'}
#' The dynamic nonparametric factors follow a random walk model with the following options for the evolution error:
#' \enumerate{
#' \item The dynamic horseshoe prior
#' \item The static horseshoe prior
#' \item The normal-inverse-gamma prior
#' }
#' In each case, the dynamic nonparametric factors are independent, and independent from the dynamic parametric factors.
#'
#'
#' @param Y the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau vector of observation points (\code{m}-dimensional)
#' @param f a function to compute the parametric component, which must return a \code{m x K_p} matrix
#' for \code{K_p} the number of parametric curves; may include a (scalar) nonlinear parameter argument
#' @param K_np the number of (nonparametric) factors; if NULL, use SVD-based proportion of variability explained
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "beta" (factors)
#' \item "fk" (loading curves)
#' \item "sigma_e" (observation error SD)
#' \item "lambda_p" (parametric nonlinear parameter)
#' \item "Wt" (evolution error variance)
#' \item "Yhat" (fitted values)
#' }
#' @param h_step integer for h-step forecasting; if NULL, do not compute any forecasts
#' @param evol_error_par string denoting the model for the parametric factor evolution;
#' must be one of
#' \itemize{
#' \item "RW": random walk with full evolution covariance matrix
#' \item "AR": independent AR(1) models for the factors
#' \item "VAR": full VAR coefficient matrix with full evolution covariance matrix
#' }
#'@param evol_error_par string denoting the model for the nonparametric factor
#' evolution error; must be one of
#' \itemize{
#' \item "DHS": the dynamic horseshoe prior
#' \item "HS": the static horseshoe prior
#' \item "NIG": the normal-inverse-gamma prior
#' }
#' @param use_obs_SV logical; when TRUE, include a stochastic volatility model
#' for the observation error variance
#' @param includeBasisInnovation logical; when TRUE, include an iid basis coefficient term for residual correlation
#' (i.e., the idiosyncratic error term for a factor model on the full basis matrix)
#' @param log_prior_lambda_p a function to evaluate the log-prior for the nonlinear
#' parametric component; if \code{NULL}, do not sample the nonlinear component
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#'
#' @details  The parametric function \code{f} should input an \code{m}-dimensional
#' vector of observation points, \code{tau}, and may include a (known or unknown) nonlinear
#' parameter, \code{lambda_p}. The function should return a \code{m x K_n} matrix, where \code{K_n} is the
#' number of parametric functions. For example, \code{f = function(tau) cbind(1,tau)}
#' includes an intercept and a linear term (\code{K_n = 2}). If the parametric function
#' includes a nonlinear term, for example, \code{f = function(tau, lambda_p) cbind(1,exp(-tau/lambda_p))},
#' then supply a (log) prior function via \code{log_prior_lambda_p} will allow for sampling of this
#' parameter. However, if \code{log_prior_lambda_p} is \code{NULL}, then the nonlinear parameter
#' will be fixed at its initialized value.
#'
#' @note The parametric factors are assumed to be independent of the nonparametric factors.
#' The evolution model for the nonparametric factors is a random walk.
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
#' f = function(tau, lambda_p) f_ns(tau, lambda_p)[,1:2]
#' mcmc_output = sfdlm(Y, tau, f, K_np = 2,
#'                    nsave = 1000, nburn = 100, nskip = 2,
#'                    mcmc_params = list("beta", "fk", "Yhat", "sigma_e", "lambda_p"),
#'                    evol_error_par = 'RW', evol_error_nonpar = "DHS")
#' # Plot the factors:
#' plot_factors(mcmc_output$beta, dates)
#'
#' # Plot the factor loading curves:
#' plot_flc(mcmc_output$fk, tau)
#'
#' # Some diagnostics: effective sample size(s) (may be slow!)
#' getEffSize(mcmc_output$beta)
#' getEffSize(mcmc_output$fk[,,-(1:length(f(1,1)))])
#' # getEffSize(mcmc_output$Yhat)
#'
#' # Check the residuals:
#' Yhat = colMeans(mcmc_output$Yhat)
#' resids = Y - Yhat
#' dev.new(); persp(tau, dates, t(resids))
#'
#' # Example: Nelson-Siegel model w/ random nonlinear parameter
#'mcmc_output = sfdlm(Y, tau, f = f_ns, K_np = 1,
#'                    nsave = 1000, nburn = 100, nskip = 2,
#'                    mcmc_params = list("beta", "fk", "Yhat", "sigma_e", "lambda_p"),
#'                    includeBasisInnovation = FALSE,
#'                    log_prior_lambda_p = function(x) dunif(x, min = 10^-4, max = 10^4, log=TRUE))
#'
#' # Example: Intercept and exponential decay
#' mcmc_output = sfdlm(Y, tau, f = function(tau0, x) cbind(1, exp(-tau0/x)), K_np = 2,
#'                    nsave = 1000, nburn = 100, nskip = 2,
#'                    mcmc_params = list("beta", "fk", "Yhat", "sigma_e", "lambda_p"))
#'
#' @import dsp
#' @import vars
#' @export
sfdlm = function(Y, tau, f, K_np = NULL,
                 nsave = 1000, nburn = 1000, nskip = 10,
                 mcmc_params = list("beta", "fk"),
                 h_step = NULL,
                 evol_error_par = "RW",
                 evol_error_nonpar = "DHS",
                 use_obs_SV = FALSE,
                 includeBasisInnovation = TRUE,
                 log_prior_lambda_p = NULL,
                 computeDIC = TRUE){

  # Check the model specifications to see if they make sense:
  if(use_obs_SV) stop("SV not yet implemented")

  # Convert to upper case, then check for matches to existing models:
  evol_error_par = toupper(evol_error_par); evol_error_nonpar = toupper(evol_error_nonpar)
  if(is.na(match(evol_error_par, c("RW", "AR", "VAR"))))
    stop("The parametric evolution error must be one of 'RW', 'AR', or 'VAR'")
  if(is.na(match(evol_error_nonpar, c("DHS", "HS", "NIG"))))
    stop("The nonparametric evolution error must be one of 'DHS', 'HS', or 'NIG'")

  # Sample the nonlinear parameter only if a (log) prior has been supplied
  sampleNonLinearParam = !is.null(log_prior_lambda_p)

  # Redefine the input function to have a silent nonlinear input, if necessary:
  f_try = try(f(tau, 1), silent = TRUE)
  if(class(f_try) == "try-error") {
    # No need to sample the nonlinear parameter, since it does not exist in this case
    sampleNonLinearParam = FALSE
    f_p = function(tau, lambda_p) f(tau)
  } else f_p = f

  # Use this to determine whether or not we need to orthogonalize the N-S matrix:
  orthogonalize = !includeBasisInnovation
  #orthogonalize = TRUE

  # Compute the dimensions:
  T = nrow(Y); m = ncol(Y)

  # Rescale observation points to [0,1]
  tau01 = (tau - min(tau))/diff(range(tau))

  # Initialize the parametric components: try state space model
  inits0 = try(par_init_ssm(Y, tau, f_p, orthogonalize), silent = TRUE)
  if(class(inits0) == "try-error") inits0 = par_init(Y, tau, f_p, orthogonalize);
  Beta_p = as.matrix(inits0$Beta_p); F_p = as.matrix(inits0$F_p); lambda_p = inits0$lambda_p

  # Number of parametric terms:
  K_p = ncol(F_p)

  # Parametric term for conditional mean
  Yhat_p = tcrossprod(Beta_p, F_p)

  # Note: if only 1 parametric term, use AR instead of VAR
  if(K_p == 1 && evol_error_par == "VAR") evol_error_par = "AR"

  # Initialize the FLC coefficients and factors:
  inits = fdlm_init(Y - Yhat_p, tau, K_np); Beta_fdlm = as.matrix(inits$Beta); Psi = as.matrix(inits$Psi); splineInfo = inits$splineInfo
  K_np = ncol(Psi)
  F_fdlm = splineInfo$Bmat%*%Psi

  # Also use the imputed data values here for initialization:
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  if(any.missing) Y = Yhat_p + inits$Y0
  BtY = tcrossprod(t(splineInfo$Bmat), Y)

  # Total number of factors:
  K = K_p + K_np

  # All factors and FLCs:
  Beta = matrix(0, nrow = T, ncol = K);
  Beta[,1:K_p] = Beta_p; Beta[,-(1:K_p)] = Beta_fdlm
  Fmat = matrix(0, nrow = m, ncol = K);
  Fmat[,1:K_p] = F_p; Fmat[,-(1:K_p)] = F_fdlm

  # Initialize the conditional expectation of the FDLM part:
  BetaPsit = tcrossprod(Beta_fdlm, Psi); Btheta = tcrossprod(BetaPsit, splineInfo$Bmat)

  # Conditional mean:
  Yhat = Yhat_p + Btheta

  # Initialize the missing obs, as well as the cross products BtY (imputed)
  #Yna = Y # The original data, including NAs
  #any.missing = any(is.na(Yna)) # Any missing obs?
  #Y.samp = fdlm_impute(Yna, Yhat, sigma_et = rep(0,T), Bmat = splineInfo$Bmat)
  #Y = Y.samp$Y; BtY = Y.samp$BtY

  # Initialize the basis coefficient residuals and the corresponding standard deviation
  BtY_p = tcrossprod(t(splineInfo$Bmat), Yhat_p)
  if(includeBasisInnovation){
    omega = t(BtY - BtY_p) - BetaPsit; sigma_w = sd(omega)
    theta = BetaPsit + omega; Btheta = tcrossprod(theta,splineInfo$Bmat)

    Yhat = Yhat_p + Btheta

  } else sigma_w = 0

  # Initialize the (time-dependent) observation error SD:
  sigma_e = sd(Y - Yhat, na.rm=TRUE)
  sigma_et = rep(sigma_e, T)

  # Initialize the FLC smoothing parameters (conditional MLE):
  lambda = apply(Psi, 2, function(x) (ncol(splineInfo$Bmat) - 2)/crossprod(x, splineInfo$Omega)%*%x)

  # Initialize the KFAS model:
  if(includeBasisInnovation){

    # In this case, project the data (and N-S curves) on Bmat
    # Note: Bmat is (m x J) and m >> J >> K, so it's faster than just using Y
    kfas_model = update_kfas_model(Y.dlm = t(BtY),
                                   Zt = cbind(crossprod(splineInfo$Bmat, F_p), Psi))

  } else {
    # The response: Y projected onto NS+FLCs (requires Fmat O.N.)
    YF = Y%*%Fmat #t(tcrossprod(t(Fmat), Y))

    # Initialize the SSModel:
    kfas_model = update_kfas_model(Y.dlm = YF, Zt = diag(K))
  }

  # Now additional the evolution parameters
  G_mu = array(diag(K), c(K,K, T)) # Evolution matrix
  Wt = array(0, c(K,K, T)) # Variance matrix
  mu_all = matrix(0, nrow = T, ncol = K)
  # Parametric terms:
  if(evol_error_par == "RW"){
    # Define the mean to be zero:
    mu_alpha = matrix(0, nrow = K_p, ncol = 1)
    # Evolution matrix is diagional:
    G_alpha = diag(K_p)
    # Full variance:
    Wt[1:K_p, 1:K_p, ] = var(diff(Beta_p))
  } else {
    # Unconditional mean (easy version)
    mu_alpha = as.matrix(colMeans(Beta_p))

    # Centered par factors:
    Beta_p_cent = Beta_p - matrix(rep(mu_alpha, each =  T), nrow = T)

    if(evol_error_par == "AR"){
      # AR(1) coefficients:
      G_alpha = diag(apply(Beta_p_cent, 2, function(x) lm(x[-1] ~ - 1+  x[-length(x)])$coef), K_p)

      # Stationarity fix:
      G_alpha[which(abs(G_alpha) > 0.99, arr.ind = TRUE)] = 0.8

      # Initialize the variance:
      Beta_p_resid = Beta_p_cent[-1,] - t(tcrossprod(G_alpha, as.matrix(Beta_p_cent[-T,])))
      Wt[1:K_p, 1:K_p, ] = diag(diag(var(Beta_p_resid)), K_p)
    } else {
      # VAR
      G_alpha = matrix(unlist(lapply(
        VAR(Beta_p_cent, p=1, "none")$varresult, coef)),
        nrow=K_p, byrow=TRUE)

      # Initialize the variance:
      Beta_p_resid = Beta_p_cent[-1,] - t(tcrossprod(G_alpha, as.matrix(Beta_p_cent[-T,])))
      Wt[1:K_p, 1:K_p, ] = var(Beta_p_resid)
    }
  }
  mu_all[,1:K_p] = rep(mu_alpha, each = T)
  G_mu[1:K_p, 1:K_p,] = G_alpha

  # Nonparametric terms:
  evolParams = initEvolParams(omega = diff(Beta_fdlm), evol_error = evol_error_nonpar)
  evolParams0 = initEvol0(mu0 = as.matrix(Beta_fdlm[1,]))
  for(k in 1:K_np) Wt[K_p + k,K_p + k,-T] = evolParams$sigma_wt[,k]^2

  # Initial variance:
  W0 = diag(10^-4, K); diag(W0)[-(1:K_p)] = evolParams0$sigma_w0^2

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('sigma_e', mcmc_params)) || computeDIC) post.sigma_e = array(NA, c(nsave, 1))
  if(!is.na(match('lambda_p', mcmc_params))) post.lambda_p = array(NA, c(nsave, 1))
  if(!is.na(match('sigma_w', mcmc_params))) post.sigma_w = array(NA, c(nsave, 1))
  if(!is.na(match('Wt', mcmc_params))) post.Wt = array(NA, c(nsave, K, K))
  if(!is.na(match('G_alpha', mcmc_params))) post.G_alpha = array(NA, c(nsave, K_p, K_p))
  if(!is.na(match('mu_alpha', mcmc_params))) post.mu_alpha = array(NA, c(nsave, K_p))
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
    if(any.missing){Y.samp = fdlm_impute(Yna, Yhat, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}

    # Sample the N-S nonlinear parameter
    if(sampleNonLinearParam){
      Yres = Y - tcrossprod(Beta_fdlm, F_fdlm)
      lambda_p = uni.slice(lambda_p, g = function(x){
        F_p_x = f_p(tau, x); if(orthogonalize) F_p_x = qr.Q(qr(F_p_x))
        sum(-0.5*rowSums((tcrossprod(Beta_p, F_p_x) - Yres)^2)/sigma_et^2) +
          log_prior_lambda_p(x)
      })

      # Redefine F_ns:
      F_p = f_p(tau, lambda_p); if(orthogonalize) F_p = qr.Q(qr(F_p))
    }

    # Sample the FLCs
    Psi = fdlm_flc(BtY = BtY - BtY_p,
                   Beta  = Beta_fdlm,
                   Psi = Psi,
                   BtB = diag(nrow(BtY)),
                   Omega = splineInfo$Omega,
                   BtCon = crossprod(splineInfo$Bmat, F_p),
                   lambda = lambda,
                   sigmat2 = rep(sigma_e^2 + sigma_w^2, T))
    # And compute the loading curves:
    F_fdlm = splineInfo$Bmat%*%Psi
    Fmat[,1:K_p] = F_p; Fmat[,-(1:K_p)] = F_fdlm

    # Sample the smoothing parameters:
    lambda = sample_lambda(lambda, Psi, Omega = splineInfo$Omega, uniformPrior = TRUE, orderLambdas = TRUE)

    # Sample the factors
    # Note: the arguments may or may not be used, depending on includeBasisInnovation (logical)
    BtFp = crossprod(splineInfo$Bmat, F_p)
    Beta = mu_all + fdlm_factor(Y = t(BtY) - tcrossprod(mu_all[,1:K_p], BtFp),
                       sigma_et = rep(sqrt(sigma_e^2 + sigma_w^2), T),
                       Wt = Wt,
                       Fmat = cbind(BtFp, Psi),
                       YF = Y%*%Fmat - mu_all,
                       Gt = G_mu,
                       W0 = W0,
                       kfas_model = kfas_model,
                       useFastImpute = !includeBasisInnovation)
    # And store the components:
    Beta_p = as.matrix(Beta[,1:K_p]); Beta_fdlm = as.matrix(Beta[,-(1:K_p)])

    # Update the Yhat term:
    Yhat_p = tcrossprod(Beta_p, F_p)
    BtY_p = tcrossprod(t(splineInfo$Bmat), Yhat_p)

    # Store this term as well:
    BetaPsit = tcrossprod(Beta_fdlm,Psi)

    # Sample the basis coefficients:
    if(includeBasisInnovation){
      chQtheta = sqrt(sigma_e^-2 + sigma_w^-2) # Chol of diag (quadratic term) is just sqrt
      linTheta = t(BtY - BtY_p)/sigma_e^2 + BetaPsit/sigma_w^2 # Linear term from the posterior
      theta = linTheta/chQtheta^2 + 1/chQtheta*rnorm(length(theta))
      Btheta = tcrossprod(theta,splineInfo$Bmat)
      sigma_w = 1/sqrt(rgamma(n = 1, shape = 0.001 + length(theta)/2, rate =0.001 + sum((theta - BetaPsit)^2)/2))

    } else {theta = BetaPsit; sigma_w = 0; Btheta = tcrossprod(theta,splineInfo$Bmat)}

    Yhat = Yhat_p + Btheta

    # Sample the observation error variance (just assume Jeffreys prior)
    sigma_e = 1/sqrt(rgamma(n = 1, shape = sum(!is.na(Y))/2, rate = sum((Y - Yhat)^2, na.rm=TRUE)/2))
    sigma_et = rep(sigma_e, T)

    # Sample the evolution parameters
    # Parametric terms:
    if(evol_error_par == "RW") {
      Wt[1:K_p, 1:K_p,] = sample_Wt(resBeta = diff(Beta_p), useDiagonal=FALSE)
    } else {
      # AR(1) or VAR(1)
      # Sample the unconditional means:
      mu_alpha = sampleARmu(yt = Beta_p, G = G_alpha, Sigma = Wt[1:K_p, 1:K_p, 1])

      # Centered par factors:
      Beta_p_cent = Beta_p - matrix(rep(mu_alpha, each =  T), nrow = T)

      # Special case for AR(1) coefs and evolution error variance
      if(evol_error_par == "AR"){
        # Sample the AR(1) coefficients:
        diag(G_alpha) = sampleAR1(h_yc = Beta_p_cent,
                                  h_phi = diag(G_alpha),
                                  #h_sigma_eta_t = t(apply(Wt[1:K_p,1:K_p,-T], 3, diag)),
                                  h_sigma_eta_t = t(apply(array(Wt[1:K_p,1:K_p,-T], c(K_p, K_p, T-1)), 3, diag)),
                                  prior_dhs_phi = c(5,2))
        # And the variance:
        Beta_p_resid = Beta_p_cent[-1,] - t(tcrossprod(G_alpha, as.matrix(Beta_p_cent[-T,])))
        Wt[1:K_p, 1:K_p, ] = sample_Wt(resBeta = Beta_p_resid, useDiagonal=TRUE)
      }

      if(evol_error_par == "VAR"){
        # VAR Sampler:
        G_alpha = sampleVAR(ytc = Beta_p_cent, Sigma = Wt[1:K_p, 1:K_p, 1])

        # Sample the variance:
        Beta_p_resid = Beta_p_cent[-1,] - t(tcrossprod(G_alpha, as.matrix(Beta_p_cent[-T,])))
        Wt[1:K_p, 1:K_p, ] = sample_Wt(resBeta = Beta_p_resid, useDiagonal=FALSE)
      }
      # Update these terms:
      mu_all[,1:K_p] = rep(mu_alpha, each = T); G_mu[1:K_p, 1:K_p,] = G_alpha
    }

    # Nonparametric terms:
    evolParams = sampleEvolParams(omega = diff(Beta_fdlm),
                                  evolParams, 1/sqrt(T), evol_error = evol_error_nonpar)
    evolParams0 = sampleEvol0(mu0 = as.matrix(Beta_fdlm[1,]), evolParams0, A = 1)
    for(k in 1:K_np) Wt[K_p + k,K_p + k,-T] = evolParams$sigma_wt[,k]^2
    diag(W0)[-(1:K_p)] = evolParams0$sigma_w0^2

    # Adjust the ordering:
    if(nsi == 10 && K_np > 1){adjOrder = order(lambda, decreasing = TRUE); lambda = lambda[adjOrder]; Psi = Psi[,adjOrder]; Beta_fdlm= as.matrix(Beta_fdlm[,adjOrder]); Beta[,-(1:K_p)] = Beta_fdlm}

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
        if(!is.na(match('sigma_e', mcmc_params)) || computeDIC) post.sigma_e[isave,] = sigma_e
        if(!is.na(match('lambda_p', mcmc_params))) post.lambda_p[isave,] = lambda_p
        if(!is.na(match('sigma_w', mcmc_params))) post.sigma_w[isave,] = sigma_w
        if(!is.na(match('Wt', mcmc_params))) post.Wt[isave,,] = Wt[,,1]
        if(!is.na(match('G_alpha', mcmc_params))) post.G_alpha[isave,,] = G_alpha
        if(!is.na(match('mu_alpha', mcmc_params))) post.mu_alpha[isave,] = mu_alpha
        if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat[isave,,] = Yhat # + sigma_e*rnorm(length(Y))
        # Assumes random walk for FDLM (nonparametric) factors:
        #if(!is.na(match('yfore', mcmc_params))) post.yfore[isave,] = F_p%*%(mu_alpha + G_alpha%*%(Beta_p[T,] - mu_alpha)) + F_fdlm%*%Beta_fdlm[T,]
        if(!is.null(h_step)) {
          GbetaT = G_alpha%*%(Beta_p[T,] - mu_alpha)
          post.yfore[isave,1, ] = F_p%*%(mu_alpha + GbetaT) + F_fdlm%*%Beta_fdlm[T,]
          if(h_step > 1){for(h in 2:h_step){
            GbetaT = G_alpha%*%GbetaT
            post.yfore[isave,h, ] = F_p%*%(mu_alpha + GbetaT) + F_fdlm%*%Beta_fdlm[T,]
          }}
        }
        if(computeDIC) post_loglike[isave] = sum(dnorm(matrix(Yna), mean = matrix(Yhat), sd = rep(sigma_et,m), log = TRUE), na.rm = TRUE)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 500)
  }

  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('lambda_p', mcmc_params))) mcmc_output$lambda_p = post.lambda_p
  if(!is.na(match('sigma_w', mcmc_params))) mcmc_output$sigma_w = post.sigma_w
  if(!is.na(match('Wt', mcmc_params))) mcmc_output$Wt = post.Wt
  if(!is.na(match('G_alpha', mcmc_params))) mcmc_output$G_alpha = post.G_alpha
  if(!is.na(match('mu_alpha', mcmc_params))) mcmc_output$mu_alpha = post.mu_alpha
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  #if(!is.na(match('yfore', mcmc_params))) mcmc_output$yfore = post.yfore
  if(!is.null(h_step)) mcmc_output$yfore = post.yfore

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(matrix(Yna),
                            mean = matrix(colMeans(post.Yhat)),
                            sd = colMeans(post.sigma_e),
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
#' MCMC Sampling Algorithm for the Parametric Functional Dynamic Linear Model
#'
#' Runs the MCMC for the (univariate) parametric FDLM. The options for the dynamic (parametric) factors are:
#' \enumerate{
#' \item Random walk model with full evolution covariance matrix
#' \item Independent AR(1) models (i.e., diagional covariance matrix)
#' \item Full vector autoregression (VAR) coefficient matrix with full evolution covariance matrix
#'}
#'
#' @param Y the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)
#' @param tau vector of observation points (\code{m}-dimensional)
#' @param f a function to compute the parametric component, which must return a \code{m x K_p} matrix
#' for \code{K_p} the number of parametric curves; may include a (scalar) nonlinear parameter argument
#' @param nsave number of MCMC iterations to record
#' @param nburn number of MCMC iterations to discard (burin-in)
#' @param nskip number of MCMC iterations to skip between saving iterations,
#' i.e., save every (nskip + 1)th draw
#' @param mcmc_params named list of parameters for which we store the MCMC output;
#' must be one or more of
#' \itemize{
#' \item "beta" (factors)
#' \item "fk" (loading curves)
#' \item "sigma_e" (observation error SD)
#' \item "lambda_p" (parametric nonlinear parameter)
#' \item "Wt" (evolution error variance)
#' \item "Yhat" (fitted values)
#' }
#' @param h_step integer for h-step forecasting; if NULL, do not compute any forecasts
#' @param evol_error_par string denoting the model for the parametric factor evolution;
#' must be one of
#' \itemize{
#' \item "RW": random walk with full evolution covariance matrix
#' \item "AR": independent AR(1) models for the factors
#' \item "VAR": full VAR coefficient matrix with full evolution covariance matrix
#' }
#' @param use_obs_SV logical; when TRUE, include a stochastic volatility model
#' for the observation error variance
#' @param orthogonalize logical; when TRUE, orthogonalize the parametric loadings matrix
#' @param log_prior_lambda_p a function to evaluate the log-prior for the nonlinear
#' parametric component; if \code{NULL}, do not sample the nonlinear component
#' @param computeDIC logical; if TRUE, compute the deviance information criterion \code{DIC}
#' and the effective number of parameters \code{p_d}
#'
#' @return A named list of the \code{nsave} MCMC samples for the parameters named in \code{mcmc_params}
#'
#'
#' @details  The parametric function \code{f} should input an \code{m}-dimensional
#' vector of observation points, \code{tau}, and may include a (known or unknown) nonlinear
#' parameter, \code{lambda_p}. The function should return a \code{m x K_n} matrix, where \code{K_n} is the
#' number of parametric functions. For example, \code{f = function(tau) cbind(1,tau)}
#' includes an intercept and a linear term (\code{K_n = 2}). If the parametric function
#' includes a nonlinear term, for example, \code{f = function(tau, lambda_p) cbind(1,exp(-tau/lambda_p))},
#' then supply a (log) prior function via \code{log_prior_lambda_p} will allow for sampling of this
#' parameter. However, if \code{log_prior_lambda_p} is \code{NULL}, then the nonlinear parameter
#' will be fixed at its initialized value.
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
#' f = function(tau, lambda_p) f_ns(tau, lambda_p)[,1:2]
#' mcmc_output = pfdlm(Y, tau, f,
#'                    nsave = 1000, nburn = 100, nskip = 2,
#'                    mcmc_params = list("beta", "fk", "Yhat"),
#'                    evol_error_par = 'VAR')
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
#' @import dsp
#' @import vars
#' @export
pfdlm = function(Y, tau, f,
                 nsave = 1000, nburn = 1000, nskip = 10,
                 mcmc_params = list("beta", "fk"),
                 h_step = NULL,
                 evol_error_par = "RW",
                 use_obs_SV = FALSE,
                 orthogonalize = TRUE,
                 log_prior_lambda_p = NULL,
                 computeDIC = TRUE){

  # Check the model specifications to see if they make sense:
  if(use_obs_SV) stop("SV not yet implemented")

  # Convert to upper case, then check for matches to existing models:
  evol_error_par = toupper(evol_error_par);

  if(is.na(match(evol_error_par, c("RW", "AR", "VAR"))))
    stop("The parametric evolution error must be one of 'RW', 'AR', or 'VAR'")

  # Sample the nonlinear parameter only if a (log) prior has been supplied
  sampleNonLinearParam = !is.null(log_prior_lambda_p)

  # Redefine the input function to have a silent nonlinear input, if necessary:
  f_try = try(f(tau, 1), silent = TRUE)
  if(class(f_try) == "try-error") {
    # No need to sample the nonlinear parameter, since it does not exist in this case
    sampleNonLinearParam = FALSE
    f_p = function(tau, lambda_p) f(tau)
  } else f_p = f

  # Compute the dimensions:
  T = nrow(Y); m = ncol(Y)

  # Rescale observation points to [0,1]
  tau01 = (tau - min(tau))/diff(range(tau))

  # Initialize the parametric components: try state space model
  inits0 = try(par_init_ssm(Y, tau, f_p, orthogonalize), silent = TRUE)
  if(class(inits0) == "try-error") inits0 = par_init(Y, tau, f_p, orthogonalize);
  Beta_p = as.matrix(inits0$Beta_p); F_p = as.matrix(inits0$F_p); lambda_p = inits0$lambda_p

  # Necessary term:
  splineInfo = getSplineInfo(tau01)

  # Number of parametric terms:
  K_p = ncol(F_p)

  # Parametric term for conditional mean
  Yhat_p = tcrossprod(Beta_p, F_p)

  # Note: if only 1 parametric term, use AR instead of VAR
  if(K_p == 1 && evol_error_par == "VAR") evol_error_par = "AR"

  # Total number of factors:
  K = K_p

  # All factors and FLCs:
  Beta = Beta_p; Fmat = F_p

  # Conditional mean:
  Yhat = Yhat_p

  # Initialize the missing obs, as well as the cross products BtY (imputed)
  Yna = Y # The original data, including NAs
  any.missing = any(is.na(Yna)) # Any missing obs?
  Y.samp = fdlm_impute(Yna, Yhat, sigma_et = rep(0,T), Bmat = splineInfo$Bmat)
  Y = Y.samp$Y; BtY = Y.samp$BtY

  # Initialize the (time-dependent) observation error SD:
  sigma_e = sd(Y - Yhat, na.rm=TRUE)
  sigma_et = rep(sigma_e, T)

  # Initialize the KFAS model:
  if(!orthogonalize){

    # In this case, project the data (and N-S curves) on Bmat
    # Note: Bmat is (m x J) and m >> J >> K, so it's faster than just using Y
    kfas_model = update_kfas_model(Y.dlm = t(BtY),
                                   Zt = crossprod(splineInfo$Bmat, F_p))

  } else {
    # The response: Y projected onto NS+FLCs (requires Fmat O.N.)
    YF = Y%*%Fmat #t(tcrossprod(t(Fmat), Y))

    # Initialize the SSModel:
    kfas_model = update_kfas_model(Y.dlm = YF, Zt = diag(K))
  }

  # Now additional the evolution parameters
  G_mu = array(diag(K), c(K,K, T)) # Evolution matrix
  Wt = array(0, c(K,K, T)) # Variance matrix
  mu_all = matrix(0, nrow = T, ncol = K)
  # Parametric terms:
  if(evol_error_par == "RW"){
    # Define the mean to be zero:
    mu_alpha = matrix(0, nrow = K_p, ncol = 1)
    # Evolution matrix is diagional:
    G_alpha = diag(K_p)
    # Full variance:
    Wt[1:K_p, 1:K_p, ] = var(diff(Beta_p))
  } else {
    # Unconditional mean (easy version)
    mu_alpha = as.matrix(colMeans(Beta_p))

    # Centered par factors:
    Beta_p_cent = Beta_p - matrix(rep(mu_alpha, each =  T), nrow = T)

    if(evol_error_par == "AR"){
      # AR(1) coefficients:
      G_alpha = diag(apply(Beta_p_cent, 2, function(x) lm(x[-1] ~ - 1+  x[-length(x)])$coef), K_p)
      # Stationarity fix:
      G_alpha[which(abs(G_alpha) > 1, arr.ind = TRUE)] = 0.8

      # Initialize the variance:
      Beta_p_resid = Beta_p_cent[-1,] - t(tcrossprod(G_alpha, as.matrix(Beta_p_cent[-T,])))
      Wt[1:K_p, 1:K_p, ] = diag(diag(var(Beta_p_resid)), K_p)
    } else {
      # VAR
      G_alpha = matrix(unlist(lapply(
        VAR(Beta_p_cent, p=1, "none")$varresult, coef)),
        nrow=K_p, byrow=TRUE)

      # Initialize the variance:
      Beta_p_resid = Beta_p_cent[-1,] - t(tcrossprod(G_alpha, as.matrix(Beta_p_cent[-T,])))
      Wt[1:K_p, 1:K_p, ] = var(Beta_p_resid)
    }
  }
  mu_all[,1:K_p] = rep(mu_alpha, each = T)
  G_mu[1:K_p, 1:K_p,] = G_alpha

  # Initial variance:
  W0 = diag(10^-4, K)

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('beta', mcmc_params))) post.beta = array(NA, c(nsave, T, K))
  if(!is.na(match('fk', mcmc_params))) post.fk = array(NA, c(nsave, m, K))
  if(!is.na(match('sigma_e', mcmc_params)) || computeDIC) post.sigma_e = array(NA, c(nsave, 1))
  if(!is.na(match('lambda_p', mcmc_params))) post.lambda_p = array(NA, c(nsave, 1))
  if(!is.na(match('Wt', mcmc_params))) post.Wt = array(NA, c(nsave, K, K))
  if(!is.na(match('G_alpha', mcmc_params))) post.G_alpha = array(NA, c(nsave, K_p, K_p))
  if(!is.na(match('mu_alpha', mcmc_params))) post.mu_alpha = array(NA, c(nsave, K_p))
  if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat = array(NA, c(nsave, T, m))
  #if(!is.na(match('yfore', mcmc_params))) post.yfore = array(NA, c(nsave, m))
  if(!is.null(h_step)) post.yfore = array(NA, c(nsave, h_step, m))
  if(computeDIC) post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    # Impute the data, Y:
    if(any.missing){Y.samp = fdlm_impute(Yna, Yhat, sigma_et = sigma_et, Bmat = splineInfo$Bmat); Y = Y.samp$Y; BtY = Y.samp$BtY}

    # Sample the N-S nonlinear parameter
    if(sampleNonLinearParam){
      Yres = Y
      lambda_p = uni.slice(lambda_p, g = function(x){
        F_p_x = f_p(tau, x); if(orthogonalize) F_p_x = qr.Q(qr(F_p_x))
        sum(-0.5*rowSums((tcrossprod(Beta_p, F_p_x) - Yres)^2)/sigma_et^2) +
          log_prior_lambda_p(x)
      })

      # Redefine F_ns:
      F_p = f_p(tau, lambda_p); if(orthogonalize) F_p = qr.Q(qr(F_p))
    }
    Fmat = F_p

    # Sample the factors
    # Note: the arguments may or may not be used, depending on includeBasisInnovation (logical)
    BtFp = crossprod(splineInfo$Bmat, F_p)
    Beta = mu_all + fdlm_factor(Y = t(BtY) - tcrossprod(mu_all[,1:K_p], BtFp),
                                sigma_et = rep(sigma_e, T),
                                Wt = Wt,
                                Fmat = BtFp,
                                YF = Y%*%Fmat - mu_all,
                                Gt = G_mu,
                                W0 = W0,
                                kfas_model = kfas_model,
                                useFastImpute = orthogonalize)
    # And store the components:
    Beta_p = as.matrix(Beta)

    # Update the Yhat term:
    Yhat_p = tcrossprod(Beta_p, F_p)
    BtY_p = tcrossprod(t(splineInfo$Bmat), Yhat_p)

    Yhat = Yhat_p

    # Sample the observation error variance (just assume Jeffreys prior)
    sigma_e = 1/sqrt(rgamma(n = 1, shape = sum(!is.na(Y))/2, rate = sum((Y - Yhat)^2, na.rm=TRUE)/2))
    sigma_et = rep(sigma_e, T)

    # Sample the evolution parameters
    # Parametric terms:
    if(evol_error_par == "RW") {
      Wt[1:K_p, 1:K_p,] = sample_Wt(resBeta = diff(Beta_p), useDiagonal=FALSE)
    } else {
      # AR(1) or VAR(1)
      # Sample the unconditional means:
      mu_alpha = sampleARmu(yt = Beta_p, G = G_alpha, Sigma = Wt[1:K_p, 1:K_p, 1])

      # Centered par factors:
      Beta_p_cent = Beta_p - matrix(rep(mu_alpha, each =  T), nrow = T)

      # Special case for AR(1) coefs and evolution error variance
      if(evol_error_par == "AR"){
        # Sample the AR(1) coefficients:
        diag(G_alpha) = sampleAR1(h_yc = Beta_p_cent,
                                  h_phi = diag(G_alpha),
                                  #h_sigma_eta_t = t(apply(Wt[1:K_p,1:K_p,-T], 3, diag)),
                                  h_sigma_eta_t = t(apply(array(Wt[1:K_p,1:K_p,-T], c(K_p, K_p, T-1)), 3, diag)),
                                  prior_dhs_phi = c(5,2))
        # And the variance:
        Beta_p_resid = Beta_p_cent[-1,] - t(tcrossprod(G_alpha, as.matrix(Beta_p_cent[-T,])))
        Wt[1:K_p, 1:K_p, ] = sample_Wt(resBeta = Beta_p_resid, useDiagonal=TRUE)
      }

      if(evol_error_par == "VAR"){
        # VAR Sampler:
        G_alpha = sampleVAR(ytc = Beta_p_cent, Sigma = Wt[1:K_p, 1:K_p, 1])

        # Sample the variance:
        Beta_p_resid = Beta_p_cent[-1,] - t(tcrossprod(G_alpha, as.matrix(Beta_p_cent[-T,])))
        Wt[1:K_p, 1:K_p, ] = sample_Wt(resBeta = Beta_p_resid, useDiagonal=FALSE)
      }
      # Update these terms:
      mu_all[,1:K_p] = rep(mu_alpha, each = T); G_mu[1:K_p, 1:K_p,] = G_alpha
    }

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
        if(!is.na(match('sigma_e', mcmc_params)) || computeDIC) post.sigma_e[isave,] = sigma_e
        if(!is.na(match('lambda_p', mcmc_params))) post.lambda_p[isave,] = lambda_p
        if(!is.na(match('Wt', mcmc_params))) post.Wt[isave,,] = Wt[,,1]
        if(!is.na(match('G_alpha', mcmc_params))) post.G_alpha[isave,,] = G_alpha
        if(!is.na(match('mu_alpha', mcmc_params))) post.mu_alpha[isave,] = mu_alpha
        if(!is.na(match('Yhat', mcmc_params)) || computeDIC) post.Yhat[isave,,] = Yhat # + sigma_e*rnorm(length(Y))
        #if(!is.na(match('yfore', mcmc_params))) post.yfore[isave,] = F_p%*%(mu_alpha + G_alpha%*%(Beta_p[T,] - mu_alpha)) + F_fdlm%*%Beta_fdlm[T,]
        if(!is.null(h_step)) {
          GbetaT = G_alpha%*%(Beta_p[T,] - mu_alpha)
          post.yfore[isave,1, ] = F_p%*%(mu_alpha + GbetaT) + F_fdlm%*%Beta_fdlm[T,]
          if(h_step > 1){for(h in 2:h_step){
            GbetaT = G_alpha%*%GbetaT
            post.yfore[isave,h, ] = F_p%*%(mu_alpha + GbetaT) + F_fdlm%*%Beta_fdlm[T,]
          }}
        }
        if(computeDIC) post_loglike[isave] = sum(dnorm(matrix(Yna), mean = matrix(Yhat), sd = rep(sigma_et,m), log = TRUE), na.rm = TRUE)

        # And reset the skip counter:
        skipcount = 0
      }
    }
    computeTimeRemaining(nsi, timer0, nstot, nrep = 500)
  }

  if(!is.na(match('beta', mcmc_params))) mcmc_output$beta = post.beta
  if(!is.na(match('fk', mcmc_params))) mcmc_output$fk = post.fk
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('lambda_p', mcmc_params))) mcmc_output$lambda_p = post.lambda_p
  if(!is.na(match('Wt', mcmc_params))) mcmc_output$Wt = post.Wt
  if(!is.na(match('G_alpha', mcmc_params))) mcmc_output$G_alpha = post.G_alpha
  if(!is.na(match('mu_alpha', mcmc_params))) mcmc_output$mu_alpha = post.mu_alpha
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  #if(!is.na(match('yfore', mcmc_params))) mcmc_output$yfore = post.yfore
  if(!is.null(h_step)) mcmc_output$yfore = post.yfore

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(matrix(Yna),
                            mean = matrix(colMeans(post.Yhat)),
                            sd = colMeans(post.sigma_e),
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
#' \item "sigma_e" (observation error SD)
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

  # Check the model specifications to see if they make sense:
  if(use_obs_SV) stop("SV not yet implemented")

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
  sigma_e = sd(Y - Btheta, na.rm=TRUE)
  sigma_et = rep(sigma_e, T)

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
  if(!is.na(match('sigma_e', mcmc_params)) || computeDIC) post.sigma_e = array(NA, c(nsave, 1))
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
                   sigmat2 = rep(sigma_e^2 + sigma_w^2, T))
    # And update the loading curves:
    Fmat = splineInfo$Bmat%*%Psi; if(useFastImpute) YF = crossprod(BtY, Psi)

    # Sample the smoothing parameters:
    lambda = sample_lambda(lambda, Psi, Omega = splineInfo$Omega, uniformPrior = TRUE, orderLambdas = TRUE)

    # Sample the factors (note: some of these arguments are unnecessary)
    Beta = fdlm_factor(Y = Y,
                       sigma_et = rep(sqrt(sigma_e^2 + sigma_w^2), T),
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
        if(!is.na(match('sigma_e', mcmc_params)) || computeDIC) post.sigma_e[isave,] = sigma_e
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
  if(!is.na(match('sigma_e', mcmc_params))) mcmc_output$sigma_e = post.sigma_e
  if(!is.na(match('sigma_w', mcmc_params))) mcmc_output$sigma_w = post.sigma_w
  if(!is.na(match('Wt', mcmc_params))) mcmc_output$Wt = post.Wt
  if(!is.na(match('Yhat', mcmc_params))) mcmc_output$Yhat = post.Yhat
  #if(!is.na(match('yfore', mcmc_params))) mcmc_output$yfore = post.yfore
  if(!is.null(h_step)) mcmc_output$yfore = post.yfore

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(dnorm(matrix(Yna),
                            mean = matrix(colMeans(post.Yhat)),
                            sd = colMeans(post.sigma_e),
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
