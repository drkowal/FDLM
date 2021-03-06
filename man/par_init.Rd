% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{par_init}
\alias{par_init}
\title{Initialize the parametric terms}
\usage{
par_init(Y, tau, f_p, orthogonalize = TRUE)
}
\arguments{
\item{Y}{the \code{T x m} data observation matrix, where \code{T} is the number of time points and \code{m} is the number of observation points (\code{NA}s allowed)}

\item{tau}{vector of observation points (\code{m}-dimensional)}

\item{f_p}{a function to compute the parametric component, which must return a \code{m x K_p} matrix
for \code{K_p} the number of parametric curves; may include a (scalar) nonlinear parameter argument}

\item{orthogonalize}{logical; when TRUE, orthogonalize the loading curve matrix}
}
\value{
a list containing
\itemize{
\item \code{Beta_p} the \code{T x K_p} matrix of parametric factors
\item \code{F_p} the \code{m x K_p} matrix of parametric loading curves
\item \code{lambda_p} the scalar nonlinear parameter
}
}
\description{
Compute initial values for the factors and nonlinear parameter (if necessary)
of the parametric component.
}
\details{
Compute initial values via the following algorithm:
\enumerate{
\item Impute missing values in \code{Y}
\item Initialize \code{lambda_p = 1} and compute \code{F_p}; orthogonalize if specified
\item Estimate \code{Beta_p} via least squares
\item Estimate an initial standard deviation \code{sigma_0} via conditional MLE
\item Estimate \code{lambda_p} via conditional MLE and recompute \code{F_p}; orthogonalize if specified
}
}
