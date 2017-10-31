#' Calculate the corrected AIC for a NLS fit
#' 
#' Calculates the corrected AIC for a curve fit with a Gaussian mixture model
#' by nonlinear least squares optimization, using the base R \code{nls}
#' function. 
#' 
#' @param fit the curve fit by the \code{stats::nls} function
#' 
#' @return the corrected AIC of the fit model
nls_aicc <- function(fit) {
  n <- length(fitted(fit))
  loglik <- as.numeric(logLik(fit))
  k <- length(coef(fit)) + 1
  AIC <- 2 * k - 2 * loglik
  AICc <- AIC + (2 * k * (k + 1)) / (n - k - 1)
  return(AICc)
}
