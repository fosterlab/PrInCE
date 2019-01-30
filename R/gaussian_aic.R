#' Model selection for Gaussian mixture models
#' 
#' Calculate the AIC, corrected AIC, or BIC for a curve fit with a Gaussian 
#' mixture model by nonlinear least squares optimization. This function 
#' permits the calculation of the AIC/AICc/BIC after rejecting some Gaussians
#' in the model, for example because their centres are outside the bounds of 
#' the profile.
#'
#' @param coefs the coefficients of the Gaussian mixture model, output by 
#'   \code{\link{fit_gaussians}}
#' @param chromatogram the raw elution profile
#' 
#' @return the AIC, corrected AIC, or BIC of the fit model
#' 
#' @name aic
NULL

#' @rdname aic
gaussian_aic <- function(coefs, chromatogram) {
  # first, calculate log likelihood
  N <- length(chromatogram)
  indices <- seq_len(N)
  fit <- fit_curve(coefs, indices)
  res <- chromatogram - fit
  w <- rep_len(1, N)
  zw <- w == 0
  loglik <- -N * (log(2 * pi) + 1 - log(N) - sum(log(w + zw)) + 
                    log(sum(w * res^2)))/2
  # next, calculate parameters
  k <- length(unlist(coefs)) + 1
  AIC <- 2 * k - 2 * loglik
  return(AIC)
}

#' @rdname aic
gaussian_aicc <- function(coefs, chromatogram) {
  # first, calculate AIC
  AIC <- gaussian_aic(coefs, chromatogram)
  # second, calculate AICc
  N <- length(chromatogram)
  k <- length(unlist(coefs)) + 1
  AICc <- AIC + (2 * k * (k + 1)) / (N - k - 1)
  return(AICc)
}

#' @rdname aic
gaussian_bic <- function(coefs, chromatogram) {
  # first, calculate log likelihood
  N <- length(chromatogram)
  indices <- seq_len(N)
  fit <- fit_curve(coefs, indices)
  res <- chromatogram - fit
  w <- rep_len(1, N)
  zw <- w == 0
  loglik <- -N * (log(2 * pi) + 1 - log(N) - sum(log(w + zw)) + 
                    log(sum(w * res^2)))/2
  # next, calculate parameters
  k <- length(unlist(coefs)) + 1
  BIC <- log(N) * k - 2 * loglik
  return(BIC)
}


