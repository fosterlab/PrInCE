#' Calculate AIC for a Gaussian mixture model
#' 
#' Calculates the AIC for a curve fit with a Gaussian mixture model
#' by nonlinear least squares optimization. This function permits the 
#' calculation of the AIC after rejecting some Gaussians in the model, 
#' for example because their centres are outside the bounds of the profile.
#' 
#' @param coefs the coefficients of the Gaussian mixture model, output by 
#'   \code{\link{fit_gaussians}}
#' @param chromatogram the raw elution profile
#' 
#' @return the AIC of the fit model
gaussian_aic <- function(coefs, chromatogram) {
  # first, calculate log likelihood
  indices = which(!is.na(chromatogram))
  N = length(indices)
  fit = fit_curve(coefs, indices)
  res = chromatogram[indices] - fit
  w = rep_len(1, N)
  zw = w == 0
  loglik = -N * (log(2 * pi) + 1 - log(N) - sum(log(w + zw)) + 
                    log(sum(w * res^2)))/2
  # next, calculate parameters
  k = length(unlist(coefs)) + 1
  AIC = 2 * k - 2 * loglik
  return(AIC)
}
