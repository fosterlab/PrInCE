#' Calculate corrected AIC for a Gaussian mixture model
#' 
#' Calculates the corrected AIC for a curve fit with a Gaussian mixture model
#' by nonlinear least squares optimization. This function permits the 
#' calculation of the AICc after rejecting some Gaussians in the model, 
#' for example because their centres are outside the bounds of the profile.
#' 
#' @param coefs the coefficients of the Gaussian mixture model, output by 
#'   \code{\link{fit_gaussians}}
#' @param chromatogram the raw elution profile
#' 
#' @return the corrected AIC of the fit model
gaussian_aicc <- function(coefs, chromatogram) {
  # first, calculate AIC
  AIC = gaussian_aic(coefs, chromatogram)
  # second, calculate AICc
  N = sum(!is.na(chromatogram))
  k = length(unlist(coefs)) + 1
  AICc = AIC + (2 * k * (k + 1)) / (N - k - 1)
  return(AICc)
}
