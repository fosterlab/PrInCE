#' Output the fit curve for a given mixture of Gaussians
#' 
#' For a Gaussian mixture model fit to a curve by \code{\link{fit_gaussians}},
#' output the fit curve using the coefficients rather than the \code{nls} 
#' object. This allows individual Gaussians to be removed from the fit model:
#' for example, if their height is below a certain threshold, or their 
#' centres are outside the bounds of the chromatogram. 
#' 
#' @param coef numeric vector of coefficients for a Gaussian mixture model fit 
#' by \code{\link{fit_gaussians}}. This function assumes that the heights of the 
#' Gaussians are specified by coefficients beginning with "A" 
#' ("A1", "A2", "A3", etc.), centres are specified by coefficients beginning
#' with "mu", and standard deviations are specified by coefficients beginning
#' with "sigma".  
#' @param indices the indices, or x-values, to predict a fitted curve for 
#' (for example, the fractions in a given chromatogram)
#' 
#' @return the fitted curve
#' 
#' @export
fit_curve <- function(coef, indices) {
  A <- coef[["A"]]
  mu <- coef[["mu"]]
  sigma <- coef[["sigma"]]
  gaussians <- length(A)
  rowSums(sapply(seq_len(gaussians), function(i) 
    A[i] * exp(-((indices - mu[i])/sigma[i])^2)))
}
