#' Fit a Gaussian mixture model to a co-elution profile
#'
#' Fit mixtures of one or more Gaussians to the curve formed by a chromatogram
#' profile, and choose the best fitting model using an information criterion
#' of choice. 
#' 
#' @param chromatogram a numeric vector corresponding to the chromatogram trace
#' @param points optional, the number of non-NA points in the raw data
#' @param max_gaussians the maximum number of Gaussians to fit; defaults to 5.
#' Note that Gaussian mixtures with more parameters than observed (i.e., 
#' non-zero or NA) points will not be fit. 
#' @param criterion the criterion to use for model selection;
#'  one of "AICc" (corrected AIC, and default), "AIC", or "BIC"
#' @param max_iterations the number of times to try fitting the curve with
#' different initial conditions; defaults to 10
#' @param min_R_squared the minimum R-squared value to accept when fitting the
#' curve with different initial conditions; defaults to 0.5
#' @param method the method used to select the initial conditions for
#' nonlinear least squares optimization (one of "guess" or "random"); 
#' see \code{\link{make_initial_conditions}} for details
#' @param filter_gaussians_center true or false: filter Gaussians whose centres
#' fall outside the bounds of the chromatogram 
#' @param filter_gaussians_height Gaussians whose heights are below this 
#' fraction of the chromatogram height will be filtered. Setting this value to
#' zero disables height-based filtering of fit Gaussians
#' @param filter_gaussians_variance_min Gaussians whose variance is below this 
#' threshold will be filtered. Setting this value to zero disables filtering.
#' @param filter_gaussians_variance_max Gaussians whose variance is above this 
#' threshold will be filtered. Setting this value to zero disables filtering.
#'
#' @return a list with five entries: the number of Gaussians used to fit
#' the curve; the R^2 of the fit; the number of iterations used to 
#' fit the curve with different initial conditions; the coefficients of the 
#' fit model; and the curve predicted by the fit model.
#' 
#' @examples
#' data(scott)
#' chrom = clean_profile(scott[1, ])
#' gauss = choose_gaussians(chrom, max_gaussians = 3)
#'
#' @export
choose_gaussians <- function(chromatogram, 
                             indices = NULL,
                             points = NULL,
                             max_gaussians = 5, 
                             criterion = c("AICc", "AIC", "BIC"),
                             max_iterations = 50, min_R_squared = 0.5,
                             method = c("guess", "random"),
                             filter_gaussians_center = TRUE,
                             filter_gaussians_height = 0.15,
                             filter_gaussians_variance_min = 0.1,
                             filter_gaussians_variance_max = 50) {
  criterion <- match.arg(criterion)
  
  # don't fit mixtures with more parameters than (experimental) points
  if (!is.null(points)) {
    max_gaussians <- min(max_gaussians, floor(points / 3))
  }
  
  # fit and choose models 
  fits <- list()
  for (n_gaussians in seq_len(max_gaussians))
    fits[[n_gaussians]] <- fit_gaussians(
      chromatogram, 
      n_gaussians = n_gaussians,
      indices = indices, 
      max_iterations = max_iterations, 
      min_R_squared = min_R_squared,
      method = method, 
      filter_gaussians_center = filter_gaussians_center, 
      filter_gaussians_height = filter_gaussians_height, 
      filter_gaussians_variance_min = filter_gaussians_variance_min,
      filter_gaussians_variance_max = filter_gaussians_variance_max)
  
  # remove any models that failed to fit
  models <- purrr::map(fits, "coefs")
  drop <- purrr::map_lgl(models, is.null)
  fits <- fits[!drop]
  coefs <- purrr::map(fits, "coefs")
  
  if (criterion == "AICc") {
    criteria <- lapply(coefs, gaussian_aicc, chromatogram)
  } else if (criterion == "AIC") { 
    criteria <- lapply(coefs, gaussian_aic)
  } else if (criterion == "BIC") {
    criteria <- lapply(coefs, gaussian_bic)  
  }
  best <- which.min(criteria)
  if (length(best) == 0) {
    return(NULL)
  } else {
    return(fits[[best]])
  }
} 
