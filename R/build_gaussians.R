#' Deconvolve profiles into Gaussian mixture models
#' 
#' Identify peaks in co-fractionation profiles by deconvolving peaks in 
#' Gaussian mixture models. Models are mixtures of between 1 and 5 Gaussians.
#' Profiles are pre-processed prior to building Gaussians by filtering and 
#' cleaning. By default, profiles with fewer than 5 non-missing points, or 
#' fewer than 5 consecutive points after imputation of single missing values,
#' are removed. Profiles are cleaned by replacing missing values with 
#' near-zero noise, imputing single missing values as the mean of neighboring 
#' points, and smoothing with a moving average filter. 
#' 
#' @param profile_matrix a numeric matrix of co-elution profiles, with proteins
#' in rows
#' @param min_points filter profiles without at least this many total, 
#' non-missing points; passed to \code{\link{filter_profiles}}
#' @param min_consecutive filter profiles without at least this many 
#' consecutive, non-missing points; passed to \code{\link{filter_profiles}}
#' @param impute_NA if true, impute single missing values with the average of
#' neighboring values; passed to \code{\link{clean_profiles}}
#' @param smooth if true, smooth the chromatogram with a moving average filter;
#' passed to \code{\link{clean_profiles}}
#' @param smooth_width width of the moving average filter, in fractions;
#' passed to \code{\link{clean_profiles}}
#' @param max_gaussians the maximum number of Gaussians to fit; defaults to 5.
#' Note that Gaussian mixtures with more parameters than observed (i.e., 
#' non-zero or NA) points will not be fit. Passed to  
#' \code{\link{choose_gaussians}}
#' @param criterion the criterion to use for model selection;
#' one of "AICc" (corrected AIC, and default), "AIC", or "BIC". Passed to
#' \code{\link{choose_gaussians}}
#' @param max_iterations the number of times to try fitting the curve with
#' different initial conditions; defaults to 10. Passed to 
#' \code{\link{fit_gaussians}}
#' @param min_R_squared the minimum R-squared value to accept when fitting the
#' curve with different initial conditions; defaults to 0.5. Passed to 
#' \code{\link{fit_gaussians}}
#' @param method the method used to select the initial conditions for
#' nonlinear least squares optimization (one of "guess" or "random"); 
#' see \code{\link{make_initial_conditions}} for details. Passed to 
#' \code{\link{fit_gaussians}}
#' @param filter_gaussians_center true or false: filter Gaussians whose centres
#' fall outside the bounds of the chromatogram. Passed to 
#' \code{\link{fit_gaussians}}
#' @param filter_gaussians_height Gaussians whose heights are below this 
#' fraction of the chromatogram height will be filtered. Setting this value to
#' zero disables height-based filtering of fit Gaussians. Passed to 
#' \code{\link{fit_gaussians}}
#' @param filter_gaussians_variance_min Gaussians whose variance falls below 
#' this number of fractions will be filtered. Setting this value to
#' zero disables filtering. Passed to 
#' \code{\link{fit_gaussians}}
#' @param filter_gaussians_variance_max Gaussians whose variance is above
#' this number of fractions will be filtered. Setting this value to
#' zero disables filtering. Passed to 
#' \code{\link{fit_gaussians}}
#' 
#' @return a list of fit Gaussian mixture models, where each item in the 
#' list contains the following five fields: the number of Gaussians used to fit
#' the curve; the R^2 of the fit; the number of iterations used to 
#' fit the curve with different initial conditions; the coefficients of the 
#' fit model; and the curve predicted by the fit model. Profiles that 
#' could not be fit by a Gaussian mixture model above the minimum R-squared 
#' cutoff will be absent from the returned list.
#' 
#' @export
build_gaussians <- function(profile_matrix,
                            min_points = 1, min_consecutive = 5,
                            impute_NA = T, smooth = T, smooth_width = 4,
                            max_gaussians = 5, 
                            criterion = c("AICc", "AIC", "BIC"),
                            max_iterations = 10, min_R_squared = 0.5,
                            method = c("guess", "random"),
                            filter_gaussians_center = T,
                            filter_gaussians_height = 0.15,
                            filter_gaussians_variance_min = 0.1,
                            filter_gaussians_variance_max = 50) {
  # preprocess chromatograms: filter and clean
  filtered <- filter_profiles(profile_matrix,
                              min_points = min_points,
                              min_consecutive = min_consecutive)
  cleaned <- clean_profiles(filtered,
                            impute_NA = impute_NA,
                            smooth = smooth,
                            smooth_width = 4)
  
  # fit Gaussians, displaying progress bar
  gaussians <- list()
  proteins <- rownames(cleaned)
  P <- length(proteins)
  message(".. fitting Gaussian mixture models to ", P, " profiles")
  pb <- progress::progress_bar$new(
    format = "fitting :what [:bar] :percent eta: :eta",
    clear = F, total = P, width = 100)
  max_len <- max(nchar(proteins))
  for (i in seq_len(P)) {
    protein <- proteins[i]
    pb$tick(tokens = list(what = sprintf(paste0("%-", max_len, "s"), protein)))
    chromatogram <- cleaned[protein,]
    points <- sum(!is.na(profile_matrix[protein,]))
    gaussian <- choose_gaussians(chromatogram, points,
                                 max_gaussians, criterion,
                                 max_iterations, min_R_squared,
                                 method, filter_gaussians_center,
                                 filter_gaussians_height,
                                 filter_gaussians_variance_min,
                                 filter_gaussians_variance_max)
    gaussians[[protein]] <- gaussian
  }
  
  return(gaussians)
}
