#' Fit a mixture of Gaussians to a chromatogram curve
#'
#' Fit mixtures of one or more Gaussians to the curve formed by a chromatogram
#' profile, using nonlinear least-squares. 
#' 
#' @param chromatogram a numeric vector corresponding to the chromatogram trace
#' @param n_gaussians the number of Gaussians to fit
#' @param indices optionally, specify the indices of the fractions for each
#'   point in the chromatogram; by default, these are defined by the base 
#'   R function \code{\link[base]{seq_along}}
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
#' @param filter_gaussians_variance_min Gaussians whose variance falls below 
#' this number of fractions will be filtered. Setting this value to
#' zero disables filtering.
#' @param filter_gaussians_variance_max Gaussians whose variance is above
#' this number of fractions will be filtered. Setting this value to
#' zero disables filtering.
#'
#' @return a list with seven entries: the number of Gaussians used to fit
#' the curve; the R^2 of the fit; the number of iterations used to 
#' fit the curve with different initial conditions; the coefficients of the 
#' fit model; the log-likelihood of the fit model; and the fited curve predicted
#' by the model.
#' 
#' @examples
#' data(scott)
#' chrom = clean_profile(scott[1, ])
#' fit = fit_gaussians(chrom, n_gaussians = 1)
#' 
#' @importFrom stats coef cor setNames
#' @importFrom minpack.lm nlsLM
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' 
#' @export
fit_gaussians <- function(chromatogram, 
                          n_gaussians,
                          indices = NULL,
                          max_iterations = 50, 
                          min_R_squared = 0.5,
                          method = c("guess", "random"),
                          filter_gaussians_center = TRUE,
                          filter_gaussians_height = 0.01,
                          filter_gaussians_variance_min = 0.1,
                          filter_gaussians_variance_max = 50) {
  if (is.null(indices)) {
    indices <- seq_along(chromatogram)
  }
  iter = 0
  bestR2 = -Inf
  bestLogLik = -Inf
  bestCoefs = NULL
  while (iter < max_iterations & bestR2 < min_R_squared) {
    # increment iteration counter
    iter <- iter + 1
    
    # create formula
    ## (p1 = A, p2 = mu, p3 = sigma)
    get_equation = function(i, var_char = 'p') {
      ## general form: A[i] * exp(-((x - mu[i])/sigma[i])^2)
      paste0(var_char, (3 * i - 2), " * exp(-((index - ", var_char, 
             (3 * i - 1), ")/", var_char, (3 * i), ")^2)")
    }
    # get the gaussian mixture part of the equation
    gaussian_eqn = paste(get_equation(seq_len(n_gaussians)), collapse = " + ")
    # add LHS to construct formula
    model = as.formula(paste("abundance ~ ", gaussian_eqn))
    
    # make initial conditions
    init = make_initial_conditions(chromatogram, 
                                   n_gaussians = n_gaussians,
                                   indices = indices)
    start = unlist(map(seq_len(n_gaussians), ~ c(
      init$A[.], init$mu[.], init$sigma[.]))) %>%
      setNames(paste0("p", seq_len(n_gaussians * 3)))
    
    # fit the model
    df = data.frame(abundance = chromatogram,
                    index = indices)
    fit <- tryCatch({
      suppressWarnings(
        minpack.lm::nlsLM(formula = model, 
                          data = df, 
                          start = start,
                          trace = FALSE))
    }, error = function(e) { 
      e 
    }, simpleError = function(e) { 
      e
    })
    if ("error" %in% class(fit))
      next
    
    # split up fit coefficients vector into list with three entries
    coefs <- coef(fit)
    coefs <- split(coefs, rep(seq_len(3), each = n_gaussians))
    coefs <- setNames(coefs, c("A", "mu", "sigma"))
    
    if (filter_gaussians_variance_min > 0) {
      # remove Gaussians with negative variances 
      sigmas <- coefs[["sigma"]]
      drop <- which(sigmas < filter_gaussians_variance_min)
      if (length(drop) > 0)
        coefs <- lapply(coefs, `[`, -drop)
    }
    
    if (filter_gaussians_variance_max > 0) {
      # remove Gaussians with extremely large amplitude
      sigmas <- coefs[["sigma"]]
      drop <- which(sigmas > filter_gaussians_variance_max)
      if (length(drop) > 0)
        coefs <- lapply(coefs, `[`, -drop)
    }
    
    if (filter_gaussians_center) {
      # remove Gaussians outside bounds of chromatogram
      means <- coefs[["mu"]]
      drop <- which(means < 0 | means > length(chromatogram))
      if (length(drop) > 0)
        coefs <- lapply(coefs, `[`, -drop)
    }
    
    if (filter_gaussians_height > 0) {
      # remove Gaussians less than 15% of height
      minHeight <- max(chromatogram) * filter_gaussians_height
      heights <- coefs[["A"]]
      drop <- which(heights < minHeight)
      if (length(drop) > 0)
        coefs <- lapply(coefs, `[`, -drop)
    }
    
    # calculate R^2
    if (length(dplyr::first(coefs)) == 0)
      next
    curveFit <- fit_curve(coefs, indices)
    R2 <- cor(chromatogram, curveFit)^2
    # replace best fit with this model?
    if (R2 > bestR2 & R2 > min_R_squared) {
      bestR2 <- R2
      bestLogLik = logLik(fit)
      bestCoefs <- coefs
    }
  }
  if (!is.null(bestCoefs)) {
    curveFit <- fit_curve(bestCoefs, indices)
  } else {
    curveFit <- NULL
  }
  results <- list(n_gaussians = n_gaussians,
                  R2 = bestR2, 
                  iterations = iter, 
                  coefs = bestCoefs, 
                  logLik = bestLogLik,
                  curveFit = curveFit)
  return(results)
}
