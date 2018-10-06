#' Identify significantly changing proteins using a likelihood ratio test
#' 
#' This function identifies proteins with a significant change between two
#' conditions in a protein correlation profiling experiment, by fitting two 
#' different mixtures of Gaussians: one with the condition as a covariate, and
#' one without. The two models are then compared with a likelihood ratio test
#' to evaluate whether taking condition into account provides a significantly 
#' better fit.
#' 
#' @param chrom1 the chromatogram of the protein of interest in the first 
#'    condition (e.g., the M/L channel)
#' @param chrom2 the chromatogram of the protein of interest in the first 
#'    condition (e.g., the H/L channel)
#'
#' @export
fold_change = function(chrom1, chrom2) {
  # length must be equal
  if (length(chrom1) != length(chrom2))
    stop("chromatogram length must be equal")
  
  # create data frame
  df = data.frame(
    condition = c(rep('1', length(chrom1)), rep('2', length(chrom2))),
    fraction = c(seq_along(chrom1), seq_along(chrom2)),
    ratio = c(chrom1, chrom2))
  
  # fit model 1 (without condition)
  ## try for different numbers of Gaussians
  for (n_gauss in seq(min_gaussians, max_gaussians)) {
    
  }
  
  # fit model 2 (with condition)
  
  # likelihood ratio test
  # beta = coef(full)["severity"],
  # loglik_full = logLik(full),
  # loglik_red = logLik(reduced),
  # lrt_stat = 2 * (loglik_full - loglik_red),
  # lrt_pval = pchisq(abs(lrt_stat), df = 1, lower.tail = F)) 
}

#' Fit the mixture of Gaussians used by the \code{fold_change} function
#' 
#' 
#' @param df a data frame containing the variables \code{fraction}, 
#'   \code{ratio}, and \code{condition}
#' @param fit_condition whether to include condition as a covariate in the model
#' @param min_gaussians the inimum number of Gaussians to fit; defaults to 1.
#' @param max_gaussians the maximum number of Gaussians to fit; defaults to 5.
#'   Gaussian mixtures with more parameters than observed (i.e., 
#'   non-zero or NA) points will not be fit. 
#' @param max_iterations the number of times to try fitting the curve with
#' different initial conditions; defaults to 10
#' @param criterion the criterion to use for model selection;
#'  one of \code{AICc} (corrected AIC, and default), \code{AIC}, or \code{BIC}
#' @param min_R_squared the minimum R-squared value to accept when fitting the
#'   curve with different initial conditions; defaults to 0.5
#' @param filter_gaussians_center true or false: filter Gaussians whose centres
#'   fall outside the bounds of the chromatogram 
#' @param filter_gaussians_height Gaussians whose heights are below this 
#'   fraction of the chromatogram height will be filtered. Setting this value to
#'   zero disables height-based filtering of fit Gaussians
#' @param filter_gaussians_variance_min Gaussians whose variance falls below 
#'   this number of fractions will be filtered. Setting this value to
#'   zero disables filtering.
#' @param filter_gaussians_variance_max Gaussians whose variance is above
#'   this number of fractions will be filtered. Setting this value to
#'   zero disables filtering.
#' @param sigma_default the default mean initial value of sigma
#'   
#' @importFrom tidyr %>% drop_na
#' @importFrom dplyr group_by summarise pull mutate
#' @importFrom purrr map
#' @importFrom nlme gnls
#' @importFrom magrittr %<>% 
fit_model = function(df, 
                     fit_condition = F, 
                     min_gaussians = 1,
                     max_gaussians = 5, 
                     max_iterations = 10,
                     criterion = c("AICc", "AIC", "BIC"),
                     min_R_squared = 0.5,
                     filter_gaussians_center = T,
                     filter_gaussians_height = 0.15,
                     filter_gaussians_variance_min = 0.1,
                     filter_gaussians_variance_max = 50,
                     sigma_default = 2) {
  criterion <- match.arg(criterion)
  
  # don't fit mixtures with more parameters than (experimental) points
  points = df %>%
    dplyr::group_by(condition) %>%
    dplyr::summarise(n = sum(is.finite(ratio))) %>%
    dplyr::pull(n)
  min_points = min(points)
  max_gaussians = min(max_gaussians, floor(min_points / 3))

  # fix any issues with input data 
  df %<>% 
    tidyr::drop_na(ratio) %>%
    dplyr::mutate(condition = factor(condition))
  
  # fit models with different numbers of Gaussians
  fits = list()
  for (n_gaussians in seq(min_gaussians, max_gaussians)) {
    iter = 0
    bestR2 = 0
    while (iter < max_iterations & bestR2 < min_R_squared) {
      # create formula
      ## (x1 = A, x2 = mu, x3 = sigma)
      get_equation = function(i, var_char = 'p') {
        paste0(var_char, (3 * i - 2), " * exp(-((fraction - ", var_char, 
               (3 * i - 1), "/", var_char, (3 * i), ")^2))")
      }
      # get the gaussian mixture part of the equation
      gaussian_eqn = paste(get_equation(seq_len(n_gaussians), 'p'), 
                           collapse = " + ")
      # add LHS and (optionally) condition
      if (fit_condition) {
        model = as.formula(paste("ratio ~ condition +", gaussian_eqn))
      } else {
        model = as.formula(paste("ratio ~ ", gaussian_eqn))
      }
      
      # make initial conditions
      mean = df %>%
        dplyr::group_by(fraction) %>%
        dplyr::summarise(ratio = mean(ratio, na.rm = T)) %>%
        dplyr::pull(ratio)
      init = PrInCE:::make_initial_conditions(mean, n_gaussians = n_gaussians)
      start = unlist(purrr::map(seq_len(n_gaussians), ~ c(
        init$A[.], init$mu[.], init$sigma[.]))) %>% 
        setNames(paste0("p", seq_len(n_gaussians * 3)))
      
      # fit the model
      fit = nlme::gnls(model = model, data = df, start = start,
                       control = nlme::gnlsControl(nlsTol = 0.1))
    }
  }
  
  # choose best model
  
}


