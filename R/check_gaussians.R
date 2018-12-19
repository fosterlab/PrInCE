#' Check the format of a list of Gaussians
#' 
#' Test whether an input list of Gaussians conforms to the format expected
#' by PrInCE: that is, a named list with five fields for each entry, i.e., the
#' number of Gaussians in the mixture model, the {\eqn{r^2}} value, the 
#' number of iterations used by \code{\link[stats]{nls}}, the coefficients of
#' each model, and the fitted curve. 
#' 
#' Optionally, some extra checks will be 
#' done on the fraction of proteins in the complete dataset for which a 
#' Gaussian mixture model could be fit, if provided. In particular, the function
#' will throw an error if fewer than \code{n_error} proteins have a fitted 
#' Gaussian, and emit a warning if fewer than \code{pct_warning} do. 
#'  
#' @param gaussians the list of Gaussians
#' @param proteins the complete set of input proteins
#' @param replicate_idx the replicate being analyzed, if input proteins are 
#'   provided; used to throw more informative error messages
#' @param n_error minimum number of proteins that can have fitted Gaussians
#'   without throwing an error
#' @param pct_warning minimum fraction of proteins that can have fitted
#'   Gaussians without giving a warning
#' 
#' @return \code{TRUE} if all conditions are met, but throws an error if any 
#'   is not
#'   
#' @export
check_gaussians = function(gaussians, proteins = NULL, replicate_idx = NULL,
                           n_error = 3, pct_warning = 0.1) {
  if (!is.null(gaussians)) {
    if (!is.list(gaussians)) {
      stop("`gaussians` object is not a list")
    }
    if (is.null(names(gaussians))) {
      stop("`gaussians` object is not a named list")
    }
    if (any(!purrr::map_lgl(gaussians, ~ "n_gaussians" %in% names(.)))) {
      stop("not all gaussians have the `n_gaussians` field")
    }
    if (any(!purrr::map_lgl(gaussians, ~ "R2" %in% names(.)))) {
      stop("not all gaussians have the `R2` field")
    }
    if (any(!purrr::map_lgl(gaussians, ~ "iterations" %in% names(.)))) {
      stop("not all gaussians have the `iterations` field")
    }
    if (any(!purrr::map_lgl(gaussians, ~ "coefs" %in% names(.)))) {
      stop("not all gaussians have the `coefs` field")
    }
    if (any(!purrr::map_lgl(gaussians, ~ "curveFit" %in% names(.)))) {
      stop("not all gaussians have the `curveFit` field")
    }
  } else {
    stop("`gaussians` object is null")
  }
  if (!is.null(proteins)) {
    repl_string = ifelse(!is.null(replicate_idx),
                         paste(" in replicate ", replicate_idx), "")
    # error if too few proteins have Gaussians
    n_gaussians = sum(proteins %in% names(gaussians))
    if (n_gaussians <= n_error) {
      stop("only ", n_gaussians, " proteins have fitted Gaussians", repl_string, 
           " (minimum: ", n_error, ")")
    }
    # warn if a low % of proteins have Gaussians
    pct_gaussians = mean(proteins %in% names(gaussians))
    if (pct_gaussians < pct_warning) {
      warning("only ", format(100 * pct_gaussians, digits = 2), 
              " of proteins have fitted Gaussians", repl_string, " (minimum: ", 
              format(100 * pct_warning, digits = 2), ")")
    }
  }
  return(T)
}
