#' Calculate the corrected AIC for a Gaussian mixture model
#' 
#' Calculates the corrected AIC for a curve fit with a Gaussian mixture model
#' by nonlinear least squares optimization.
#' 
#' @param fit the curve fit by the \code{stats::nls} function
#' 
#' @return the corrected AIC of the fit model
gaussian_aicc <- function(fit, indices) {
  # first, calculate log likelihood
  
  
  
  n <- length(fitted(fit))
  loglik <- as.numeric(logLik(fit))
  k <- length(coef(fit)) + 1
  AIC <- 2 * k - 2 * loglik
  AICc <- AIC + (2 * k * (k + 1)) / (n - k - 1)
  return(AICc)
}

# 
# res <- object$m$resid()
# N <- length(res)
# if (is.null(w <- object$weights)) 
#   w <- rep_len(1, N)
# zw <- w == 0
# val <- -N * (log(2 * pi) + 1 - log(N) - sum(log(w + zw)) + 
#                log(sum(w * res^2)))/2
# attr(val, "df") <- 1L + length(coef(object))
# attr(val, "nobs") <- attr(val, "nall") <- sum(!zw)
# class(val) <- "logLik"
# val