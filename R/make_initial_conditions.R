#' Make initial conditions for curve fitting with a mixture of Gaussians
#' 
#' Construct a set of initial conditions for curve fitting using nonlinear
#' least squares using a mixture of Gaussians. The "guess" method ports 
#' code from the Matlab release of PrInCE. This method finds local maxima 
#' within the chromatogram, orders them by their separation (in number of 
#' fractions) from the previous local maxima, and uses the positions and heights
#' of these local maxima (+/- some random noise) as initial conditions for 
#' Gaussian curve-fitting. The "random" method simply picks random values 
#' within the fraction and intensity intervals as starting points for Gaussian
#' curve-fitting. The initial value of sigma is set by default to a random 
#' number within +/- 0.5 of two for both modes; this is based on our manual 
#' inspection of a large number of chromatograms. 
#' 
#' @param chromatogram a numeric vector corresponding to the chromatogram trace
#' @param n_gaussians the number of Gaussians being fit 
#' @param clean if \code{TRUE}, use the cleaned version of the profile to guess
#' initial conditions
#' @param method one of "guess" or "random", discussed above
#' @param sigma_default the default mean initial value of sigma
#' @param sigma_noise the amount of random noise to add or subtract from 
#' the default mean initial value of sigma
#' @param mu_noise the amount of random noise to add or subtract from the 
#' Gaussian centers in "guess" mode
#' @param A_noise the amount of random noise to add or subtract from the 
#' Gaussian heights in "guess" mode 
#' 
#' @return a list of three numeric vectors (A, mu, and sigma), each having
#' a length equal to the maximum number of Gaussians to fit 
#' 
#' @examples
#' data(scott)
#' set.seed(0)
#' start = make_initial_conditions(chrom, n_gaussians = 2, method = "guess")
#' 
#' @importFrom purrr map_dbl map_lgl
#' 
#' @export
make_initial_conditions <- function(chromatogram, n_gaussians, 
                                    clean = TRUE, 
                                    indices = NULL,
                                    method = c("guess", "random"),
                                    sigma_default = 2, 
                                    sigma_noise = 0.5, mu_noise = 1.5, 
                                    A_noise = 0.5) {
  method <- match.arg(method)
  if (is.null(indices)) {
    indices = seq_along(chromatogram)
  }
  if (clean) {
    ## identify points (if any) where chromatogram should be split
    pos = which(map_lgl(seq_len(length(indices) - 1), ~
                          indices[. + 1] < indices[.])) + 1
    ## clean profiles and reassemble
    split = unname(split(chromatogram, cumsum(
      seq_along(chromatogram) %in% pos)))
    chromatogram = unlist(map(split, clean_profile))
  }
  if (method == "guess") {
    # find fractions that represent local maxima
    peaksX = indices[which(diff(sign(diff(chromatogram))) == -2) + 1]
    peaksX = peaksX[is.finite(peaksX)]
    # catch local minima at start or end
    peaksX = c(peaksX, 1)
    peaksX = c(peaksX, max(indices))
    # drop redundant peaks
    peaksX = unique(peaksX)
    # drop peaks with missing values
    peaksY = suppressWarnings(
      map_dbl(peaksX, ~ max(chromatogram[indices == .], na.rm = T)))
    peaksX = peaksX[is.finite(peaksY)]
    # if there is still nothing, just pick the maximum value
    if (length(peaksX) == 0) {
      peaksX = which.max(chromatogram)
    }
    # get intensities for sorted local maxima
    ## (picking max value at each idx if multiple exist)
    peaksY = suppressWarnings(
      map_dbl(peaksX, ~ max(chromatogram[indices == .], na.rm = T)))
    # order by height
    peaksX <- peaksX[order(-peaksY)]
    peaksY <- peaksY[order(-peaksY)]
    # use local maxima as initial conditions for Gaussian fitting
    A <- numeric(0)
    mu <- numeric(0)
    sigma <- numeric(0)
    for (i in seq_len(n_gaussians)) {
      A[i] <- abs(peaksY[i] + runif(1) - 0.5)
      mu[i] <- peaksX[i] + runif(1, max = 3) - 1.5
      sigma[i] <- sigma_default + runif(1) - 0.5
    }
    # if there are not enough peaks, pick additional random values
    fractions <- length(chromatogram)
    A[is.na(A)] <- runif(sum(is.na(A)), max = fractions)
    mu[is.na(mu)] <- runif(sum(is.na(mu)), max = max(peaksY))
    initial_conditions <- list(A = A, mu = mu, sigma = sigma)
    return(initial_conditions)
  } else if (method == "random") {
    minHeight <- min(chromatogram)
    maxHeight <- max(chromatogram)
    A <- runif(n_gaussians, min = minHeight, max = maxHeight)
    mu <- runif(n_gaussians, min = 1, max = length(chromatogram))
    sigma <- sigma_default + runif(n_gaussians) - 0.5
    initial_conditions <- list(A = A, mu = mu, sigma = sigma)
    return(initial_conditions)
  }
}
