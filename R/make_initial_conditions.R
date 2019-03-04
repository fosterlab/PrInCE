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
#' chrom <- clean_profile(scott[16, ])
#' set.seed(0)
#' start <- make_initial_conditions(chrom, n_gaussians = 2, method = "guess")
#' 
#' @importFrom dplyr first last nth
#' 
#' @export
make_initial_conditions <- function(chromatogram, n_gaussians, 
                                    method = c("guess", "random"),
                                    sigma_default = 2, 
                                    sigma_noise = 0.5, mu_noise = 1.5, 
                                    A_noise = 0.5) {
  method <- match.arg(method)
  if (method == "guess") {
    # find fractions that represent local maxima
    peaksX <- which(diff(sign(diff(chromatogram))) == -2) + 1
    # catch local minima at start or end
    if (first(chromatogram) > nth(chromatogram, 2))
      peaksX <- c(peaksX, 1)
    if (last(chromatogram) > nth(chromatogram, -2))
      peaksX <- c(peaksX, length(chromatogram))
    # order local maxima by distance (on one side, only)
    ## distances <- diff(peaksX)
    ## peaksX <- peaksX[order(-distances)]
    # get intensities for sorted local maxima
    peaksY <- chromatogram[peaksX]
    # order by height
    peaksX <- peaksX[order(-peaksY)]
    # use local maxima as initial conditions for Gaussian fitting
    A <- numeric(0)
    mu <- numeric(0)
    sigma <- numeric(0)
    for (i in seq_len(n_gaussians)) {
      A[i] <- peaksY[i] + runif(1) - 0.5
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
