#' Calculate the co-apex score for every protein pair
#' 
#' Calculate the co-apex score for every pair of proteins. This is defined as 
#' the minimum Euclidean distance between any two Gaussians fit to each 
#' profile. 
#' 
#' @param gaussians a list of Gaussian mixture models fit to the profile matrix
#' by \code{link{build_gaussians}}
#' @param proteins all proteins being scored, optionally including those 
#' without Gaussian fits 
#' 
#' @return a matrix of co-apex scores 
#' 
#' @examples 
#' data(scott_gaussians)
#' gauss <- scott_gaussians[seq_len(25)]
#' CA <- co_apex(gauss)
#' 
#' @importFrom stats dist
#' @importFrom purrr map
#' 
#' @export
co_apex <- function(gaussians, proteins = NULL) {
  if (is.null(proteins)) {
    proteins <- names(gaussians)
  }
  n_proteins <- length(proteins)
  n_gaussians <- length(gaussians)
  gaussian_names <- names(gaussians)
  ## first, calculate Euclidean distance between every pair of Gaussians 
  gaussian_centers <- map(gaussians, c("coefs", "mu"))
  gaussian_sigmas <- map(gaussians, c("coefs", "sigma"))
  gaussian_matrix <- cbind(unlist(gaussian_centers), unlist(gaussian_sigmas))
  CA <- as.matrix(dist(gaussian_matrix))
  ## get indices for each protein 
  gaussian_indices <- rep(gaussian_names, lengths(gaussian_centers))
  ## calculate min co-apex score 
  co_apex <- matrix(NA, nrow = n_proteins, ncol = n_proteins,
                    dimnames = list(proteins, proteins))
  gaussian_names <- intersect(names(gaussians), proteins)
  for (protein_A in gaussian_names) {
    idxs_A <- which(gaussian_indices == protein_A)
    for (protein_B in gaussian_names) {
      idxs_B <- which(gaussian_indices == protein_B)
      co_apex[protein_A, protein_B] <- min(CA[idxs_A, idxs_B]) 
    }
  }
  return(co_apex)
}