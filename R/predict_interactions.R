#' Predict interactions for a profile matrix
#'
#' Predict pairwise interactions between proteins in a profile matrix, using a
#' naive Bayes classifier trained on dataset-derived features including:
#' * one minus the Pearson correlation coefficient and 
#' * the corresponding P value, 
#' * the Euclidean distance between profiles, 
#' * the difference (in fractions) between the maximum values of each profile, 
#' (the peak score) and
#' * the co-apex score, defined as the minimum Euclidean distance between any
#' two Gaussians
#' 
#' @param profile_matrix a numeric matrix of co-elution profiles, with proteins
#' in rows
#' @param gaussians a list of Gaussian mixture models fit to the profile matrix
#' by \code{link{build_gaussians}}
#' @param gold_standard a matrix or data frame of "gold standard" interactions
#' used to train a naive Bayes classifier
#' 
#' @return a data frame containing the values of these five features for each
#' protein pair and the score output by the naive Bayes classifier
#' 
#' @export
predict_interactions <- function(profile_matrix, gaussians, gold_standard) {
  # get cleaned chromatograms
  cleaned <- clean_profiles(profile_matrix)
  proteins <- rownames(cleaned)
  n_proteins <- length(proteins)
  
  # calculate Pearson correlation and P-value distances
  cors <- Hmisc::rcorr(t(cleaned))
  cor.R <- 1 - cors$r
  cor.P <- cors$P
  # calculate Euclidean distance
  eucl <- as.matrix(dist(cleaned, method = 'euclidean'))
  # calculate co-peak distance
  maxes <- apply(cleaned, 1, which.max)
  co_peak <- as.matrix(dist(maxes))
  # calculate co-apex (Gaussian) score
  CA <- co_apex(gaussians, proteins)
  
  # collapse features
  feature_matrices <- list(cor.R, cor.P, eucl, co_peak, co_apex)
  ## make sure all dimensions are identical
  if (!all(purrr::map_int(feature_matrices, nrow) == n_proteins) |
      !all(purrr::map_int(feature_matrices, ncol) == n_proteins))
    stop("at least one feature matrix did not have correct dimensions")
  tri <- upper.tri(co_peak)
  idxs <- which(tri, arr.ind = T)
  input <- data.frame(protein_A = rownames(co_peak)[idxs[, 1]], 
                   protein_B = rownames(co_peak)[idxs[, 2]]) 
  input <- cbind(input, purrr::map(feature_matrices, ~ .[tri]))
  colnames(input)[3:7] <- c("cor_R", "cor_P", "euclidean_distance",
                            "co_peak", "co_apex")

  # identify true positives
  
  # run naive Bayes 
}
