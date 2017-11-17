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
#' @return a ranked data frame of pairwise interactions, with the 
#' classifier score, label, and cumulative precision for each interaction 
#' 
#' @export
predict_interactions <- function(profile_matrix, gaussians, gold_standard) {
  # replace missing values with near-zero noise
  cleaned <- clean_profiles(profile_matrix, impute_NA = F, smooth = F,
                            noise_floor = 0.05)
  proteins <- rownames(cleaned)
  n_proteins <- length(proteins)
  
  # calculate Pearson correlation and P-value distances
  cor_R_raw <- 1 - cor(t(profile_matrix), use = 'pairwise.complete.obs')
  cor_R_cleaned <- 1 - cor(t(cleaned))
  cor_P <- Hmisc::rcorr(t(profile_matrix))$P
  ## set P-values with 2 pairwise observations to zero 
  pairs <- crossprod(t(!is.na(profile_matrix)))
  cor_P[pairs <= 2] <- NA
  # calculate Euclidean distance
  eucl <- as.matrix(dist(cleaned, method = 'euclidean'))
  # calculate co-peak distance
  maxes <- apply(cleaned, 1, which.max)
  co_peak <- as.matrix(dist(maxes))
  # calculate co-apex (Gaussian) score
  CA <- co_apex(gaussians, proteins)
  
  # collapse features
  feature_matrices <- list(cor_R_raw, cor_R_cleaned, cor_P, eucl, co_peak, CA)
  ## make sure all dimensions are identical
  if (!all(purrr::map_int(feature_matrices, nrow) == n_proteins) |
      !all(purrr::map_int(feature_matrices, ncol) == n_proteins))
    stop("at least one feature matrix did not have correct dimensions")
  tri <- upper.tri(co_peak)
  idxs <- which(tri, arr.ind = T)
  input <- data.frame(protein_A = rownames(co_peak)[idxs[, 1]], 
                      protein_B = rownames(co_peak)[idxs[, 2]]) 
  input <- cbind(input, purrr::map(feature_matrices, ~ .[tri]))
  colnames(input)[3:8] <- c("cor_R_raw", "cor_R_cleaned", "cor_P", 
                            "euclidean_distance", "co_peak", "co_apex")

  # identify true positives
  label_mat <- make_label_matrix(gold_standard, co_peak)
  tri <- upper.tri(label_mat)
  labels <- label_mat[tri]
  
  # predict with an ensemble of naive Bayes classifiers
  predictions <- predict_NB_ensemble(input, labels)
  
  # create ranked data frame
  interactions <- cbind(input[, 1:2], score = predictions, label = labels)
  interactions <- dplyr::arrange(interactions, -score)
  
  # calculate precision
  interactions$precision <- calculate_precision(interactions$label)

  return(interactions)
}
