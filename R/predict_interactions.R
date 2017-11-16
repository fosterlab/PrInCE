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
#' @param min_precision calculate precision for each interaction in a ranked
#' list until this value is reached 
#' 
#' @return a data frame containing the values of these five features for each
#' protein pair and the score output by the naive Bayes classifier
#' 
#' @export
predict_interactions <- function(profile_matrix, gaussians, gold_standard,
                                 min_precision = 0.5) {
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
  
  # train naive Bayes
  training_idxs <- which(!is.na(labels))
  training_labels <- as.factor(labels[training_idxs])
  training <- input[training_idxs, -c(1:2)]
  laplace <- mean(labels[training_idxs])
  nb <- naivebayes::naive_bayes(training, training_labels, laplace = laplace)
  
  # predictions <- naivebayes:::predict.naive_bayes(
  #   nb, training, type = 'prob')
  
  nb <- e1071::naiveBayes(training, training_labels)
  predictions <- e1071:::predict.naiveBayes(nb, training, type = 'raw')
  interactions <- cbind(input[training_idxs, 1:2], posterior = predictions[, "1"])
  interactions$label <- training_labels
  interactions <- dplyr::arrange(interactions, -posterior)
   
  
  # predict all data 
  predictions <- naivebayes:::predict.naive_bayes(
    nb, input[, -c(1:2)], type = 'prob')
  
  # create ranked data frame
  interactions <- cbind(input[, 1:2], posterior = predictions[, "1"])
  interactions <- dplyr::arrange(interactions, -posterior)
  
  # calculate precision
  interactions <- calculate_precision(interactions, labels, min_precision)
  
  return(interactions)
}
