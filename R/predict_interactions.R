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
  # calculate features
  input <- calculate_features(profile_matrix, gaussians)
  
  # identify true positives
  labels <- make_labels(gold_standard, input)

  # predict with an ensemble of naive Bayes classifiers
  interactions <- predict_NB_ensemble(input, labels)
  
  # calculate precision
  interactions$precision <- calculate_precision(interactions$label)

  return(interactions)
}
