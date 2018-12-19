#' Calculate precision at each point in a sequence
#' 
#' Calculate the precision of a list of interactions at each point in the list,
#' given a set of labels.
#' 
#' @param labels a vector of zeroes (FPs) and ones (TPs)
#' 
#' @return a vector of the same length giving the precision at each point in the
#' input vector
#' 
#' @examples
#' ## calculate features
#' data(scott)
#' data(scott_gaussians)
#' subset = scott[seq_len(500), ] ## limit to first 500 proteins
#' gauss = scott_gaussians[names(scott_gaussians) %in% rownames(subset)]
#' features = calculate_features(subset, gauss)
#' ## make training labels
#' data(gold_standard)
#' ref = adjacency_matrix_from_list(gold_standard)
#' labels = make_labels(ref, features)
#' ## predict interactions with naive Bayes classifier
#' ppi = predict_NB_ensemble(features, labels, cv_folds = 3, models = 1)
#' ## tag precision of each interaction
#' ppi$precision = calculate_precision(ppi$label)
#' 
#' @export
calculate_precision <- function(labels) {
  if (is.factor(labels))
    labels <- as.integer(labels) - 1
  TPs <- cumsum(replace(labels, is.na(labels), 0))
  FPs <- cumsum(replace(labels == 0, is.na(labels), 0))
  precision <- TPs / (TPs + FPs)
}
