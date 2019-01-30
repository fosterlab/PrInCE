#' Predict interactions given a set of features and examples
#' 
#' Discriminate interacting from non-interacting protein pairs by training a
#' machine learning model on a set of labelled examples, given a set of 
#' features derived from a co-elution profile matrix (see 
#' \code{\link[PrInCE]{calculate_features}}. 
#' 
#' PrInCE implements four different classifiers (naive Bayes, support vector
#' machine, random forest, and logistic regression). Naive Bayes is used as a
#' default. The classifiers are trained
#' on the gold standards using a ten-fold cross-validation procedure, training
#' on 90% of the data and predicting on the remaining 10%. For protein pairs
#' that are part of the training data, the held-out split is used to assign
#' a classifier score, whereas for the remaining protein pairs, the median of
#' all ten folds is used. Furthermore, to ensure the results are not sensitive
#' to the precise classifier split used, an ensemble of multiple classifiers
#' (ten, by default) is trained, and the classifier score is subsequently
#' averaged across classifiers. 
#' 
#' PrInCE can also ensemble across multiple different types of classifiers, 
#' by supplying the \code{"ensemble"} option to the \code{classifier} argument. 
#' 
#' @param features a data frame with proteins in the first two columns, and 
#'   features to be passed to the classifier in the remaining columns
#' @param gold_standard an adjacency matrix of "gold standard" interactions
#'  used to train the classifier 
#' @param classifier the type of classifier to use: one of \code{"NB"} (naive
#'   Bayes), \code{"SVM"} (support vector machine), \code{"RF"} (random forest),
#'   \code{"LR"} (logistic regression), or \code{"ensemble"} (an ensemble of
#'   all four)
#' @param verbose if \code{TRUE}, print a series of messages about the stage
#'   of the analysis
#' @param models the number of classifiers to train and average across, each
#'   with a different k-fold cross-validation split
#' @param cv_folds the number of folds to use for k-fold cross-validation
#' @param trees for random forests only, the number of trees in the forest
#' 
#' @return a ranked data frame of pairwise interactions, with the 
#' classifier score, label, and cumulative precision for each interaction 
#' 
#' @examples 
#' ## calculate features
#' data(scott)
#' data(scott_gaussians)
#' subset <- scott[seq_len(500), ] ## limit to first 500 proteins
#' gauss <- scott_gaussians[names(scott_gaussians) %in% rownames(subset)]
#' features <- calculate_features(subset, gauss)
#' ## load training data
#' data(gold_standard)
#' ref <- adjacency_matrix_from_list(gold_standard)
#' ## predict interactions
#' ppi <- predict_interactions(features, ref, cv_folds = 3, models = 1)
#' 
#' @importFrom dplyr starts_with group_by mutate_if mutate ungroup arrange
#'   full_join
#' @importFrom magrittr %>%
#' 
#' @export
predict_interactions <- function(
  features, gold_standard, classifier = c("NB", "SVM", "RF", "LR", "ensemble"),
  verbose = FALSE, models = 10, cv_folds = 10, trees = 500) {
  classifier <- match.arg(classifier)
  
  ## define global variables to prevent check complaining
  protein_A <- NULL; protein_B <- NULL; score.x <- NULL; score.y <- NULL;
  score.x.x <- NULL; score.y.y <- NULL
  
  # make labels
  if (verbose) {
    message("making labels ...")
  }
  labels <- make_labels(gold_standard, features)
  
  if (classifier %in% c("NB", "SVM", "RF", "LR")) {
    if (verbose) {
      message("training classifiers ...")
    }
    # predict with a classifier ensemble
    interactions = predict_ensemble(features, 
                                    labels, 
                                    classifier = classifier, 
                                    models = models, 
                                    cv_folds = cv_folds, 
                                    trees = trees)
  } else if (classifier == "ensemble") {
    # predict all four separately 
    if (verbose) {
      message("training naive Bayes classifiers ...")
    }
    interactions_NB <- predict_ensemble(
      features, labels, classifier = "NB", models, cv_folds) 
    if (verbose) {
      message("training random forest classifiers ...")
    }
    interactions_RF <- predict_ensemble(
      features, labels, classifier = "RF", models, cv_folds, trees) 
    if (verbose) {
      message("training support vector machine classifiers ...")
    }
    interactions_SVM <- predict_ensemble(
      features, labels, classifier = "SVM", models, cv_folds)
    if (verbose) {
      message("training logistic regression classifiers ...")
    }
    interactions_LR <- predict_ensemble(
      features, labels, classifier = "LR", models, cv_folds)
    if (verbose) {
      message("ensembling predictions ...")
    }
    interactions_list <- list(interactions_NB, interactions_LR,
                             interactions_RF, interactions_SVM)
    interactions <- Reduce(function(x, y) full_join(
      x, y, by = c('protein_A', 'protein_B')), interactions_list) %>%
      dplyr::select(-starts_with("label")) %>%
      mutate_if(is.numeric, ~ rank(-.)) %>%
      group_by(protein_A, protein_B) %>%
      mutate(mean = mean(c(score.x.x, score.y.y, score.x, score.y), 
                         na.rm = TRUE)) %>%
      ungroup() %>%
      arrange(mean)
    interactions$label <- make_labels(gold_standard, interactions)
  }
  
  # calculate precision
  interactions$precision <- calculate_precision(interactions$label)
  
  return(interactions)
}
