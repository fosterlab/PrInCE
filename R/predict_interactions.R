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
#' @param seed the seed for random number generation, to ensure reproducible
#'   results
#' 
#' @return a ranked data frame of pairwise interactions, with the 
#' classifier score, label, and cumulative precision for each interaction 
#' 
#' @importFrom dplyr starts_with group_by mutate_if mutate ungroup arrange
#'   full_join
#' @importFrom magrittr %>%
#' 
#' @export
predict_interactions = function(
  features, gold_standard, classifier = c("NB", "SVM", "RF", "LR", "ensemble"),
  verbose = F, models = 10, cv_folds = 10, trees = 500, seed = 0) {
  classifier = match.arg(classifier)
  
  ## define global variables to prevent check complaining
  protein_A = NULL; protein_B = NULL; score.x = NULL; score.y = NULL;
  score.x.x = NULL; score.y.y = NULL
  
  # make labels
  if (verbose) {
    message("making labels ...")
  }
  labels = make_labels(gold_standard, features)
  
  if (classifier %in% c("NB", "SVM", "RF", "LR")) {
    if (verbose) {
      message("training classifiers ...")
    }
    # predict with a classifier ensemble
    if (classifier == "NB") {
      interactions = predict_NB_ensemble(
        features, labels, models, cv_folds, seed)
    } else if (classifier == "RF") {
      interactions = predict_RF_ensemble(
        features, labels, models, cv_folds, seed, trees) 
    } else if (classifier == "SVM") {
      interactions = predict_SVM_ensemble(
        features, labels, models, cv_folds, seed)
    } else if (classifier == "LR") {
      interactions = predict_logistic_regression_ensemble(
        features, labels, models, cv_folds, seed)
    }
  } else if (classifier == "ensemble") {
    # predict all four separately 
    if (verbose) {
      message("training naive Bayes classifiers ...")
    }
    interactions_NB = predict_NB_ensemble(
      features, labels, models, cv_folds, seed) 
    if (verbose) {
      message("training random forest classifiers ...")
    }
    interactions_RF = predict_RF_ensemble(
      features, labels, models, cv_folds, seed, trees) 
    if (verbose) {
      message("training support vector machine classifiers ...")
    }
    interactions_SVM = predict_SVM_ensemble(
      features, labels, models, cv_folds, seed)
    if (verbose) {
      message("training logistic regression classifiers ...")
    }
    interactions_LR = predict_logistic_regression_ensemble(
      features, labels, models, cv_folds, seed)
    if (verbose) {
      message("ensembling predictions ...")
    }
    interactions_list = list(interactions_NB, interactions_LR,
                             interactions_RF, interactions_SVM)
    interactions = Reduce(function(x, y) full_join(
      x, y, by = c('protein_A', 'protein_B')), interactions_list) %>%
      dplyr::select(-starts_with("label")) %>%
      mutate_if(is.numeric, ~ rank(-.)) %>%
      group_by(protein_A, protein_B) %>%
      mutate(mean = mean(c(score.x.x, score.y.y, score.x, score.y), 
                         na.rm = T)) %>%
      ungroup() %>%
      arrange(mean)
    interactions$label = make_labels(gold_standard, interactions)
  }
  
  # calculate precision
  interactions$precision = calculate_precision(interactions$label)
  
  return(interactions)
}
