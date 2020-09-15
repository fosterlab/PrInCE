#' Predict interactions using an ensemble of classifiers 
#' 
#' Use an ensemble of classifiers to predict interactions from
#' co-elution dataset features. The ensemble approach ensures that 
#' results are robust to the partitioning of the dataset into folds. For each
#' model, the median of classifier scores across all folds is calculated.
#' Then, the median of all such medians across all models is calculated. 
#' 
#' @param dat a data frame containing interacting gene/protein pairs in the
#'   first two columns, and the features to use for classification in the 
#'   remaining columns
#' @param labels labels for each interaction in \code{dat}: 0 for negatives,
#'   1 for positives, and NA for interactions outside the reference set 
#' @param classifier the type of classifier to use; one of \code{"NB"} 
#'   (naive Bayes), \code{"SVM"} (support vector machine), \code{"RF"}
#'   (random forest), or \code{"LR"} (logistic regression)
#' @param models the number of classifiers to train
#' @param cv_folds the number of folds to split the reference dataset into 
#'   when training each classifier. By default, each 
#'   classifier uses ten-fold cross-validation, i.e., the classifier is trained
#'   on 90\% of the dataset and used to classify the remaining 10\%
#' @param trees for random forest classifiers only, the number of trees to 
#'   grow for each fold
#' @param node_columns a vector of length two, denoting either the indices 
#'   (integer vector) or column names (character vector) of the columns within 
#'   the input data frame containing the nodes participating in pairwise 
#'   interactions; defaults to the first two columns of the data frame 
#'   (\code{c(1, 2)})
#'  
#' @return the input data frame of pairwise interactions, ranked by the 
#' median of classifier scores across all ensembled models
#' 
#' @examples
#' ## calculate features
#' data(scott)
#' data(scott_gaussians)
#' subset <- scott[seq_len(500), ] ## limit to first 500 proteins
#' gauss <- scott_gaussians[names(scott_gaussians) %in% rownames(subset)]
#' features <- calculate_features(subset, gauss)
#' ## make training labels
#' data(gold_standard)
#' ref <- adjacency_matrix_from_list(gold_standard)
#' labels <- make_labels(ref, features)
#' ## predict interactions with naive Bayes classifier
#' ppi <- predict_ensemble(features, labels, classifier = "NB", 
#'                         cv_folds = 3, models = 1)
#' 
#' @importFrom stats predict binomial
#' @importFrom robustbase rowMedians
#' @importFrom progress progress_bar
#' @importFrom naivebayes naive_bayes
#' @importFrom LiblineaR LiblineaR
#' @importFrom ranger ranger
#' @importFrom speedglm speedglm
#' @importFrom dplyr arrange
#' 
#' @export
predict_ensemble <- function(dat, 
                             labels, 
                             classifier = c("NB", "SVM", "RF", "LR"), 
                             models = 1, 
                             cv_folds = 10,
                             trees = 500,
                             node_columns = c(1, 2)) {
  classifier <- match.arg(classifier)
  # length of node columns must be exactly two (pairwise interactions)
  if (length(node_columns) != 2) {
    stop("length of `node_columns` must be exactly 2")
  }
  
  ## define global variables to prevent check complaining
  score <- NULL
  
  # scale all features
  if (is.numeric(node_columns)) {
    node_colnames <- colnames(dat)[node_columns]
  } else if (is.character(node_columns)) {
    node_colnames <- node_columns
  } else {
    stop("`node_columns` must be an integer or character vector")
  }
  if (classifier == "SVM") {
    dat[, !colnames(dat) %in% node_colnames] <- sapply(
      dat[, !colnames(dat) %in% node_colnames], scale)
  }
  
  # extract training data
  training_idxs <- which(!is.na(labels))
  training_labels <- as.factor(labels[training_idxs])
  training <- dat[training_idxs, !colnames(dat) %in% node_colnames]
  
  # create matrix to hold medians from each model
  n_interactions <- nrow(dat)
  col1 <- node_columns[1]
  col2 <- node_columns[2]
  interaction_names <- paste0(dat[[col1]], "_", dat[[col2]])
  ensembled <- matrix(NA, ncol = models, nrow = n_interactions,
                      dimnames = list(interaction_names))
  
  # create progress bar
  total_models <- models * cv_folds
  pb <- progress_bar$new(
    format = "running fold :what [:bar] :percent eta: :eta",
    clear = FALSE, total = total_models, width = 80)
  counter <- 0
  
  # create models
  for (i in seq_len(models)) {
    folds <- cut(seq_len(nrow(training)), breaks = cv_folds, labels = FALSE)
    folds <- sample(folds) ## randomize
    clf_scores <- matrix(NA, ncol = cv_folds, nrow = n_interactions,
                         dimnames = list(interaction_names))
    for (fold in seq_len(cv_folds)) {
      # print message
      counter <- counter + 1
      pb$tick(tokens = list(what = sprintf(
        paste0("%-", nchar(total_models), "s"), counter)))
      
      # train model
      clf_data <- training[which(folds != fold),]
      clf_labels <- as.factor(training_labels[which(folds != fold)])
      clf_data_labeled <- cbind(clf_data, label = clf_labels)
      clf <- switch(classifier,
                    NB = naive_bayes(clf_data, clf_labels),
                    SVM = LiblineaR(clf_data, clf_labels, type = 2),
                    RF = ranger(data = clf_data_labeled, 
                                dependent.variable.name = "label",
                                probability = TRUE,
                                num.trees = trees,
                                num.threads = 1),
                    LR = speedglm(label ~ ., clf_data_labeled, 
                                  family = binomial()))
      
      # classify
      withheld_idxs <- as.integer(rownames(training))[folds == fold]
      test_data <- dat[-withheld_idxs, !colnames(dat) %in% node_colnames]
      predictions <- switch(
        classifier, 
        NB = predict(clf, test_data, type = 'prob', threshold = 5e-324),
        SVM = predict(clf, test_data, decisionValues = TRUE),
        RF = predict(clf, test_data, num.threads = 1),
        LR = predict(clf, test_data, type = 'response'))
      ## extract predictions as numeric vector
      predictions <- switch(
        classifier,
        NB = predictions[, "1"],
        SVM = -1.0 * predictions$decisionValues[, "0"],
        RF = predictions[[1]][, "1"],
        LR = predictions)
      # assign scores
      clf_scores[-withheld_idxs, fold] <- predictions
      
      # call GC
      gc()
    }
    medians <- setNames(rowMedians(clf_scores, na.rm = TRUE),
                        rownames(clf_scores))
    ensembled[, i] <- medians
  }
  
  # calculate median of medians across ensembled models
  ensembled_medians <- setNames(rowMedians(ensembled, na.rm = TRUE),
                                rownames(ensembled))
  
  # create ranked data frame
  interactions <- cbind(dat[, node_columns], score = ensembled_medians, 
                        label = labels)
  interactions <- arrange(interactions, -score)
  
  return(interactions)
}