#' Predict interactions using an ensemble of naive Bayes classifiers 
#' 
#' Use an ensemble of naive Bayes classifiers to predict interactions from
#' co-elution dataset features. The ensemble approach ensures that 
#' results are robust to the partitioning of the dataset into folds. For each
#' model, the median of classifier scores across all folds is calculated.
#' Then, the median of all such medians across all models is calculated. 
#' 
#' @param input a data frame containing interacting gene/protein pairs in the
#' first two columns, and the features to use for classification in the 
#' remaining columns
#' @param labels labels for each interaction in \code{input}: 0 for negatives,
#' 1 for positives, and NA for interactions outside the reference set 
#' @param models the number of naive Bayes classifiers to train
#' @param cv_folds the number of folds to split the reference dataset into 
#' when training each naive Bayes classifier. By default, each naive Bayes
#' classifier uses ten-fold cross-validation, i.e., the classifier is trained
#' on 90\% of the dataset and used to classify the remaining 10\%
#' 
#' @return the input data frame of pairwise interactions, ranked by the 
#' median of classifier scores across all ensembled models
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
#' 
#' @importFrom stats predict
#' 
#' @export
predict_NB_ensemble <- function(input, labels, models = 10, cv_folds = 10) {
  ## define global variables to prevent check complaining
  score = NULL
  
  # extract training data
  training_idxs <- which(!is.na(labels))
  training_labels <- as.factor(labels[training_idxs])
  training <- input[training_idxs, -c(1, 2)]
  
  # create matrix to hold medians from each model
  n_interactions = nrow(input)
  interaction_names <- paste0(input[[1]], "_", input[[2]])
  ensembled <- matrix(NA, ncol = models, nrow = n_interactions,
                      dimnames = list(interaction_names))
  
  # create progress bar
  total_models <- models * cv_folds
  pb <- progress::progress_bar$new(
    format = "running fold :what [:bar] :percent eta: :eta",
    clear = FALSE, total = total_models, width = 80)
  counter <- 0
  
  # create models
  for (i in seq_len(models)) {
    folds <- cut(seq_len(nrow(training)), breaks = cv_folds, labels = FALSE)
    folds <- sample(folds) ## randomize
    nb_scores <- matrix(NA, ncol = cv_folds, nrow = n_interactions,
                        dimnames = list(interaction_names))
    for (fold in seq_len(cv_folds)) {
      # print message
      counter <- counter + 1
      pb$tick(tokens = list(what = sprintf(
        paste0("%-", nchar(total_models), "s"), counter)))

      # train model
      nb_data <- training[which(folds != fold),]
      nb_labels <- as.factor(training_labels[which(folds != fold)])
      nb = naivebayes::naive_bayes(nb_data, nb_labels)
      
      # classify
      withheld_idxs = as.integer(rownames(training))[folds == fold]
      predictions = predict(
        nb, input[-withheld_idxs, -c(1, 2)], type = 'prob', threshold = 1e-10)
      predictions = predictions[, "1"]
      nb_scores[-withheld_idxs, fold] <- predictions
      
      # call GC
      gc()
    }
    medians <- setNames(robustbase::rowMedians(nb_scores, na.rm = TRUE),
                        rownames(nb_scores))
    ensembled[, i] <- medians
  }
  
  # calculate median of medians across ensembled models
  ensembled_medians <- setNames(robustbase::rowMedians(ensembled, na.rm = TRUE),
                                rownames(ensembled))
  
  # create ranked data frame
  interactions <- cbind(input[, c(1, 2)], score = ensembled_medians, 
                        label = labels)
  interactions <- dplyr::arrange(interactions, -score)
  
  return(interactions)
}