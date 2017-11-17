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
#' @return the median of classifier scores across all models for each row 
#' in the input data fame 
#' 
#' @export
predict_NB_ensemble <- function(input, labels, models = 10, cv_folds = 10,
                                seed = 0) {
  # set seed
  set.seed(seed)
  
  # extract training data
  training_idxs <- which(!is.na(labels))
  training_labels <- as.factor(labels[training_idxs])
  training <- input[training_idxs, -c(1:2)]
  
  # create matrix to hold medians from each model
  n_interactions = nrow(input)
  interaction_names <- paste0(input[[1]], "_", input[[2]])
  ensembled <- matrix(NA, ncol = models, nrow = n_interactions,
                      dimnames = list(interaction_names))
  
  # create progress bar
  total_models <- models * cv_folds
  pb <- progress::progress_bar$new(
    format = "running fold :what [:bar] :percent eta: :eta",
    clear = F, total = total_models, width = 100)
  counter <- 0
  
  # create models
  for (i in seq_len(models)) {
    folds <- cut(seq_len(nrow(training)), breaks = cv_folds, labels = F)
    folds <- sample(folds) ## randomize
    nb_scores <- matrix(NA, ncol = cv_folds, nrow = n_interactions,
                        dimnames = list(interaction_names))
    for (fold in seq_len(cv_folds)) {
      # print message
      counter <- counter + 1
      pb$tick(tokens = list(what = sprintf(
        paste0("%-", nchar(total_models), "s"), counter)))

      # train model
      nb = naivebayes::naive_bayes(training[which(folds == fold),],
                                   as.factor(labels[which(folds == fold)]))
      
      # classify
      withheld_idxs = as.integer(rownames(training))[folds != fold]
      predictions = naivebayes:::predict.naive_bayes(
        nb, input[-withheld_idxs, -c(1:2)], type = 'prob')
      predictions = predictions[, "1"]
      nb_scores[-withheld_idxs, fold] <- predictions
    }
    medians <- setNames(robustbase::rowMedians(nb_scores, na.rm = T),
                        rownames(nb_scores))
    ensembled[, i] <- medians
  }
  
  ensembled_medians <- setNames(robustbase::rowMedians(ensembled, na.rm = T),
                                rownames(ensembled))
  return(ensembled_medians)
}