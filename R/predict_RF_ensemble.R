#' Predict interactions using an ensemble of random forests
#' 
#' Use an ensemble of random forest classifiers to predict interactions from
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
#' @param models the number of random forest classifiers to train
#' @param cv_folds the number of folds to split the reference dataset into 
#' when training each random forests classifier. By default, each 
#' classifier uses ten-fold cross-validation, i.e., the classifier is trained
#' on 90\% of the dataset and used to classify the remaining 10\%
#' @param trees number of trees to grow for each fold
#' 
#' @return the input data frame of pairwise interactions, ranked by the 
#' median of classifier scores across all ensembled models
#' 
#' @export
predict_RF_ensemble <- function(input, labels, models = 1, cv_folds = 10,
                                seed = 0, trees = 500) {
  # set seed
  set.seed(seed)
  
  # replace missing data
  input <- replace_missing_data(input)
  
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
    rf_scores <- matrix(NA, ncol = cv_folds, nrow = n_interactions,
                        dimnames = list(interaction_names))
    for (fold in seq_len(cv_folds)) {
      # print message
      counter <- counter + 1
      pb$tick(tokens = list(what = sprintf(
        paste0("%-", nchar(total_models), "s"), counter)))

      # train model
      rf_data <- training[which(folds != fold),]
      rf_labels <- as.factor(training_labels[which(folds != fold)])
      rf_data_labeled <- cbind(rf_data, label = rf_labels)
      rf <- ranger::ranger(data = rf_data_labeled, 
                           dependent.variable.name = "label",
                           probability = T,
                           seed = seed,
                           num.trees = trees)
      
      # classify
      withheld_idxs = as.integer(rownames(training))[folds == fold]
      rf_test_data <- input[-withheld_idxs, -c(1:2)]
      predictions <- ranger:::predict.ranger(rf, rf_test_data, seed = seed)
      predictions <- predictions[[1]][, "1"]
      rf_scores[-withheld_idxs, fold] <- predictions
      
      # call GC
      gc()
    }
    medians <- setNames(robustbase::rowMedians(rf_scores, na.rm = T),
                        rownames(rf_scores))
    ensembled[, i] <- medians
  }
  
  # calculate median of medians across ensembled models
  ensembled_medians <- setNames(robustbase::rowMedians(ensembled, na.rm = T),
                                rownames(ensembled))
  
  # create ranked data frame
  interactions <- cbind(input[, 1:2], score = ensembled_medians, label = labels)
  interactions <- dplyr::arrange(interactions, -score)
  
  return(interactions)
}