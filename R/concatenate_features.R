#' Combine features across multiple replicates
#' 
#' Concatenate features extracted from multiple replicates to a single data 
#' frame that will be used as input to a classifier. Doing so allows the 
#' classifier to naturally weight evidence for an interaction between each
#' protein pair from each feature in each replicate in proportion to its
#' discriminatory power on known examples.
#' 
#' @param feature_list a list of feature data frames, as produced by 
#'   \code{\link[PrInCE]{calculate_features}}, with proteins in the first two 
#'   columns
#' 
#' @return a data frame containing features for all protein pairs across all 
#'   replicates 
#' 
#' @importFrom dplyr full_join
#' 
#' @export
concatenate_features = function(feature_list) {
  features = Reduce(function(x, y)
    full_join(x, y, by = colnames(x)[1:2]), 
    feature_list)
  return(features)
}