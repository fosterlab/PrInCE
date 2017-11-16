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
#' @export
calculate_precision <- function(labels) {
  if (is.factor(labels))
    labels <- as.integer(labels) - 1
  TPs <- cumsum(replace(labels, is.na(labels), 0))
  FPs <- cumsum(replace(labels == 0, is.na(labels), 0))
  precision <- TPs / (TPs + FPs)
}
