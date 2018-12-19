#' Threshold interactions at a given precision cutoff
#'  
#' @param interactions the ranked list of interactions output by 
#'   \code{\link[PrInCE]{predict_interactions}}, including a \code{precision} 
#'   column
#' @param threshold the minimum precision of the unweighted interaction 
#'   network to return
#' 
#' @return the subset of the original ranked list at the given precision
#' 
#' @export
threshold_precision = function(interactions, threshold) {
  n_keep = max(which(interactions$precision >= threshold), na.rm = T)
  interactions[seq_len(n_keep), ]
}