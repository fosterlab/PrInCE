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
#' @examples
#' data(scott)
#' data(scott_gaussians)
#' data(gold_standard)
#' # analyze only the first 100 profiles
#' subset <- scott[seq_len(500), ]
#' gauss <- scott_gaussians[names(scott_gaussians) %in% rownames(subset)]
#' ppi <- PrInCE(subset, gold_standard, gaussians = gauss, models = 1, 
#'              cv_folds = 3)
#' network <- threshold_precision(ppi, threshold = 0.5)
#' nrow(network)
#' 
#' @export
threshold_precision <- function(interactions, threshold) {
  n_keep <- max(which(interactions$precision >= threshold), na.rm = TRUE)
  interactions[seq_len(n_keep), ]
}