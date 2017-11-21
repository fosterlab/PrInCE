#' Create a feature matrix for a naive Bayes classifier from a data frame
#' 
#' Convert a data frame containing pairwise interactions, and  
#' a score or other data associated with each interaction, into a feature matrix
#' with the same dimensions and row/column names as those used to train 
#' naive Bayes classifiers for default PrInCE features.
#' 
#' @param dat a data frame with pairwise interactions in the first two columns,
#' and optionally some information about each interaction (such as a score or
#' category) in the third column
#' @param profile_matrix the profile matrix for which interactions are being 
#' predicted 
#' @param col_name the name of the column in the input data frame that contains
#' interaction features 
#' 
#' @return a square matrix with the same row and column names as the input 
#' profile matrix, for use in interaction prediction 
#' 
#' @export
feature_matrix_from_data_frame <- function(dat, profile_matrix, col_name) {
  # get all proteins in profile matrix
  proteins <- rownames(profile_matrix)
  # subset data frame to these proteins 
  overlap_idxs <- dat[[1]] %in% proteins & dat[[2]] %in% proteins
  filtered <- dat[overlap_idxs,]
  if (nrow(filtered) == 0)
    return(NULL)
  # create empty matrix
  feature_matrix = matrix(0, nrow = length(proteins), ncol = length(proteins),
                          dimnames = list(proteins, proteins))
  # fill interaction scores
  if (!is.null(col_name)) {
    feature_column <- filtered[[col_name]]
  } else {
    feature_column <- filtered[[3]]
  }
  feature_matrix[as.matrix(filtered[, 1:2])] <- feature_column
  return(feature_matrix)
}