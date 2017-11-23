#' Create a feature vector for a classifier from a data frame
#' 
#' Convert a data frame containing pairwise interactions, and  
#' a score or other data associated with each interaction, into a feature vector
#' that matches the dimensions of a data frame used as input to a classifier, 
#' such as a naive Bayes, random forests, or support vector machine classifier.
#' 
#' @param dat a data frame with pairwise interactions in the first two columns,
#' and the features to convert to a vector in some other column
#' @param input the data frame of features to be used by the classifier, 
#' with protein pairs in the first two columns
#' @param col_name the name of the column in the first data frame that 
#' contains interaction features 
#' 
#' @return a vector matching the dimensions and order of the feature data frame,
#' to use as input for a classifier in interaction prediction 
#' 
#' @export
make_feature_from_data_frame <- function(dat, input, col_name) {
  # get all proteins in feature data frame
  ref_proteins1 <- input[[1]]
  ref_proteins2 <- input[[2]]
  ref_proteins <- unique(c(ref_proteins1, ref_proteins2))
  # subset data frame to these proteins 
  overlap_idxs <- dat[[1]] %in% ref_proteins & dat[[2]] %in% ref_proteins
  if (sum(overlap_idxs) == 0)
    stop("no proteins overlap between original and new features")
  filtered <- dat[overlap_idxs,]
  # convert the first data frame to a matrix
  new_proteins1 <- filtered[[1]]
  new_proteins2 <- filtered[[2]] 
  new_proteins <- unique(c(new_proteins1, new_proteins2))
  feature_matrix <- matrix(NA, nrow = length(new_proteins), 
                           ncol = length(new_proteins),
                           dimnames = list(new_proteins, new_proteins))
  # fill interaction scores
  if (!is.null(col_name)) {
    feature_column <- filtered[[col_name]]
  } else {
    feature_column <- filtered[[3]]
  }
  feature_matrix[as.matrix(filtered[, 1:2])] <- feature_column
  # index matrix 
  feat_idxs <- ref_proteins1 %in% rownames(feature_matrix) & 
    ref_proteins2 %in% rownames(feature_matrix)
  idxing_mat <- cbind(ref_proteins1[feat_idxs], ref_proteins2[feat_idxs])
  feature <- rep(NA, nrow(input))
  feature[feat_idxs] <- feature_matrix[idxing_mat]
  return(feature)
}