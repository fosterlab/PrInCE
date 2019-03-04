#' Create a feature vector for a classifier from a data frame
#' 
#' Convert a data frame containing pairwise interactions, and  
#' a score or other data associated with each interaction, into a feature vector
#' that matches the dimensions of a data frame used as input to a classifier, 
#' such as a naive Bayes, random forests, or support vector machine classifier.
#' 
#' @param dat a data frame containing pairwise interactions and a feature to be
#'   converted to a vector in a third column
#' @param target the data frame of features that will be provided as input to a 
#'   classifier
#' @param dat_node_cols a vector of length two, denoting either the indices 
#'   (integer vector) or column names (character vector) of the columns within 
#'   the feature data frame; defaults to the first two columns of the data 
#'   frame (\code{c(1, 2)})
#' @param target_node_cols a vector of length two, denoting either the indices 
#'   (integer vector) or column names (character vector) of the columns within 
#'   the target data frame; defaults to the first two columns of the data 
#'   frame (\code{c(1, 2)})
#' @param feature_col the name or index of the column in the first data frame 
#'   that contains a feature for each interaction
#' @param default_value the default value for protein pairs that are not in the 
#'   first data frame (set, by default, to \code{NA})
#' 
#' @return a vector matching the dimensions and order of the feature data frame,
#' to use as input for a classifier in interaction prediction 
#' 
#' @export
make_feature_from_data_frame <- function(dat, 
                                         target, 
                                         dat_node_cols = c(1, 2),
                                         target_node_cols = c(1, 2),
                                         feature_col = 3, 
                                         default_value = NA) {
  # length of column names/indices must be exactly two (pairwise interactions)
  if (length(dat_node_cols) != 2) {
    stop("length of `dat_node_cols` must be exactly 2")
  }
  if (length(target_node_cols) != 2) {
    stop("length of `target_node_cols` must be exactly 2")
  }
  
  # get all proteins in feature data frame
  target_col1 <- target_node_cols[1]
  target_col2 <- target_node_cols[2]
  ref_proteins1 <- target[[target_col1]]
  ref_proteins2 <- target[[target_col2]]
  ref_proteins <- unique(c(ref_proteins1, ref_proteins2))
  
  # subset data frame to these proteins 
  query_col1 = dat_node_cols[1]
  query_col2 = dat_node_cols[2]
  overlap_idxs <- dat[[query_col1]] %in% ref_proteins & 
    dat[[query_col2]] %in% ref_proteins
  if (sum(overlap_idxs) == 0)
    stop("no proteins overlap between original and new features")
  filtered <- dat[overlap_idxs, ]
  
  # convert the first data frame to a matrix
  new_proteins1 <- filtered[[query_col1]]
  new_proteins2 <- filtered[[query_col2]] 
  new_proteins <- unique(c(new_proteins1, new_proteins2))
  feature_matrix <- matrix(default_value, nrow = length(new_proteins), 
                           ncol = length(new_proteins),
                           dimnames = list(new_proteins, new_proteins))
  
  # fill interaction scores
  feature_column <- filtered[[col_name]]
  feature_matrix[as.matrix(filtered[, dat_node_cols])] <- feature_column
  
  # index matrix to get feature vector
  feat_idxs <- ref_proteins1 %in% rownames(feature_matrix) & 
    ref_proteins2 %in% rownames(feature_matrix)
  idxing_mat <- cbind(ref_proteins1[feat_idxs], ref_proteins2[feat_idxs])
  feature <- rep(default_value, nrow(target))
  feature[feat_idxs] <- feature_matrix[idxing_mat]
  return(feature)
}
