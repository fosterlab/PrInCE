#' Create a gold-standard label matrix
#' 
#' Create a gold-standard label matrix with the same dimensions as the 
#' feature matrices being classified. 
#' 
#' @param gold_standard an adjacency matrix of gold-standard interactions
#' @param feature_matrix any of the feature matrices being used as input to a 
#' classifier 
#' 
#' @return a matrix with the same dimensions as the input feature matrix
#' being used as input to the classifier 
make_label_matrix <- function(gold_standard, feature_matrix) {
  proteins <- rownames(feature_matrix)
  label_mat <- feature_matrix
  label_mat[] <- NA
  overlap_idxs <- rownames(gold_standard) %in% proteins
  gold_standard <- gold_standard[overlap_idxs, overlap_idxs]
  tri <- upper.tri(gold_standard)
  idxs <- which(tri, arr.ind = T)
  pairs <- data.frame(protein_A = rownames(gold_standard)[idxs[, 1]], 
                      protein_B = rownames(gold_standard)[idxs[, 2]],
                      values = gold_standard[tri])
  label_mat[as.matrix(pairs[, 1:2])] <- pairs$values
  return(label_mat)
}