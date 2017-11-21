#' Create a gold-standard label matrix
#' 
#' Create a gold-standard label matrix with the same dimensions as the 
#' feature matrices being classified. 
#' 
#' @param gold_standard an adjacency matrix of gold-standard interactions
#' @param profile_matrix the profile matrix for which interactions are being
#' predicted 
#' 
#' @return a matrix with the same dimensions as the input feature matrix
#' being used as input to the classifier 
#' 
#' @export
make_label_matrix <- function(gold_standard, profile_matrix) {
  proteins <- rownames(profile_matrix)
  label_mat <- matrix(NA, nrow = length(proteins), ncol = length(proteins),
                      dimnames = list(proteins, proteins))
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