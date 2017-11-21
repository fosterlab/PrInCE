#' Create an adjacency matrix from a data frame
#' 
#' Convert a data frame containing pairwise interactions into an
#' adjacency matrix. The resulting square adjacency matrix contains ones for 
#' proteins that are found in interactions and zeroes otherwise.
#' 
#' @param dat a data frame with pairwise interactions in the first two columns
#' @param symmetric if true, interactions in both directions will be 
#' added to the adjacency matrix 
#' 
#' @return an adjacency matrix between all interacting proteins
#' 
#' @export
adjacency_matrix_from_data_frame <- function(dat, symmetric = T) {
  proteins <- unique(c(dat[[1]], dat[[2]]))
  n_proteins <- length(proteins)
  adjacency <- matrix(0, nrow = n_proteins, ncol = n_proteins,
                      dimnames = list(proteins, proteins))
  
  # identify pairwise interactions
  adjacency[as.matrix(dat[, 1:2])] <- 1
  # add symmetric interactions
  if (symmetric) {
    adjacency[as.matrix(dat[, 2:1])] <- 1
  }
  return(adjacency)
}