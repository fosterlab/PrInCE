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
#' @examples
#' ppi <- data.frame(protein_A = paste0("protein", seq_len(10)),
#'                   protein_B = paste0("protein", c(rep(3, 2), rep(5, 5), 
#'                                      rep(7, 3))))
#' adj <- adjacency_matrix_from_data_frame(ppi)
#' 
#' @export
adjacency_matrix_from_data_frame <- function(dat, symmetric = TRUE) {
  # remove self-interactions
  dat <- dat[as.character(dat[[1]]) != as.character(dat[[2]]), ]
  # create empty matrix
  proteins <- unique(c(as.character(dat[[1]]), as.character(dat[[2]])))
  n_proteins <- length(proteins)
  adjacency <- matrix(0, nrow = n_proteins, ncol = n_proteins,
                      dimnames = list(proteins, proteins))
  
  # identify pairwise interactions
  adjacency[as.matrix(dat[, c(1, 2)])] <- 1
  # add symmetric interactions
  if (symmetric) {
    adjacency[as.matrix(dat[, c(2, 1)])] <- 1
  }
  return(adjacency)
}
