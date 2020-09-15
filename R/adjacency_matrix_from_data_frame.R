#' Create an adjacency matrix from a data frame
#' 
#' Convert a data frame containing pairwise interactions into an
#' adjacency matrix. The resulting square adjacency matrix contains ones for 
#' proteins that are found in interactions and zeroes otherwise.
#' 
#' @param dat a data frame containing pairwise interactions
#' @param symmetric if true, interactions in both directions will be 
#'   added to the adjacency matrix
#' @param node_columns a vector of length two, denoting either the indices 
#'   (integer vector) or column names (character vector) of the columns within 
#'   the data frame containing the nodes participating in pairwise interactions;
#'   defaults to the first two columns of the data frame (\code{c(1, 2)})
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
adjacency_matrix_from_data_frame <- function(dat, 
                                             symmetric = TRUE,
                                             node_columns = c(1, 2)
                                             ) {
  # length of node columns must be exactly two (pairwise interactions)
  if (length(node_columns) != 2) {
    stop("length of `node_columns` must be exactly 2")
  }
    
  # remove self-interactions
  col1 <- node_columns[1]
  col2 <- node_columns[2]
  dat <- dat[dat[[col1]] != dat[[col2]], ]
  nodes1 <- as.character(dat[[col1]])
  nodes2 <- as.character(dat[[col2]])
  
  # create empty matrix
  proteins <- unique(c(nodes1, nodes2))
  n_proteins <- length(proteins)
  adjacency <- matrix(0, nrow = n_proteins, ncol = n_proteins,
                      dimnames = list(proteins, proteins))
  
  # identify pairwise interactions
  adjacency[as.matrix(dat[, node_columns])] <- 1
  # add symmetric interactions
  if (symmetric) {
    adjacency[as.matrix(dat[, rev(node_columns)])] <- 1
  }
  return(adjacency)
}
