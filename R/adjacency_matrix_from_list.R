#' Create an adjacency matrix from a list of complexes
#' 
#' Convert a list of complexes into a pairwise adjacency matrix. The resulting
#' square adjacency matrix contains ones for proteins that are found in the
#' same complex and zeroes otherwise.
#' 
#' @param complexes a list of complexes, with each entry containing complex
#' subunits as a character vector
#' 
#' @return an adjacency matrix between all complex subunits
#' 
#' @examples
#' data(gold_standard)
#' adj <- adjacency_matrix_from_list(gold_standard)
#' 
#' @importFrom purrr map_dfr
#' 
#' @export
adjacency_matrix_from_list <- function(complexes) {
  proteins <- unique(unlist(complexes))
  n_proteins <- length(proteins)
  adjacency <- matrix(0, nrow = n_proteins, ncol = n_proteins,
                      dimnames = list(proteins, proteins))
  
  # identify pairwise interactions
  pairs <- map_dfr(complexes, ~ tidyr::crossing(., .))
  adjacency[as.matrix(pairs)] <- 1
  
  return(adjacency)
}