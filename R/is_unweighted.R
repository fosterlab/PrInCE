#' Test whether a network is unweighted
#'  
#' @param network the network to analyze 
#' 
#' @return true if the input network is a square logical or numeric matrix
#' 
#' @examples 
#' data(gold_standard)
#' adj <- adjacency_matrix_from_list(gold_standard)
#' is_unweighted(adj) ## returns TRUE
#' 
#' @importFrom tester is_square_matrix is_logical_matrix is_numeric_matrix
#' 
#' @export
is_unweighted <- function(network) {
  # check input is square logical, integer, or numeric matrix
  square <- is_square_matrix(network)
  if (!square)
    stop("input is not a square matrix")
  lgl <- is_logical_matrix(network)
  dbl <- is_numeric_matrix(network)
  if (!lgl & !dbl)
    stop("input could not be converted to adjacency matrix")
  # convert to numeric 
  network <- network * 1
  # all values should be missing, zero, or one
  all(is.na(network) | network == 1 | network == 0)
}