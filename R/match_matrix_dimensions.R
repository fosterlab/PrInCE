#' Match the dimensions of a query matrix to a profile matrix
#' 
#' Match the row and column names of a square feature matrix to the row names
#' of a profile matrix, adding rows/columns containing \code{NA}s when 
#' proteins in the profile matrix are missing from the feature matrix. 
#' 
#' @param query a square matrix containing features for pairs of proteins 
#' @param profile_matrix the profile matrix for which interactions are being 
#' predicted 
#' 
#' @return a square matrix with the same row and column names as the input 
#' profile matrix, for use in interaction prediction 
#' 
#' @examples 
#' data(gold_standard)
#' subset = adjacency_matrix_from_list(gold_standard[seq(1, 200)])
#' target = adjacency_matrix_from_list(gold_standard)
#' matched = match_matrix_dimensions(query, target)
#' dim(subset)
#' dim(target)
#' dim(matched)
#' 
#' @export
match_matrix_dimensions <- function(query, profile_matrix) {
  proteins <- rownames(profile_matrix)
  # add empty rows and columns
  diff <- setdiff(proteins, colnames(query))
  new_rows <- matrix(NA, nrow = length(diff), ncol = ncol(query),
                     dimnames = list(diff, colnames(query)))
  query <- rbind(query, new_rows)
  new_cols <- matrix(NA, nrow = nrow(query), ncol = length(diff),
                     dimnames = list(rownames(query), diff))
  query <- cbind(query, new_cols)
  # same order as profile matrix
  query <- query[rownames(profile_matrix), rownames(profile_matrix)]
  return(query)
}
