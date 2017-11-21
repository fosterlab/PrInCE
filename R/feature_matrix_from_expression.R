#' Create a feature matrix for a naive Bayes classifier from expression data
#' 
#' Convert a gene or protein expression matrix into a feature matrix
#' with the same dimensions and row/column names as those used to train 
#' naive Bayes classifiers for default PrInCE features by calculating
#' the correlation between each pair of genes or proteins.  
#' 
#' @param expr a matrix containing gene or protein expression data, with 
#' genes/proteins in columns and samples in rows
#' @param profile_matrix the profile matrix for which interactions are being 
#' predicted 
#' @param ... arguments passed to \code{cor}
#' 
#' @return a square matrix with the same row and column names as the input 
#' profile matrix, for use in interaction prediction 
#' 
#' @export
feature_matrix_from_expression <- function(expr, profile_matrix, ...) {
  # get all proteins in profile matrix
  proteins <- rownames(profile_matrix)
  # subset data to these proteins 
  filtered <- expr[, colnames(expr) %in% proteins]
  if (ncol(filtered) == 0)
    stop("no proteins overlap between expression and profile matrices")
  # create empty matrix
  feature_matrix = matrix(NA, nrow = length(proteins), ncol = length(proteins),
                          dimnames = list(proteins, proteins))
  # fill coexpression
  coexpr <- cor(filtered, ...)
  # add empty rows and columns
  diff <- setdiff(proteins, colnames(filtered))
  new_rows <- matrix(NA, nrow = length(diff), ncol = ncol(coexpr),
                     dimnames = list(diff, colnames(coexpr)))
  coexpr <- rbind(coexpr, new_rows)
  new_cols <- matrix(NA, nrow = nrow(coexpr), ncol = length(diff),
                     dimnames = list(rownames(coexpr), diff))
  coexpr <- cbind(coexpr, new_cols)
  # same order as profile matrix
  coexpr <- coexpr[rownames(profile_matrix), rownames(profile_matrix)]
  return(coexpr)
}