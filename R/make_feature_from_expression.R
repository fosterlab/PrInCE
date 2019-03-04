#' Create a feature vector from expression data
#' 
#' Convert a gene or protein expression matrix into a feature vector
#' that matches the dimensions of a data frame used as input to a classifier, 
#' such as a naive Bayes, random forests, or support vector machine classifier,
#' by calculating the correlation between each pair of genes or proteins.  
#' 
#' @param expr a matrix containing gene or protein expression data, with 
#'   genes/proteins in columns and samples in rows
#' @param dat the data frame of features to be used by the classifier, 
#'    with protein pairs in the columns specified by the \code{node_columns}
#'    argument
#' @param node_columns a vector of length two, denoting either the indices 
#'   (integer vector) or column names (character vector) of the columns within 
#'   the data frame containing the nodes participating in pairwise interactions;
#'   defaults to the first two columns of the data frame (\code{c(1, 2)})
#' @param ... arguments passed to \code{cor}
#' 
#' @return a vector matching the dimensions and order of the feature data frame,
#' to use as input for a classifier in interaction prediction 
#' 
#' @export
make_feature_from_expression <- function(expr, dat, 
                                         node_columns = c(1, 2),
                                         ...) {
  # length of node columns must be exactly two (pairwise interactions)
  if (length(node_columns) != 2) {
    stop("length of `node_columns` must be exactly 2")
  }
  
  # get all proteins in feature data frame
  col1 <- node_columns[1]
  col2 <- node_columns[2]
  proteins_1 <- dat[[col1]]
  proteins_2 <- dat[[col2]]
  proteins <- unique(c(proteins_1, proteins_2))
  
  # subset expression to these proteins 
  filtered <- 1.0 * expr[, colnames(expr) %in% proteins]
  if (ncol(filtered) == 0)
    stop("no proteins overlap between expression and feature data")
  
  # create coexpression matrix
  coexpr <- cor(filtered, ...)
  
  # index matrix to get feature vector
  feat_idxs <- proteins_1 %in% rownames(coexpr) & 
    proteins_2 %in% rownames(coexpr)
  idxing_mat <- cbind(proteins_1[feat_idxs], proteins_2[feat_idxs])
  feature <- rep(NA, nrow(dat))
  feature[feat_idxs] <- coexpr[idxing_mat]
  return(feature)
}
