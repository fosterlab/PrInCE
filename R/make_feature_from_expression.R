#' Create a feature vector from expression data
#' 
#' Convert a gene or protein expression matrix into a feature vector
#' that matches the dimensions of a data frame used as input to a classifier, 
#' such as a naive Bayes, random forests, or support vector machine classifier,
#' by calculating the correlation between each pair of genes or proteins.  
#' 
#' @param expr a matrix containing gene or protein expression data, with 
#' genes/proteins in columns and samples in rows
#' @param dat the data frame of features to be used by the classifier, 
#' with protein pairs in the first two columns
#' @param ... arguments passed to \code{cor}
#' 
#' @return a vector matching the dimensions and order of the feature data frame,
#' to use as input for a classifier in interaction prediction 
#' 
#' @export
make_feature_from_expression <- function(expr, dat, ...) {
  # get all proteins in feature data frame
  proteins_1 <- dat[[1]]
  proteins_2 <- dat[[2]]
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
