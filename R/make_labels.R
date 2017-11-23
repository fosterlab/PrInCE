#' Make labels for a classifier based on a gold standard
#' 
#' Create labels for a classifier for protein pairs in the same order as 
#' in a dataset that will be used as input to a classifier, in a 
#' memory-friendly way. 
#' 
#' @param gold_standard an adjacency matrix of gold-standard interactions
#' @param input a data frame with interacting proteins in the first two
#' columns 
#' 
#' @return a vector of the same length as the input dataset, containing 
#' \code{NA}s for protein pairs not in the gold standard and ones or zeroes
#' based on the content of the adjacency matrix
#' 
#' @export
make_labels <- function(gold_standard, input) {
  proteins_1 <- input[[1]]
  proteins_2 <- input[[2]]
  lab_idxs <- proteins_1 %in% rownames(gold_standard) &
    proteins_2 %in% rownames(gold_standard)
  if (sum(lab_idxs) == 0)
    stop("no proteins overlap between input and gold standard")
  idxing_mat <- cbind(proteins_1[lab_idxs], proteins_2[lab_idxs])
  labels <- rep(NA, nrow(input))
  labels[lab_idxs] <- gold_standard[idxing_mat]
  return(labels)
}