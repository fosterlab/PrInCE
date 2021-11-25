#' Make labels for a classifier based on a gold standard
#' 
#' Create labels for a classifier for protein pairs in the same order as 
#' in a dataset that will be used as input to a classifier, in a 
#' memory-friendly way. 
#' 
#' @param gold_standard an adjacency matrix of gold-standard interactions
#' @param dat a data frame with interacting proteins in the first two
#'   columns 
#' @param node_columns a vector of length two, denoting either the indices 
#'   (integer vector) or column names (character vector) of the columns within 
#'   the data frame containing the nodes participating in pairwise interactions;
#'   defaults to the first two columns of the data frame (\code{c(1, 2)})
#' @param protein_groups optionally, specify a list linking each protein in the 
#'   first two columns of the input data frame to a protein group 
#' 
#' @return a vector of the same length as the input dataset, containing 
#' \code{NA}s for protein pairs not in the gold standard and ones or zeroes
#' based on the content of the adjacency matrix
#' 
#' @examples 
#' data(gold_standard)
#' adj <- adjacency_matrix_from_list(gold_standard)
#' proteins <- unique(unlist(gold_standard))
#' dat <- data.frame(protein_A = sample(proteins, 10), 
#'                   protein_B = sample(proteins, 10))
#' labels <- make_labels(adj, dat)
#' 
#' @importFrom purrr map map_lgl
#' @importFrom tidyr crossing
#' 
#' @export
make_labels <- function(gold_standard, dat, 
                        node_columns = c(1, 2),
                        protein_groups = NULL) {
  # length of node columns must be exactly two (pairwise interactions)
  if (length(node_columns) != 2) {
    stop("length of `node_columns` must be exactly 2")
  }
  col1 <- node_columns[1]
  col2 <- node_columns[2]
  proteins_1 <- as.character(dat[[col1]])
  proteins_2 <- as.character(dat[[col2]])
  if (is.null(protein_groups)) {
    lab_idxs <- proteins_1 %in% rownames(gold_standard) &
      proteins_2 %in% rownames(gold_standard)
    if (sum(lab_idxs) == 0)
      stop("no proteins overlap between input and gold standard")
    idxing_mat <- cbind(proteins_1[lab_idxs], proteins_2[lab_idxs])
    labels <- rep(NA, nrow(dat))
    labels[lab_idxs] <- gold_standard[idxing_mat]
  } else {
    names(protein_groups) = unlist(sapply(protein_groups, function(x) x[1]))
    # filter groups based on presence/absence in matrix
    protein_groups <- map(protein_groups, ~ .[. %in% rownames(gold_standard)])
    # identify rows that overlap with the gold standard matrix
    groups_1 <- protein_groups[proteins_1]
    groups_2 <- protein_groups[proteins_2]
    ref_proteins <- rownames(gold_standard)
    overlap_1 <- map_lgl(groups_1, ~ any(. %in% ref_proteins))
    overlap_2 <- map_lgl(groups_2, ~ any(. %in% ref_proteins))
    lab_idxs <- overlap_1 & overlap_2
    if (sum(lab_idxs) == 0)
      stop("no protein groups overlap between input and gold standard")
    # create labels
    labels <- rep(NA, nrow(dat))
    # map groups of length 1 by matrix indexing 
    lengths_1 <- lengths(groups_1)
    lengths_2 <- lengths(groups_2)
    single_lab_idxs <- lab_idxs & lengths_1 == 1 & lengths_2 == 1
    idxing_mat <- cbind(unlist(groups_1[single_lab_idxs]), 
                        unlist(groups_2[single_lab_idxs]))
    labels[single_lab_idxs] <- gold_standard[idxing_mat]
    # map groups of length > 1
    multi_lab_idxs <- lab_idxs & !single_lab_idxs
    labels[multi_lab_idxs] <- map_lgl(which(multi_lab_idxs), ~ any(
      gold_standard[as.matrix(crossing(groups_1[[.]], groups_2[[.]]))]))
  } 
  return(labels)
}