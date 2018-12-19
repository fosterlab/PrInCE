#' Make labels for a classifier based on a gold standard
#' 
#' Create labels for a classifier for protein pairs in the same order as 
#' in a dataset that will be used as input to a classifier, in a 
#' memory-friendly way. 
#' 
#' @param gold_standard an adjacency matrix of gold-standard interactions
#' @param input a data frame with interacting proteins in the first two
#' columns 
#' @param protein_groups optionally, specify a list linking each protein in the 
#' first two columns of the input data frame to a protein group 
#' 
#' @return a vector of the same length as the input dataset, containing 
#' \code{NA}s for protein pairs not in the gold standard and ones or zeroes
#' based on the content of the adjacency matrix
#' 
#' @examples 
#' data(gold_standard)
#' adj = adjacency_matrix_from_list(gold_standard)
#' proteins = unique(unlist(gold_standard))
#' input = data.frame(protein_A = sample(proteins, 10), 
#'                    protein_B = sample(proteins, 10))
#' labels = make_labels(adj, input)
#' 
#' @export
make_labels <- function(gold_standard, input, protein_groups = NULL) {
  proteins_1 <- as.character(input[[1]])
  proteins_2 <- as.character(input[[2]])
  if (is.null(protein_groups)) {
    lab_idxs <- proteins_1 %in% rownames(gold_standard) &
      proteins_2 %in% rownames(gold_standard)
    if (sum(lab_idxs) == 0)
      stop("no proteins overlap between input and gold standard")
    idxing_mat <- cbind(proteins_1[lab_idxs], proteins_2[lab_idxs])
    labels <- rep(NA, nrow(input))
    labels[lab_idxs] <- gold_standard[idxing_mat]
  } else {
    # filter groups based on presence/absence in matrix
    protein_groups <- purrr::map(protein_groups,
                                 ~ .[. %in% rownames(gold_standard)])
    # identify rows that overlap with the gold standard matrix
    groups_1 <- protein_groups[proteins_1]
    groups_2 <- protein_groups[proteins_2]
    ref_proteins <- rownames(gold_standard)
    overlap_1 <- purrr::map_lgl(groups_1, ~ any(. %in% ref_proteins))
    overlap_2 <- purrr::map_lgl(groups_2, ~ any(. %in% ref_proteins))
    lab_idxs <- overlap_1 & overlap_2
    if (sum(lab_idxs) == 0)
      stop("no protein groups overlap between input and gold standard")
    # create labels
    labels <- rep(NA, nrow(input))
    # map groups of length 1 by matrix indexing 
    lengths_1 <- lengths(groups_1)
    lengths_2 <- lengths(groups_2)
    single_lab_idxs <- lab_idxs & lengths_1 == 1 & lengths_2 == 1
    idxing_mat <- cbind(unlist(groups_1[single_lab_idxs]), 
                        unlist(groups_2[single_lab_idxs]))
    labels[single_lab_idxs] <- gold_standard[idxing_mat]
    # map groups of length > 1
    multi_lab_idxs <- lab_idxs & !single_lab_idxs
    labels[multi_lab_idxs] <- purrr::map_lgl(which(multi_lab_idxs), ~ any(
      gold_standard[as.matrix(tidyr::crossing(groups_1[[.]], groups_2[[.]]))]))
  } 
  return(labels)
}