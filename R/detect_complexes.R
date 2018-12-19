#' Detect significantly interacting complexes in a chromatogram matrix
#' 
#' Use a permutation testing approach to identify complexes that show a 
#' significant tendency to interact, relative to random sets of complexes
#' of equivalent size. The function begins by calculating the Pearson 
#' correlation or Euclidean distance between all proteins in the matrix, and 
#' 
#' @param mat a matrix of chromatograms, with proteins in the columns and
#'   fractions in the rows
#' @param complexes a named list of protein complexes, where the name is the
#'   complex name and the entries are proteins within that complex 
#' @param min_pairs the minimum number of pairwise observations to count a 
#'   correlation or distance towards the z score 
#' @param method method to use to calculate edge weights;
#'   one of \code{pearson} or \code{euclidean}
#' @param bootstraps number of bootstraps to execute to estimate z scores
#' @param progress whether to show the progress of the function
#' 
#' @return a named vector of z scores for each complex in the input list
#' 
#' @importFrom stats cor dist na.omit median sd
#' @importFrom utils combn
#' @importFrom progress progress_bar
#' @importFrom purrr map_dbl
#' @export
detect_complexes = function(mat, complexes, 
                            method = c("pearson", "euclidean"),
                            min_pairs = 10, 
                            bootstraps = 100,
                            progress = T) {
  method = match.arg(method)
  
  # construct network
  n = crossprod(!is.na(mat))
  if (method == "pearson") {
    network = cor(mat, use = 'pairwise.complete.obs')
  } else if (method == "euclidean") {
    network = as.matrix(dist(t(mat)))  
  }
  network[n < min_pairs] = NA
  
  # do permutation testing
  z_scores = numeric(0)
  pb = progress::progress_bar$new(
    format = "complex :what [:bar] :percent eta: :eta",
    clear = F, total = length(complexes), width = 80)
  for (i in seq_along(complexes)) {
    set.seed(i)
    complex_name = names(complexes)[i]
    complex = complexes[[i]]
    
    # subset complex to proteins present in this network 
    nodes = colnames(network)
    overlap = intersect(complex, nodes)
    
    # abort if overlap size is < 3
    if (length(overlap) < 3) {
      z_scores[complex_name] = NA
    } else {
      # calculate median PCC for intra-complex interactions
      idxing_mat = t(combn(overlap, 2))
      edge_weights = na.omit(network[idxing_mat])
      obs = median(edge_weights, na.rm = T)
      
      # compare to random expectation
      rnd = purrr::map_dbl(seq_len(bootstraps), ~ {
        # generate random set of proteins equivalent to # of complex subunits
        # present in this network 
        rnd_complex = sample(nodes, length(overlap))
        # calculate observed D statistic
        idxing_mat = t(combn(rnd_complex, 2))
        rnd_weights = as.numeric(na.omit(network[idxing_mat]))
        median(rnd_weights, na.rm = T)
      })
      
      # calculate Z score
      z = (obs - mean(rnd, na.rm = T)) / sd(rnd, na.rm = T)
      z_scores[complex_name] = z
    }
    
    # tick progress bar
    if (progress) {
      pb$tick(tokens = list(what = sprintf(
        paste0("%-", nchar(length(complexes)), "s"), i)))
    }
  }
  
  return(z_scores)
}