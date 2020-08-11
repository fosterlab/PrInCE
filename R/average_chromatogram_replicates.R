#' Average values of chromatogram matrix over several replicates
#'
#' Find all proteins that are detected between all replicates and generates
#' the mean value of each fraction of each protein that occur in at least
#' a certain number of replicates, adding \code{NA}s otherwise.
#'
#' @param mats a list of chromatogram matricies
#' @param min_replicates the minimum number of replicates the protein fraction
#' must appear in.
#'
#' @return a matrix of averaged chromatogram values between replicates
#'
#' @importFrom dplyr na_if
#' @importFrom purrr map
#' @importFrom purrr reduce
#' @importFrom magrittr %>%
#'
#' @export
average_chromatogram_replicates <- function(mats, min_replicates = 1) {
  mats <- map(mats, ~ na_if(., 0))

  # get all proteins throughout replicates
  proteins <- map(mats, ~ rownames(.)) %>%
    reduce(union)

  # create an empty matrix with all proteins and columns
  emptymatrix <- matrix(
    0,
    nrow = length(proteins),
    ncol = ncol(mats[[1]]),
    dimnames = list(proteins, colnames(mats[[1]]))
  )

  # count in union matrix if rep matrix has value
  countchroms <- function(mat) {
    mat1 <- emptymatrix
    mat1[rownames(mat), colnames(mat)] <- !is.na(mat)
    return(mat1)
  }

  # write in union matrix if rep matrix's value
  sumchroms <- function(mat) {
    mat1 <- emptymatrix
    mat1[rownames(mat), colnames(mat)] <- (mat)
    mat1[is.na(mat1)] <- 0
    return(mat1)
  }

  num_replicates <- map(mats, ~ countchroms(.)) %>%
    reduce(`+`)

  sum_replicates <- map(mats, ~ sumchroms(.)) %>%
    reduce(`+`)

  num_replicates[num_replicates < min_replicates] <- NA

  mean <- sum_replicates / num_replicates
  return(mean)
}