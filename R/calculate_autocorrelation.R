#' Calculate the autocorrelation for each protein between a pair of co-elution 
#' experiments.
#'
#' For a given protein, the correlation coefficient to all other proteins 
#' in the first condition is calculated, yielding a vector of correlation 
#' coefficients. The same procedure is repeated for the second condition,
#' and the two vectors of correlation coefficients are themselves correlated,
#' yielding a metric whereby higher values reflect proteins with unchanging
#' interaction profiles between conditions, while lower values reflect proteins 
#' with substantially changing interaction profiles.
#'
#' Note that all of zero, \code{NA}, \code{NaN}, and infinite values are
#' all treated equivalently as missing values when applying the 
#' \code{min_fractions} and \code{min_pairs} filters, but different handling of 
#' missing values will produce different autocorrelation scores.
#'
#' @param profile1 a numeric matrix or data frame with proteins in rows and 
#'   fractions in columns, or a \code{\linkS4class{MSnSet}} object, representing
#'   the first co-elution condition
#' @param profile2 a numeric matrix or data frame with proteins in rows and 
#'   fractions in columns, or a \code{\linkS4class{MSnSet}} object, representing
#'   the second co-elution condition
#' @param cor_method the correlation method to use; one of \code{"pearson"},
#'   \code{"spearman"}, or \code{"kendall"}).
#' @param min_fractions filter proteins not quantified in at least this many
#'   fractions
#' @param min_pairs remove correlations between protein pairs not co-occuring in
#'   at least this many fractions from the autocorrelation calculation
#'
#' @return a named vector of autocorrelation scores for all proteins found in
#' both matrices.
#'
#' @importFrom purrr map_dbl
#'
#' @export
calculate_autocorrelation <- function(profile1,
                                      profile2,
                                      cor_method = c("pearson", "spearman", 
                                                     "kendall"),
                                      min_replicates = 1,
                                      min_fractions = 1,
                                      min_pairs = 0) {
  cor_method <- match.arg(cor_method)
  
  # coerce profiles to matrices
  profile1 <- coerce_to_matrix(profile1)
  profile2 <- coerce_to_matrix(profile2)
  
  # identify 'missing' values as zeroes, NAs, or NaN
  missing1 <- !is.finite(profile1) | profile1 == 0
  missing2 <- !is.finite(profile2) | profile2 == 0
  
  # remove proteins that appear in less than min_fractions
  n_fractions1 <- rowSums(!missing1)
  n_fractions2 <- rowSums(!missing2)
  profile1 <- profile1[n_fractions1 >= min_fractions, ]
  profile2 <- profile2[n_fractions2 >= min_fractions, ]
  
  # filter to proteins found in both conditions
  overlap <- intersect(rownames(profile1), rownames(profile2))
  profile1 <- profile1[overlap, ]
  profile2 <- profile2[overlap, ]
  ## also filter missing
  missing1 <- missing1[overlap, ]
  missing2 <- missing2[overlap, ]
  
  # calculate correlations for each chromatogram
  cor1 <- cor(t(profile1), method = cor_method, use = 'pairwise.complete.obs')
  cor2 <- cor(t(profile2), method = cor_method, use = 'pairwise.complete.obs')
  
  # censor correlations with less than min_pairs in the raw chromatograms
  n1 <- crossprod(t(!missing1))
  n2 <- crossprod(t(!missing2))
  cor1[n1 < min_pairs] <- NA
  cor2[n2 < min_pairs] <- NA
  
  # calculate autocorrelations
  autocors <- map_dbl(overlap, ~ cor(cor1[., ], cor2[., ], 
                                     use = 'pairwise.complete.obs'))
  names(autocors) <- overlap
  return(autocors)
}

#' @importFrom MSnbase exprs
coerce_to_matrix <- function(matrix_like) {
  profile_matrix <- matrix_like
  if (is(matrix_like, "MSnSet")) {
    profile_matrix <- exprs(matrix_like)
  } else {
    profile_matrix <- data.matrix(matrix_like)
  }
  if (!is.numeric(profile_matrix)) {
    stop("input could not be coerced to numeric matrix: ",
         paste(dim(profile_matrix), collapse = ' '))
  }
  return(profile_matrix)
}
