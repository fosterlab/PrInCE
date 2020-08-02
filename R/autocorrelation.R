#' Calculate autocorrelation from a pair co-elution profiles
#'
#' Calculates a protein's chromatogram to all other eluted proteins. The
#' resulting vector of correlation values is correlated with the same protein
#' in the co-elution of the second experimental condition. This gives a metric
#' whereby higher values reflect proteins that are not rewired between conditions,
#' but lower values reflect proteins that re more rewired between conditions.
#'
#' @param profile1 a matrix or dataframe, the first condition with proteins in
#' rows and fractions in columns.
#'
#' @param profile2 a matrix or dataframe, the second condition
#' with proteins in rows and fractions in columns.
#'
#' @param cor_method the correlation method to use. (One of "Pearson",
#' "Spearman", "Kendall")
#'
#' @param min_fractions the minimum number of fractions a protein must appear
#' in to be included in autocorrelation
#'
#' @param min_pairs ignores pairs of chromatograms within a profile that don't
#' co-occur in at least this many fractions.
#'
#' @return a named vector of autocorrelation scores for all proteins found in
#' both matricies
#'
#' @importFrom dplyr na_if
#' @importFrom purr map
#' @importFrom purr map_dbl
#'
#' @export
autocorrelation <- function(profile1,
                            profile2,
                            cor_method = c("pearson", "spearman", "kendall"),
                            min_fractions = 1,
                            min_pairs = 10) {
  cor_method <- match.arg(cor_method)

  # check profile inputs
  if (!(is.matrix(profile1) && is.matrix(profile2))) {
    if (is.data.frame((profile1) && is.data.frame(profile2))) {
      profile1 <- as.matrix(profile1)
      profile2 <- as.matrix(profile2)
    } else {
      stop("Input pair must be matricies or dataframes")
    }
  }

  # convert zeros to NA
  profile1 <- na_if(profile1, 0)
  profile2 <- na_if(profile2, 0)

  # remove proteins that appear in less than min_fractions
  fracnum1 <- apply(profile1, 1, function(x) {
    sum(is.finite(x))
  }) >= min_fractions
  profile1 <- profile1[fracnum1, ]

  fracnum2 <- apply(profile2, 1, function(x) {
    sum(is.finite(x))
  }) >= min_fractions
  profile2 <- profile2[fracnum2, ]

  # filter to proteins in both conditions
  overlap <- intersect(rownames(profile1), rownames(profile2))
  profile1 <- profile1[overlap, ]
  profile2 <- profile2[overlap, ]

  # calculate correlations for each chromatogram
  chroms <- list(profile1, profile2)
  cors <- map(chroms, ~ cor(
    t(.),
    method = cor_method,
    use = "pairwise.complete.obs"
  ))
  ns <- map(chroms, ~ crossprod(!is.na(t(.))))

  # censor correlations with < min_pairs in the raw chromatograms
  cors[[1]][ns[[1]] < min_pairs] <- NA
  cors[[2]][ns[[2]] < min_pairs] <- NA

  # calculate autocorrelations
  autocors <- map_dbl(seq_along(overlap), ~ cor(
    cors[[1]][., ], cors[[2]][., ],
    use = "pairwise.complete.obs"
  ))
  names(autocors) <- overlap
}