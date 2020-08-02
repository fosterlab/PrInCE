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
#' @export
autocorrelation <- function(profile1,
                            profile2,
                            cor_method = c("pearson", "spearman", "kendall"),
                            min_fractions = 1,
                            min_pairs = 0) {

}