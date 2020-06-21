#' Calculate the default features used to predict interactions in PrInCE
#' 
#' Calculate the six features that are used to discriminate interacting and
#' non-interacting protein pairs based on co-elution profiles in PrInCE, 
#' namely: raw Pearson R value, cleaned Pearson R value, raw Pearson P-value,
#' Euclidean distance, co-peak, and co-apex. Optionally, one or more of these
#' can be disabled.
#' 
#' @param profile_matrix a numeric matrix of co-elution profiles, with proteins
#' in rows, or a \code{\linkS4class{MSnSet}} object
#' @param gaussians a list of Gaussian mixture models fit to the profile matrix
#' by \code{link{build_gaussians}}
#' @param min_pairs minimum number of overlapping fractions between any given
#'   protein pair to consider a potential interaction
#' @param pearson_R_raw if true, include the Pearson correlation (R) between
#' raw profiles as a feature
#' @param pearson_R_cleaned if true, include the Pearson correlation (R) between
#' cleaned profiles as a feature
#' @param pearson_P if true, include the P-value of the Pearson correlation 
#' between raw profiles as a feature
#' @param euclidean_distance if true, include the Euclidean distance between
#' cleaned profiles as a feature
#' @param co_peak if true, include the 'co-peak score' (that is, the distance,
#' in fractions, between the single highest value of each profile) as a feature 
#' @param co_apex if true, include the 'co-apex score' (that is, the minimum
#' Euclidean distance between any pair of fit Gaussians) as a feature
#' 
#' @return a data frame containing the calculated features for all possible
#' protein pairs
#' 
#' @importFrom stats cor dist
#' @importFrom MSnbase exprs
#' @importFrom methods is
#' @importFrom Hmisc rcorr
#' @importFrom purrr map map_int
#' @importFrom dplyr filter
#' 
#' @export
calculate_features <- function(profile_matrix, gaussians,
                               min_pairs = 0,
                               pearson_R_raw = TRUE,
                               pearson_R_cleaned = TRUE,
                               pearson_P = TRUE,
                               euclidean_distance = TRUE,
                               co_peak = TRUE,
                               co_apex = TRUE,
                               n_pairs = FALSE
) {
  if (is(profile_matrix, "MSnSet")) {
    profile_matrix <- exprs(profile_matrix)
  }
  # remove at least pairs supported by only 2 points
  if (min_pairs < 3) {
    min_pairs <- 3
  }
  
  # replace missing values with near-zero noise
  cleaned <- clean_profiles(profile_matrix, impute_NA = FALSE, smooth = FALSE,
                            noise_floor = 0.05)
  proteins <- rownames(cleaned)
  n_proteins <- length(proteins)
  
  # create features list 
  feature_matrices <- list() 
  ## cor_R_raw, cor_R_cleaned, cor_P, eucl, co_peak, CA)
  
  # calculate Pearson correlation and P-value distances
  pairs <- crossprod(t(!is.na(profile_matrix)))
  if (pearson_R_raw) {
    cor_R_raw <- suppressWarnings(
      1 - cor(t(profile_matrix), use = 'pairwise.complete.obs'))
    feature_matrices[["cor_R_raw"]] <- cor_R_raw
  }
  if (pearson_R_cleaned) {
    cor_R_cleaned <- 1 - cor(t(cleaned))
    feature_matrices[["cor_R_cleaned"]] <- cor_R_cleaned
  }
  if (pearson_P) {
    cor_P <- suppressWarnings(
      rcorr(t(profile_matrix))$P)
    feature_matrices[["cor_P"]] <- cor_P  
  }
  # calculate Euclidean distance
  if (euclidean_distance) {
    eucl <- as.matrix(dist(cleaned, method = 'euclidean'))
    feature_matrices[["euclidean_distance"]] <- eucl
  }
  # calculate co-peak distance
  if (co_peak) {
    maxes <- apply(cleaned, 1, which.max)
    co_peak <- as.matrix(dist(maxes))
    feature_matrices[["co_peak"]] <- co_peak
  }
  # calculate co-apex (Gaussian) score
  if (co_apex) {
    CA <- co_apex(gaussians, proteins)
    feature_matrices[["co_apex"]] <- CA
  }
  
  # add number of pairs to the feature matrix list
  feature_matrices[["n_pairs"]] <- pairs
  
  ## make sure all dimensions are identical
  if (!all(map_int(feature_matrices, nrow) == n_proteins) |
      !all(map_int(feature_matrices, ncol) == n_proteins))
    stop("at least one feature matrix did not have correct dimensions")
  if (length(feature_matrices) == 0)
    stop("no features were calculated")
  first <- feature_matrices[[1]]
  tri <- upper.tri(first)
  idxs <- which(tri, arr.ind = TRUE)
  dat <- data.frame(protein_A = rownames(first)[idxs[, 1]], 
                    protein_B = rownames(first)[idxs[, 2]],
                    stringsAsFactors = FALSE) 
  dat <- cbind(dat, map(feature_matrices, ~ .[tri]))
  
  ## cap Euclidean distance at the 99.5th %ile
  max_euclidean_quantile = 0.995
  if ("euclidean_distance" %in% colnames(dat)) {
    vec = dat$euclidean_distance
    threshold = quantile(vec, probs = max_euclidean_quantile, na.rm = T)
    dat$euclidean_distance[vec >= threshold] <- threshold
  }
  
  # filter features based on n_pairs
  dat <- filter(dat, n_pairs >= min_pairs)
  # remove n_pairs from the feature list
  if (!n_pairs) {
    dat$n_pairs <- NULL
  }
  
  return(dat)
}
