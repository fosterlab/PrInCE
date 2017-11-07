#' Preprocess a co-elution profile matrix
#' 
#' Clean a matrix of co-elution/co-fractionation profiles by 
#' (1) imputing single missing
#' values with the average of neighboring values, (2) replacing missing values
#' with random, near-zero noise, and (3) smoothing with a moving average
#' filter. 
#' 
#' @param profile_matrix a numeric matrix of co-elution profiles, with proteins
#' in rows
#' @param impute_NA if true, impute single missing values with the average of
#' neighboring values 
#' @param smooth if true, smooth the chromatogram with a moving average filter
#' @param smooth_width width of the moving average filter, in fractions 
#' 
#' @return a cleaned matrix
#' 
#' @export
clean_profiles <- function(profile_matrix, impute_NA = T, smooth = T, 
                         smooth_width = 4) {
  if (impute_NA | smooth)
    message(".. cleaning ", nrow(profile_matrix), " profiles")
  t(apply(profile_matrix, 1, clean_profile, impute_NA, smooth, smooth_width))
}
