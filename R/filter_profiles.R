#' Filter a co-elution profile matrix
#' 
#' Filter a matrix of co-elution/co-fractionation profiles by removing
#' profiles without a certain number of non-mising or consecutive points.
#' 
#' @param profile_matrix a numeric matrix of co-elution profiles, with proteins
#' in rows
#' @param min_points filter profiles without at least this many total, 
#' non-missing points
#' @param min_consecutive filter profiles without at least this many 
#' consecutive, non-missing points
#' 
#' @export
filter_profiles <- function(profile_matrix, min_points = 5, 
                            min_consecutive = 5) {
  nas <- is.na(profile_matrix)

  # filter profiles without a certain number of points
  if (!is.na(min_points) & !is.null(min_points) & min_points > 0) {
    profile_matrix <- profile_matrix[rowSums(!nas) >= min_points,]
  }
  
  # filter profiles without a certain number of consecutive points 
  if (!is.na(min_consecutive) & !is.null(min_consecutive) &
      min_consecutive > 0) {
    imputed <- t(apply(profile_matrix, 1, impute_neighbors))
    imputed_nas <- is.na(imputed)
    profile_matrix <- profile_matrix[rowSums(!imputed_nas) >= min_consecutive,]
  }
  
  return(profile_matrix)
}