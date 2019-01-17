#' Filter a co-elution profile matrix
#' 
#' Filter a matrix of co-elution/co-fractionation profiles by removing
#' profiles without a certain number of non-mising or consecutive points.
#' 
#' @param profile_matrix a numeric matrix of co-elution profiles, with proteins
#'   in rows, or a \code{\linkS4class{MSnSet}} object
#' @param min_points filter profiles without at least this many total, 
#'   non-missing points
#' @param min_consecutive filter profiles without at least this many 
#'   consecutive, non-missing points
#' 
#' @return the filtered profile matrix
#' 
#' @examples
#' data(scott)
#' nrow(scott)
#' filtered = filter_profiles(scott)
#' nrow(scott)
#' 
#' @importFrom MSnbase exprs
#' @importFrom Biobase exprs<-
#' @importFrom methods is
#' 
#' @export
filter_profiles <- function(profile_matrix, min_points = 1, 
                            min_consecutive = 5) {
  if (is(profile_matrix, "MSnSet")) {
    msn = profile_matrix
    profile_matrix = exprs(profile_matrix)
  }

  # filter profiles without N non-missing points
  nas <- is.na(profile_matrix)
  if (!is.na(min_points) & !is.null(min_points) & min_points > 0) {
    profile_matrix <- profile_matrix[rowSums(!nas) >= min_points,]
  } else {
    # need to at least filter profiles without any points
    profile_matrix <- profile_matrix[rowSums(!nas) >= 1,]
  }
  
  # filter profiles without N consecutive points after imputation
  if (!is.na(min_consecutive) & !is.null(min_consecutive) &
      min_consecutive > 0) {
    imputed <- t(apply(profile_matrix, 1, impute_neighbors))
    imputed_nas <- is.na(imputed)
    rles <- apply(imputed_nas, 1, rle)
    max_consecutive <- purrr::map_int(rles, ~ max(.$lengths[.$values == FALSE]))
    profile_matrix <- profile_matrix[max_consecutive >= min_consecutive,]
  }
  
  if (exists("msn")) {
    exprs(msn) = profile_matrix
    return(msn)
  } else {
    return(profile_matrix)
  }
}
