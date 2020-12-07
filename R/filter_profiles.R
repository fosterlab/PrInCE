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
#' filtered <- filter_profiles(scott)
#' nrow(scott)
#' 
#' @importFrom MSnbase exprs
#' @importFrom Biobase exprs<-
#' @importFrom methods is
#' @importFrom purrr map_int
#' 
#' @export
filter_profiles <- function(profile_matrix, min_points = 5, 
                            min_consecutive = 1) {
  # extract the matrix and impute individual missing points
  if (is(profile_matrix, "MSnSet")) {
    expr <- exprs(profile_matrix)
  } else {
    expr <- profile_matrix
  }
  imputed <- t(apply(expr, 1, impute_neighbors))
  
  # filter profiles without N non-missing points after imputation
  nas <- is.na(imputed)
  if (!is.na(min_points) & !is.null(min_points) & min_points > 0) {
    profile_matrix <- profile_matrix[rowSums(!nas) >= min_points, ]
  } else {
    # need to at least filter profiles without any points
    profile_matrix <- profile_matrix[rowSums(!nas) >= 1, ]
  }
  
  # filter profiles without N consecutive points after imputation
  if (!is.na(min_consecutive) & !is.null(min_consecutive) &
      min_consecutive > 0) {
    rles <- apply(nas, 1, rle)
    max_consecutive <- map_int(rles, ~ ifelse(FALSE %in% .$values, max(.$lengths[.$values == FALSE]), 0))
    profile_matrix <- profile_matrix[max_consecutive >= min_consecutive, ]
  }
  
  return(profile_matrix)
}
