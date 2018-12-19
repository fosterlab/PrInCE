#' Impute single missing values
#' 
#' Impute single missing values within a chromatogram profile as the average
#' of their neighbors.
#'
#' @param chromatogram a numeric vector corresponding to the chromatogram trace
#' 
#' @return the imputed chromatogram 
#' 
#' @examples
#' data(scott)
#' chrom = scott[16, ]
#' imputed = impute_neighbors(chrom)
#' 
#' @export
impute_neighbors <- function(chromatogram) {
  fractions <- length(chromatogram)
  nas <- which(is.na(chromatogram) | chromatogram == 0)
  for (i in nas) {
    prevIdx <- i - 1
    nextIdx <- i + 1
    if (prevIdx < 1 | nextIdx > fractions - 1)
      next
    neighbors <- c(chromatogram[prevIdx], chromatogram[nextIdx])
    if (!any(is.na(neighbors))) {
      chromatogram[i] <- mean(neighbors)
    }
  }
  
  return(chromatogram)
}
