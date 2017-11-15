#' Preprocess a co-elution profile
#' 
#' Clean a co-elution/co-fractionation profile by (1) imputing single missing
#' values with the average of neighboring values, (2) replacing missing values
#' with random, near-zero noise, and (3) smoothing with a moving average
#' filter. 
#' 
#' @param chromatogram a numeric vector corresponding to the chromatogram trace
#' @param impute_NA if true, impute single missing values with the average of
#' neighboring values 
#' @param smooth if true, smooth the chromatogram with a moving average filter
#' @param smooth_width width of the moving average filter, in fractions
#' @param noise_floor mean value of the near-zero noise to add  
#' 
#' @return a cleaned profile
#' 
#' @export
clean_profile <- function(chromatogram, impute_NA = T, smooth = T,
                          smooth_width = 4, noise_floor = 0.001) {
  # copy chromatogram
  cleaned <- chromatogram
    
  # impute missing values using mean of neighbors
  if (impute_NA) {
    impute_neighbors(cleaned)
  }
  
  # replace remaining NAs with near-zero noise
  nas <- is.na(cleaned) | cleaned == 0
  cleaned[nas] <- runif(sum(nas), min = 0, max = noise_floor)
  
  # smooth with a moving average filter
  if (smooth) {
    smooth_input <- c(runif(smooth_width, min = 0, max = noise_floor),
                      cleaned, 
                      runif(smooth_width, min = 0, max = noise_floor))
    smoothed <- forecast::ma(smooth_input, smooth_width)
    smoothed <- smoothed[seq(smooth_width + 1, length(smoothed) - smooth_width)]
    cleaned <- setNames(smoothed, names(cleaned))
  }
  
  return(cleaned)
}
