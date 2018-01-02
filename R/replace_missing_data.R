#' Replace missing data with median ± random noise 
#' 
#' Replace missing data within each numeric column of a data frame with 
#' the column median, plus or minus some random noise, in order to train 
#' classifiers that do not easily ignore missing data (e.g. random forests or
#' support vector machines).
#' 
#' @param input the data frame to 
#' @param noise_pct the standard deviation of the random normal 
#' distribution from which to draw added noise, expressed as a 
#' percentage of the standard deviation of the non-missing values in each 
#' column
#' 
#' @return a data frame with missing values in each numeric column replaced
#' by the column median, plus or minus some random noise
#' 
#' @export
replace_missing_data <- function(input, noise_pct = 0.05) {
  for (col_name in colnames(input)) {
    column <- input[[col_name]]
    if (!is.numeric(column))
      next
    missing <- is.na(column)
    input[[col_name]][missing] <- 
      median(column, na.rm = T) + rnorm(sum(missing)) * 
      sd(column, na.rm = T) * noise_pct 
  }
  return(input)
}