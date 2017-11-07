#' Predicting Interactomes from Co-Elution
#' 
#' <description>
#' 
PrInCE <- function(profile_matrix) {
  # check input 
  profile_matrix <- data.matrix(profile_matrix)
  if (!is.numeric(profile_matrix))
    stop("Input could not be coerced to numeric matrix")
  
  # build Gaussians
  gaussians <- build_gaussians(profile_matrix)

  # 
}