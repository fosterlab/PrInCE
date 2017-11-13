#' Predicting Interactomes from Co-Elution
#' 
#' <description>
#' 
PrInCE <- function(profile_matrix, gold_standard) {
  # check input 
  profile_matrix <- data.matrix(profile_matrix)
  if (!is.numeric(profile_matrix))
    stop("Input could not be coerced to numeric matrix")
  
  # build Gaussians
  gaussians <- build_gaussians(profile_matrix)

  # alignment
  
  # predict interactions
  interactions <- predict_interactions(profile_matrix, gaussians, gold_standard)

  ## also: fold-changes, complexes
}