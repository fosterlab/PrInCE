#' PrInCE: Prediction of Interactomes from Co-Elution
#' 
#' PrInCE is a pipeline 
#' 
#' @param profile_matrix a numeric matrix of co-elution profiles, or something
#' that can be coerced to it, with proteins in rows
#' @param gold_standard an adjacency matrix of gold-standard interactions.
#' Zeroes in the adjacency matrix are interpreted by PrInCE as "true negatives"
#' when calculating precision. 
#' 
#' @return a ranked data frame of interacting proteins, with the precision
#' at each point in the list
#' 
#' @references 
#' \insertRef{stacey2017}{PrInCE}
#' 
#' \insertRef{scott2015}{PrInCE}
#' 
#' \insertRef{kristensen2012}{PrInCE}
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