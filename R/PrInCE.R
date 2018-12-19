#' PrInCE: Prediction of Interactomes from Co-Elution
#' 
#' PrInCE is a computational approach to infer protein-protein interaction 
#' networks from co-elution proteomics data, also called co-migration,
#' co-fractionation, or protein correlation profiling. This family of methods
#' separates interacting protein complexes on the basis of their diameter
#' or biochemical properties. Protein-protein interactions can then be 
#' inferred for pairs of proteins with similar elution profiles. PrInCE
#' implements a machine-learning approach to identify protein-protein 
#' interactions given a set of labelled examples, using features derived
#' exclusively from the data. This allows PrInCE to infer high-quality 
#' protein interaction networks from raw proteomics data, without bias
#' towards known interactions or functionally associated proteins, making 
#' PrInCE a unique resource for discovery.
#' 
#' PrInCE takes as input a co-elution matrix, with detected proteins in rows and
#' fractions as columns, and a set of 'gold standard' true positives and true
#' negatives. If replicate experiments were performed, a list of co-elution
#' matrices can be provided as input. PrInCE will construct features for each
#' replicate separately and use features from all replicates as input to the
#' classifier. The 'gold standard' can be either a data frame or adjacency
#' matrix of known interactions (and non-interactions), or a list of protein
#' complexes. For computational convenience, Gaussian mixture models can be
#' pre-fit to every profile and provided separately to the \code{PrInCE} 
#' function.
#' 
#' PrInCE implements three different types of classifiers to predict 
#' protein-protein interaction networks, including naive Bayes (the default),
#' random forests, and support vector machines. The classifiers are trained
#' on the gold standards using a ten-fold cross-validation procedure, training
#' on 90% of the data and predicting on the remaining 10%. For protein pairs
#' that are part of the training data, the held-out split is used to assign
#' a classifier score, whereas for the remaining protein pairs, the median of
#' all ten folds is used. Furthermore, to ensure the results are not sensitive
#' to the precise classifier split used, an ensemble of multiple classifiers
#' (ten, by default) is trained, and the classifier score is subsequently
#' averaged across classifiers. PrInCE can also ensemble across a set of
#' classifiers. 
#' 
#' By default, PrInCE calculates six features from each pair of co-elution
#' profiles as input to the classifier, including conventional similarity 
#' metrics but also several features specifically adapted to co-elution 
#' proteomics. For example, one such feature is derived
#' from fitting a Gaussian mixture model to each elution profile, then 
#' calculating the smallest Euclidean distance between any pair of fitted 
#' Gaussians. The complete set of features includes: 
#' 
#' \enumerate{
#'   \item the Pearson correlation between raw co-elution profiles;
#'   \item the p-value of the Pearson correlation between raw co-elution 
#'     profiles;
#'   \item the Pearson correlation between cleaned profiles, which are generated
#'     by imputing single missing values with the mean of their neighbors,
#'     replacing remaining missing values with random near-zero noise, and 
#'     smoothing the profiles using a moving average filter (see
#'     \code{\link[PrInCE]{clean_profile}});
#'  \item the Euclidean distance between cleaned profiles;
#'  \item the 'co-peak' score, defined as the distance, in fractions, between 
#'    the maximum values of each profile; and
#'  \item the 'co-apex' score, defined as the minimum Euclidean distance between
#'    any pair of fit Gaussians
#' }
#' 
#' The output of PrInCE is a ranked data frame, containing the classifier score
#' for every possible protein pair. PrInCE also calculates the precision at
#' every point in this ranked list, using the 'gold standard' set of protein
#' complexes or binary interactions. Our recommendation is to select a 
#' threshold for the precision and use this to construct an unweighted 
#' protein interaction network. 
#' 
#' @param profiles the co-elution profile matrix, or a list of profile matrices
#'   if replicate experiments were performed. Can be a single numeric matrix,
#'   with proteins in rows and fractions in columns, or a list of matrices  
#' @param gold_standard a set of 'gold standard' interactions, used to train the
#'   classifier. Can be provided either as an adjacency matrix, in which 
#'   both rows and columns correspond to protein IDs in the co-elution matrix
#'   or matrices, or as a list of proteins in the same complex, which will be
#'   converted to an adjacency matrix by PrInCE. Zeroes in the adjacency matrix 
#'   are interpreted by PrInCE as "true negatives" when calculating precision. 
#' @param gaussians optionally, provide Gaussian mixture models fit by 
#'   the \code{\link[PrInCE]{build_gaussians}} function. If \code{profiles} is
#'   a numeric matrix, this should be the named list output by 
#'   \code{\link[PrInCE]{build_gaussians}} for that matrix; if \code{profiles} 
#'   is a list of numeric matrices, this should be a list of named lists
#' @param verbose if \code{TRUE}, print a series of messages about the stage
#'   of the analysis
#' 
#' @return a ranked data frame of interacting proteins, with the precision
#'   at each point in the list
#' 
#' @export
#' 
#' @references 
#' \insertRef{stacey2017}{PrInCE}
#' 
#' \insertRef{scott2015}{PrInCE}
#' 
#' \insertRef{kristensen2012}{PrInCE}
#' 
#' \insertRef{skinnider2018}{PrInCE}
PrInCE = function(profiles, gold_standard, gaussians = NULL, verbose = F) {
  # check profile input 
  if (is.list(profiles)) {
    for (replicate_idx in seq_along(profiles)) {
      profile_matrix = profiles[[replicate_idx]]
      profile_matrix = data.matrix(profile_matrix)
      profiles[[replicate_idx]] = profile_matrix
      if (!is.numeric(profile_matrix))
        stop("list input (item #", replicate_idx, 
             ") could not be coerced to numeric matrix")
      # also check Gaussians
      if (!is.null(gaussians)) {
        if (length(gaussians) < replicate_idx) {
          stop("fewer Gaussians than profiles provided")
        }
        check_gaussians(gaussians[[replicate_idx]], rownames(profile_matrix),
                        replicate_idx)
      }
    }
  } else {
    profile_matrix = data.matrix(profiles)
    if (!is.numeric(profile_matrix))
      stop("input could not be coerced to numeric matrix")
    # wrap in a list
    profiles = list(profile_matrix)
    # also check Gaussians
    if (!is.null(gaussians)) {
      check_gaussians(gaussians)
      gaussians = list(gaussians)
    }
  }
  
  # check gold standard input
  if (is.data.frame(gold_standard)) {
    # convert to adjacency matrix
    gold_standard = PrInCE::adjacency_matrix_from_data_frame(gold_standard)
  } else if (is.list(gold_standard)) {
    gold_standard = PrInCE::adjacency_matrix_from_list(gold_standard)
  }
  if (!PrInCE::is_unweighted(gold_standard)) {
    stop("could not convert supplied gold standards to adjacency matrix")
  }

  # get features for each matrix separately 
  features = list()
  for (replicate_idx in seq_along(profiles)) {
    if (verbose) {
      message("generating features for replicate ", replicate_idx, " ...")
    }
    mat = profiles[[replicate_idx]]
    
    # read Gaussians, or fit if they haven't been yet
    if (!is.null(gaussians)) {
      gauss = gaussians[[replicate_idx]]
    } else {
      if (verbose) {
        message("  fitting Gaussians ...")
        gauss = build_gaussians(mat, max_iterations = 50)
      }
    }
    
    # filter matrix based on Gaussians
    before = nrow(mat)
    mat = mat[names(gauss), ]
    after = nrow(mat)
    if (verbose) {
      message("  fit mixtures of Gaussians to ", after, " of ", before, 
              " profiles")
    }
    
    # calculate features
    feat = calculate_features(mat, gauss)
    features[[replicate_idx]] = feat
  }
  
  # collapse into a single data frame 
  if (verbose) {
    message("concatenating features across replicates ...")
  }
  input = concatenate_features(features)
  
  # predict interactions
  interactions = predict_interactions(input, gold_standard, verbose = verbose)
  
  # return
  return(interactions)
}
