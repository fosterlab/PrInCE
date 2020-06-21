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
#' function. The matrix, or matrices, can be provided to PrInCE either as
#' numeric matrices or as \code{\linkS4class{MSnSet}} objects.
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
#'   with proteins in rows and fractions in columns, or a list of matrices.
#'   Alternatively, can be provided as a single 
#'   \code{\linkS4class{MSnSet}} object or a list of objects.
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
#' @param precision optionally, return only interactions above the given 
#'   precision; by default, all interactions are returned and the user can
#'   subsequently threshold the list using the  
#'   \code{\link[PrInCE]{threshold_precision}} function
#' @param verbose if \code{TRUE}, print a series of messages about the stage
#'   of the analysis
#' @param min_points filter profiles without at least this many total, 
#'   non-missing points; passed to \code{\link{filter_profiles}}
#' @param min_consecutive filter profiles without at least this many 
#'   consecutive, non-missing points; passed to \code{\link{filter_profiles}}
#' @param min_pairs minimum number of overlapping fractions between any given
#'   protein pair to consider a potential interaction
#' @param impute_NA if true, impute single missing values with the average of
#'   neighboring values; passed to \code{\link{clean_profiles}}
#' @param smooth if true, smooth the chromatogram with a moving average filter;
#'   passed to \code{\link{clean_profiles}}
#' @param smooth_width width of the moving average filter, in fractions;
#'   passed to \code{\link{clean_profiles}}
#' @param max_gaussians the maximum number of Gaussians to fit; defaults to 5.
#'   Note that Gaussian mixtures with more parameters than observed (i.e., 
#'   non-zero or NA) points will not be fit. Passed to  
#'   \code{\link{choose_gaussians}}
#' @param criterion the criterion to use for model selection;
#'   one of "AICc" (corrected AIC, and default), "AIC", or "BIC". Passed to
#'   \code{\link{choose_gaussians}}
#' @param max_iterations the number of times to try fitting the curve with
#'   different initial conditions; defaults to 50. Passed to 
#'   \code{\link{fit_gaussians}}
#' @param min_R_squared the minimum R-squared value to accept when fitting the
#'   curve with different initial conditions; defaults to 0.5. Passed to 
#'   \code{\link{fit_gaussians}}
#' @param method the method used to select the initial conditions for
#'   nonlinear least squares optimization (one of "guess" or "random"); 
#'   see \code{\link{make_initial_conditions}} for details. Passed to 
#'   \code{\link{fit_gaussians}}
#' @param pearson_R_raw if true, include the Pearson correlation (R) between
#'   raw profiles as a feature
#' @param pearson_R_cleaned if true, include the Pearson correlation (R) between
#'   cleaned profiles as a feature
#' @param pearson_P if true, include the P-value of the Pearson correlation 
#'   between raw profiles as a feature
#' @param euclidean_distance if true, include the Euclidean distance between
#'   cleaned profiles as a feature
#' @param co_peak if true, include the 'co-peak score' (that is, the distance,
#'   in fractions, between the single highest value of each profile) as a 
#'   feature 
#' @param co_apex if true, include the 'co-apex score' (that is, the minimum
#'   Euclidean distance between any pair of fit Gaussians) as a feature
#' @param n_pairs if \code{TRUE}, include the number of fractions in which
#'   both of a given pair of proteins were detected as a feature
#' @param classifier the type of classifier to use: one of \code{"NB"} (naive
#'   Bayes), \code{"SVM"} (support vector machine), \code{"RF"} (random forest),
#'   \code{"LR"} (logistic regression), or \code{"ensemble"} (an ensemble of
#'   all four)
#' @param models the number of classifiers to train and average across, each
#'   with a different k-fold cross-validation split
#' @param cv_folds the number of folds to use for k-fold cross-validation
#' @param trees for random forests only, the number of trees in the forest
#' 
#' @return a ranked data frame of interacting proteins, with the precision
#'   at each point in the list
#'   
#' @examples
#' data(scott)
#' data(scott_gaussians)
#' data(gold_standard)
#' # analyze only the first 100 profiles
#' subset <- scott[seq_len(500), ]
#' gauss <- scott_gaussians[names(scott_gaussians) %in% rownames(subset)]
#' ppi <- PrInCE(subset, gold_standard, gaussians = gauss, models = 1, 
#'               cv_folds = 3)
#' 
#' @importFrom Rdpack reprompt
#' @importFrom MSnbase exprs
#' @importFrom methods is
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
PrInCE = function(profiles, gold_standard, 
                  gaussians = NULL, 
                  precision = NULL,
                  verbose = FALSE,
                  ## build_gaussians
                  min_points = 1,
                  min_consecutive = 5,
                  min_pairs = 3,
                  impute_NA = TRUE,
                  smooth = TRUE,
                  smooth_width = 4,
                  max_gaussians = 5,
                  max_iterations = 50,
                  min_R_squared = 0.5,
                  method = c("guess", "random"),
                  criterion = c("AICc", "AIC", "BIC"),
                  ## calculate_features
                  pearson_R_raw = TRUE,
                  pearson_R_cleaned = TRUE,
                  pearson_P = TRUE,
                  euclidean_distance = TRUE,
                  co_peak = TRUE,
                  co_apex = TRUE,
                  n_pairs = FALSE,
                  ## predict_interactions
                  classifier = c("NB", "SVM", "RF", "LR", "ensemble"), 
                  models = 10, 
                  cv_folds = 10,
                  trees = 500
) {
  method <- match.arg(method)
  criterion <- match.arg(criterion)
  classifie <- match.arg(classifier)
  
  # check profile input 
  if (is.list(profiles) & !is.data.frame(profiles)) {
    for (replicate_idx in seq_along(profiles)) {
      replicate <- profiles[[replicate_idx]]
      if (is(replicate, "MSnSet")) {
        profile_matrix <- exprs(replicate)
      } else {
        profile_matrix <- data.matrix(replicate)
      }
      if (!is.numeric(profile_matrix)) {
        stop("list input (item #", replicate_idx, 
             ") could not be coerced to numeric matrix")
      }
      profiles[[replicate_idx]] <- profile_matrix
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
    if (is(profiles, "MSnSet")) {
      profile_matrix <- exprs(profiles)
    } else {
      profile_matrix <- data.matrix(profiles)
    }
    if (!is.numeric(profile_matrix))
      stop("input could not be coerced to numeric matrix")
    # wrap in a list
    profiles <- list(profile_matrix)
    # also check Gaussians
    if (!is.null(gaussians)) {
      check_gaussians(gaussians)
      gaussians <- list(gaussians)
    }
  }
  
  # check gold standard input
  if (is.data.frame(gold_standard)) {
    # convert to adjacency matrix
    gold_standard <- adjacency_matrix_from_data_frame(gold_standard)
  } else if (is.list(gold_standard)) {
    gold_standard <- adjacency_matrix_from_list(gold_standard)
  }
  if (!is_unweighted(gold_standard)) {
    stop("could not convert supplied gold standards to adjacency matrix")
  }
  
  # get features for each matrix separately 
  features <- list()
  for (replicate_idx in seq_along(profiles)) {
    if (verbose) {
      message("generating features for replicate ", replicate_idx, " ...")
    }
    mat <- profiles[[replicate_idx]]
    
    # read Gaussians, or fit if they haven't been yet
    if (!is.null(gaussians)) {
      gauss <- gaussians[[replicate_idx]]
    } else {
      if (verbose) {
        message("  fitting Gaussians ...")
      }
      gauss <- build_gaussians(mat,
                               min_points = min_points,
                               min_consecutive = min_consecutive,
                               impute_NA = impute_NA,
                               smooth = smooth,
                               smooth_width = smooth_width,
                               max_gaussians = max_gaussians,
                               max_iterations = max_iterations,
                               min_R_squared = min_R_squared,
                               method = method)
    }
    
    # filter matrix based on Gaussians
    before <- nrow(mat)
    mat <- mat[names(gauss), ]
    after <- nrow(mat)
    if (verbose) {
      message("  fit mixtures of Gaussians to ", after, " of ", before, 
              " profiles")
    }
    
    # calculate features
    feat <- calculate_features(mat, gauss,
                               min_pairs = min_pairs,
                               pearson_R_raw = pearson_R_raw,
                               pearson_R_cleaned = pearson_R_cleaned,
                               pearson_P = pearson_P,
                               euclidean_distance = euclidean_distance,
                               co_peak = co_peak,
                               co_apex = co_apex,
                               n_pairs = n_pairs)
    features[[replicate_idx]] <- feat
  }
  
  # collapse into a single data frame 
  if (verbose) {
    message("concatenating features across replicates ...")
  }
  input <- concatenate_features(features)
  
  # predict interactions
  interactions <- predict_interactions(input, gold_standard, 
                                       classifier = classifier,
                                       models = models,
                                       cv_folds = cv_folds,
                                       trees = trees,
                                       verbose = verbose)
  
  # optionally threshold based on precison
  if (!is.null(precision)) {
    before <- nrow(interactions)
    interactions <- threshold_precision(interactions, precision)
    if (nrow(interactions) == 0) {
      warning("none of ", before, " ranked protein pairs had precision >= ", 
              precision)
    }
  }
  
  return(interactions)
}
