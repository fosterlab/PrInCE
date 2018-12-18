#' Fitted Gaussian mixture models for the \code{scott} dataset
#' 
#' The \code{\link[PrInCE]{scott}} dataset consists of protein
#' co-migration profiles derived from size exclusion chromatography (SEC) of
#' cytoplasmic fractions from Jurkat T cells, 4 hours following Fas stimulation.
#' The \code{scott_gaussians} object contains Gaussian mixture models fit by
#' the function \code{\link[PrInCE]{build_gaussians}}; this is bundled with the
#' R package in order to expedite the demonstration code, as the process of 
#' Gaussian fitting is one of the more time-consuming aspects of the package. 
#' 
#' As with the \code{\link[PrInCE]{scott}} dataset, the code used to generate
#' this data object is provided in the \code{data-raw} directory of the package
#' source.
#'
#' @docType data
#' @usage data(scott_gaussians)
#' @format a named list with 970 entries; names are proteins, and list items 
#'   conain information about fitted Gaussians in the format that PrInCE expects
"scott_gaussians"