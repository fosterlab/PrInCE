#' Cytoplasmic interactome of Jurkat T cells during apoptosis
#' 
#' Co-migration profiles derived from size exclusion chromatography (SEC) of
#' cytoplasmic fractions from Jurkat T cells, 4 hours following Fas stimulation.
#' 
#' Protein quantitation was accomplished by SILAC (stable isotopic labelling
#' by amino acids in cell culture), and is ratiometric, i.e., it reflects 
#' the ratio between the intensity of the heavy isotope and the light isotope 
#' ("H/L"). The dataset was initially described in Scott et al., 
#' \emph{Mol. Syst. Biol.} 2017. The heavy isotope channel from replicate
#' 1 is included in the \code{PrInCE} package. The R script used to generate
#' this matrix from the supplementary materials of the paper is provided in the
#' \code{data-raw} directory of the package source code.
#'
#' @docType data
#' @usage data(scott)
#' @format a data frame with 1560 rows and 55 columns, with proteins in rows and
#'   SEC fractions in columns
#' @source \url{http://msb.embopress.org/content/13/1/906}
"scott"