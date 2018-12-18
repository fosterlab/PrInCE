#' Reference set of human protein complexes
#' 
#' A reference set of 477 experimentally confirmed human protein complexes,
#' derived from the EBI Complex Portal database. 
#' 
#' 477 protein complexes, ranging in size from 2 to 44 proteins and involving
#' 877 proteins in total, to provide a reference set of true positive and true 
#' negative interactions (intra- and inter-complex interactions, respectively)
#' for demonstration in PrInCE analysis of a co-elution dataset.
#' Other "gold standards" are possible in practice, most notably the CORUM
#' database; however, the Complex Portal reference set is included in this
#' package due to its CC-BY license. 
#' 
#'
#' @docType data
#' @usage data(gold_standard)
#' @format a list containing 477 entries (character vectors)
#' @source \url{https://www.ebi.ac.uk/complexportal/complex/organisms}
"gold_standard"