#' Metadata on samples from OTU (or other lowest taxonomic level) table of human diet data
#'
#' Metadata providing subject names and collection dates of samples from the
#' OTU table.
#'
#' @docType data
#'
#' @usage data(metadatasample_diet)
#'
#' @format A data frame with 483 rows and 4 variables:
#' \describe{
#'   \item{subject}{subject name}
#'   \item{sample}{sample ID corresponding to a variable name from \code{otusample_diet}}
#'   \item{collection_date}{number of days into the study the sample was collected}
#'   \item{Supplement}{metadata variable indicating if subject was in the EVOO or MCT group}
#' }
#'
#' @keywords datasets
"metadatasample_diet"
