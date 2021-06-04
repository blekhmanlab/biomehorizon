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
#'   \item{subject}{subject name in character format}
#'   \item{sample}{sample ID in character format corresponding to a variable name from \code{otusample_diet}}
#'   \item{collection_date}{number of days in numeric format into the study the sample was collected}
#'   \item{supplement}{metadata variable in character format indicating if subject was given EVOO or MCT as a dietary supplement on days 10-17}
#' }
#'
#' @keywords datasets
"metadatasample_diet"
