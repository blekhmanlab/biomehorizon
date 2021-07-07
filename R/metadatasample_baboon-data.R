#' Metadata on samples from OTU table of wild baboon data
#'
#' Metadata providing subject names and collection dates of samples from the
#' OTU table.
#'
#' @docType data
#'
#' @usage data(metadatasample_baboon)
#'
#' @format A data frame with 276 rows and 7 variables:
#' \describe{
#'   \item{subject}{subject name}
#'   \item{sample}{sample ID corresponding to a variable name from \code{otusample_baboon}}
#'   \item{collection_date}{number of days into the study the sample was collected}
#'   \item{sex}{if subject was male or female}
#'   \item{season}{if sample was collected in the wet or dry season}
#'   \item{rain_month_mm}{amount of rainfall in mm for month prior to sample collection}
#'   \item{diet_PC1}{a measure of dietary composition for a sample}
#' }
#'
#' @keywords datasets
"metadatasample_baboon"
