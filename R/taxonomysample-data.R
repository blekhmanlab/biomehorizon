#' Taxonomy information for OTUs
#'
#' Taxonomy information for each OTU listed in the sample OTU table. The first
#' variable contains OTU IDs, as listed in \code{otusample}, and subsequent
#' columns provide taxonomic clasification up to Genus, or the most specific
#' level possible for a given OTU. OTUs that are classified more broadly have
#' \code{NA} values for narrower taxonomic levels that do not apply.
#'
#' @docType data
#'
#' @usage data(taxonomysample)
#'
#' @format A data frame with 8814 rows and 7 variables.
#'
#' @keywords datasets
"taxonomysample"
