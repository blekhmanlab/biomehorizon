#' Taxonomy information for OTUs
#'
#' Taxonomy information for each OTU ID listed in the sample OTU table. The first
#' variable contains OTU IDs, as listed in \code{otusample_baboon}, and subsequent
#' columns provide taxonomic clasification up to Genus, or the most specific
#' level possible for a given taxon. Taxo IDs that are classified more broadly have
#' \code{NA} values for narrower taxonomic levels that do not apply.
#'
#' @docType data
#'
#' @usage data(taxonomysample_baboon)
#'
#' @format A data frame with 2922 rows and 7 variables.
#'
#' @keywords datasets
"taxonomysample_baboon"
