#' Taxonomy information for OTUs
#'
#' Taxonomy information for each taxon ID listed in the sample OTU table. The first
#' variable contains taxon IDs, as listed in \code{otusample_diet}, and the subsequent
#' column provides taxonomic clasification up to Genus, or the most specific
#' level possible for a given taxon. Taxo IDs that are classified more broadly have
#' \code{NA} values for narrower taxonomic levels that do not apply.
#'
#' @docType data
#'
#' @usage data(taxonomysample_diet)
#'
#' @format A data frame with 4583 rows and 2 variables.
#'
#' @keywords datasets
"taxonomysample_diet"
