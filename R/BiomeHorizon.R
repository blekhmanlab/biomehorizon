#' Preliminary Data Cleaning and Preperation
#'
#' This function prepares the OTU table and additional datasets for analysis
#' with the \code{horizonplot()} function.
#'
#' The \code{prepanel()} function has 6 main purposes in preparing data sets and
#' other parameters for the main \code{horizonplot()} function:
#'
#' 1) Filter the OTU table to the OTUs displayed on the final horizon plot, and
#' to the samples of just one individual (for datasets with multiple subjects).
#' By default, the "most important" OTUs are selected using four filtering
#' thresholds: \code{thresh_prevalence}, \code{thresh_abundance},
#' \code{thresh_abundance_override}, and \code{thresh_NA}. They can also be
#' manually specified as a vector of OTU IDs using \code{otulist}.
#'
#' 2) If single OTU analysis is enabled, convert the OTU table to values by
#' subject for the OTU being analyzed
#'
#' 3) Ensure data sets are formatted correctly
#'
#' 4) Set the functions for finding the \code{origin} and horizon band thickness
#' (\code{band.thickness}) of each OTU panel, if the default (\code{NA}) or a
#' constant is entered.
#'
#' 5) Set other parameters to their defaults, and ensure correct data types are
#' entered. For boolean values, \code{NA} is converted to \code{FALSE}.
#'
#' 6) Check for common user errors, such as entering ".8" rather than "80" as a
#' percentage filtering threshold (this will leave a warning message).
#'
#' By default, OTUs are filtered automatically using two thresholds. An
#' abundance threshold (\code{thresh_abundance}) sets the minimum average
#' proportion an OTU must represent across all samples, and a prevalence
#' threshold (\code{thresh_prevalence}) sets the minimum proportion of all
#' samples where this OTU must be present (at least 1 sample read). These
#' thresholds can be used in combination, or alone by setting one of them to
#' \code{0} or \code{NA}.
#'
#' In addition, you can set a second abundance threshold that overrides the
#' prevalence threshold if it is reached, using
#' \code{thresh_abundance_override}. This is useful for catching OTUs that are
#' abundant for a brief period of time, but are absent from most of the samples,
#' and are nevertheless important to include in analysis. This is disabled by
#' default (\code{thresh_abundance_override == NA}).
#'
#' Finally, a fourth filtering threshold, \code{thresh_NA}, filters out OTUs
#' with missing data in a substantial fraction of the samples. This defaults to
#' eliminating OTUs missing data in >5\% of samples.
#'
#' Alternatively, OTUs can be manually specified in \code{otulist} as a vector
#' of OTU IDs. The order in which these are specified will also determine the
#' arrangement of OTU panels on the horizon plot.
#'
#' You can also compare a single OTU across multiple subjects, by specifying the
#' OTU ID in \code{singleVarOTU}. This is useful for comparing the same
#' timepoint across multiple individuals, rather than multiple OTUs or taxa.
#'
#' @param otudata Data frame representing OTU Table. Assumes first column
#'   contains OTU IDs, and all other columns are numeric vectors containing the
#'   number of sample reads for each OTU. Values can also be represented as
#'   proportions or percentages of the total sample for each OTU.
#' @param metadata Data frame representing metadata table; matches samples to
#'   collection dates, and to subject names if applicable. If this data frame is
#'   supplemented, the columns with sample IDs, collection dates and subject
#'   names should be named "sample", "collection_date" and "subject",
#'   respectively. \code{collection_date} must be of class \code{numeric} or
#'   \code{Date}.
#' @param taxonomydata Taxonomy information for OTUs, used for labeling facets.
#'   There are two options: \itemize{ \item A data frame with columns as
#'   taxonomic levels, plus a column for OTU IDs. Assumes first column contains
#'   OTU IDs, and all other columns are character vectors. If OTUs have
#'   different levels of taxonomic classification (e.g. one is specified up to
#'   Genus and one only to Phylum), then NAs should substitute levels without
#'   specification. \item A vector where each element contains the entire
#'   taxonomy for an OTU, with taxonomic levels separated by semicolons. The
#'   order of this vector should match the order of OTU IDs specified in
#'   \code{otudata}. } Taxonomic levels should start from Kingdom and can go as
#'   far as Subspecies. Defaults to \code{NA} (do not label by taxonomy).
#' @param thresh_prevalence numeric threshold for OTU filtering. Minimum \% of
#'   total samples in which OTU must be present to be included in analysis
#'   (defaults to 80).
#' @param thresh_abundance numeric threshold for OTU filtering. Minimum \% of
#'   total sample reads the OTU must constitute to be included in analysis
#'   (defaults to 0.5).
#' @param thresh_abundance_override numeric threshold for OTU filtering. Minimum
#' \% of total sample reads the OTU must constitute to override all other
#' standards, and be included in analysis (defaults to \code{NA}: disabled).
#' @param thresh_NA numeric threshold for OTU filtering. Maximum \% of samples
#'   with missing data (defaults to 5).
#' @param otulist character vector specifying OTU IDs for manual selection. Also
#'   determines the order from top to bottom of OTU panels displayed on the
#'   horizon plot. Defaults to \code{NA} (use filtering thresholds). In this case, OTU
#'   panels will be ordered alphabetically by OTU ID.
#' @param regularInterval integer. For regularized data, this specifies the
#'   fixed interval of days separating each sample timepoint. If this value is
#'   20, for example, new timepoints will be created at 1, 21, 41, 61, etc. To
#'   leave data irregularly spaced, do not specify a number here. Defaults to
#'   \code{NA} (do not regularize).
#' @param maxGap numeric specifying the maximum number of days between the
#'   previous and subsequent irregular timepoints in order to interpolate a new
#'   timepoint. If the distance between the nearest time points exceeds the
#'   threshold specified by \code{maxGap}, all OTU values for that time point
#'   will be set to \code{NA}, and a scale break in the time axis will appear on the
#'   horizon plot. Must be an integer > 0.
#' @param minSamplesPerFacet numeric. For regularized data with breaks in the
#'   time axis, specifies the minimum number of samples required of each facet
#'   time interval. Facets without this many timepoints will be removed.
#'   Defaults to 2.
#' @param subj character, used for datasets with multiple individual
#'   microbiomes. Filter samples to this subject or subjects. In most cases, you
#'   should specify just one subject, but if single OTU analysis is enabled you
#'   can select multiple subjects. Subject names should be described in metadata
#'   under the variable "subject". Defaults to \code{NA} (assume all samples are
#'   from one individual; do not filter by subject name).
#' @param singleVarOTU character string specifying an OTU ID for facetting by
#'   subject. Facetting by subject requires metadata with columns on sample and
#'   subject, with an equal number of samples for each subject. If collection
#'   dates are provided, they must be identical for each subject. If they are
#'   not provided, the function assumes samples are ordered chronologically. A
#'   subset of subjects may be selected for analysis by supplying a vector of
#'   multiple subjects to \code{subj}.
#' @param band.thickness The height of each horizontal band (denoted by a unique
#'   color), i.e. the size of the scale of a horizon subplot. There are three
#'   options: \itemize{ \item If \code{NA}, the default, the band thickness will
#'   be evaluated using the function \code{function(y) {max((abs(y -
#'   origin(y))), na.rm=TRUE) / nbands}}. This calculates the maximum extreme
#'   (lowest or highest abundance value) divided by the number of bands. \item A
#'   \code{function} will be called with a single argument, the sample values
#'   for one OTU, to evaluate a unique band thickness for each panel based on
#'   its sample values. The return value must be numeric. \item A numeric
#'   constant, providing a fixed band thickness for all OTUs. This should be
#'   expressed as a percentage (0-100).}
#' @param origin The baseline (value=0, the base of the first positive band) for
#'   horizon subplots. There are three options: \itemize{ \item If \code{NA},
#'   the default, the origin will be evaluated separately for each OTU using the
#'   median of the sample values. \item A \code{function} will be called with a
#'   single argument, a numeric vector representing the sample values for one
#'   OTU, to evaluate a unique origin for each panel. The return value must be
#'   numeric. \item A numeric constant, providing a fixed origin value for all
#'   OTUs. This should be expressed as a percentage (0-100).}
#' @param facetLabelsByTaxonomy If \code{TRUE}, label facets by taxonomy, using
#'   \code{taxonomydata}. Facets will be labelled using the most specific
#'   classification available for each OTU. If \code{FALSE} (default), label
#'   facets by OTU ID.
#' @param customFacetLabels Use a custom character vector to label facets. Length
#'   of this vector should match the number of OTUs post-filtering, or the number of
#'   subjects if single OTU analysis is enabled. Overrides
#'   facetLabelsByTaxonomy, but if set to \code{NA} (the default),
#'   facetLabelsByTaxonomy is used instead.
#' @param interpolate_NA logical. How should \code{NA} values be dealt with? If
#'   \code{TRUE} (default), \code{NA} values are interpolated using previous and
#'   subsequent OTU values. If \code{FALSE}, they are set to value=0. Note that
#'   this only applies to sample timepoints that contain values for some OTUs;
#'   if a sample consists entirely of NAs, it will be treated as a break in the
#'   timescale (see \code{maxGap}).
#' @param formatStep If \code{FALSE} (default), horizon plot is a line graph. If
#'   \code{TRUE}, horizon plot is formatted as a step graph, with steps
#'   horizontal and then vertical.
#' @param nbands integer specifying the number of positive bands (each denoted
#'   by a unique color) on each horizon subplot. For example, if you set
#'   \code{nbands=4}, there will be four positive bands and four negative bands,
#'   with 8 total colors. Must be an integer >=3. If \code{nbands} > 5, you must
#'   supply your own color palette of length \code{2 * nbands}.
#'
#' @return Returns a list containing the appropriate arguments for the
#'   horizonplot function. This result list should then be inputted into
#'   \code{horizonplot()} to produce the graph. You should not need to alter any
#'   parameters in this list before using them in \code{horizonplot}, but this
#'   preliminary function allows you to check the refined parameters in case of
#'   an error in \code{horizonplot}.
#'
#' @examples
#' # Pass just the OTU table to prepanel, and it will assume all samples belong
#' # to the same subject.
#' prepanel(otusample = otusample_diet)
#'
#' # Supplement metadata and a subject name, and it will select samples from
#' # just one subject (this is what you should do with more than one subject).
#' prepanel(otusample = otusample_diet, metadatasample = metadatasample_diet, subj="MCTs01")
#'
#' # Pass taxonomydata to prepanel if you want to label facets by taxonomy
#' # rather than by OTU ID.
#' prepanel(otusample = otusample_diet, metadatasample = metadatasample_diet, 
#' taxonomydata = taxonomysample_diet, subj="MCTs01", facetLabelsByTaxonomy=TRUE)
#'
#' # OTU filtering using both a prevalence and an abundance standard (default)
#' prepanel(otusample = otusample_diet, metadatasample = metadatasample_diet, subj="MCTs01", 
#' thresh_prevalence=75, thresh_abundance=0.75)
#'
#' # OTU filtering using just an abundance standard
#' prepanel(otusample = otusample_diet, metadatasample = metadatasample_diet, subj="MCTs01",
#' thresh_prevalence=NA, thresh_abundance=0.75)
#'
#' # If an OTU's average abundance reaches a high enough threshold, override
#' # other standards and include it in analysis
#' prepanel(otusample = otusample_diet, metadatasample = metadatasample_diet, subj="MCTs01", 
#' thresh_prevalence=90, thresh_abundance=0.75, thresh_abundance_override=1.5)
#'
#' # Filter OTUs where >2% samples are NA values
#' prepanel(otusample = otusample_diet, metadatasample = metadatasample_diet, subj="MCTs01", 
#' thresh_NA=2)
#'
#' # You can also manually select OTUs by OTU ID
#' prepanel(otusample = otusample_diet, metadatasample = metadatasample_diet, subj="MCTs01",
#' otulist=c("taxon 1", "taxon 2", "taxon 10", "taxon 14"))
#'
#' # Manual selection can be used to specify the order OTUs will appear on
#' # the horizon plot. For example, these two datasets have identical OTUs, but
#' # they are ordered differently.
#' params <- prepanel(otusample = otusample_diet, metadatasample = metadatasample_diet, 
#' subj="MCTs01", thresh_prevalence=95, thresh_abundance=1.5, 
#' otulist=c("taxon 1", "taxon 2", "taxon 10", "taxon 14"))
#' params[[1]]$otuid
#' params <- prepanel(otusample = otusample_diet, metadatasample = metadatasample_diet, 
#' subj="MCTs01", otulist=c("taxon 10", "taxon 2", "taxon 1", "taxon 14"))
#' params[[1]]$otuid
#'
#' # The origin and band.thickness variables can be set to either a numeric
#' # constant or a function that evaluates separately for every OTU subpanel based
#' # on its sample values.
#'
#' # Use a fixed origin of 5% for all OTU subpanels
#' prepanel(otusample = otusample_diet, metadatasample = metadatasample_diet, 
#' subj="MCTs01", origin=5)
#'
#' # Evaluate a different origin for each OTU subpanel using a custom function
#' prepanel(otusample = otusample_diet, metadatasample = metadatasample_diet, 
#' subj="MCTs01", origin=function(y){mad(y, na.rm=TRUE)})
#'
#' @import dplyr
#' @importFrom magrittr %>%
#'
#' @export
prepanel <- function(otudata, metadata=NA, taxonomydata=NA,
                     thresh_prevalence=80, thresh_abundance=0.5, thresh_abundance_override=NA, thresh_NA=5,
                     regularInterval=NA, maxGap=NA, minSamplesPerFacet=2,
                     otulist=NA,
                     subj=NA, singleVarOTU=NA,
                     band.thickness=NA, origin=NA,
                     facetLabelsByTaxonomy=FALSE, customFacetLabels=NA,
                     interpolate_NA=TRUE, formatStep=FALSE,
                     nbands=4)
{
  # Basic otudata format checks
  if(!is.data.frame(otudata)) {
    stop("`otudata` must be a data frame")
  }
  if(ncol(otudata) < 2) {
    stop("`otudata` must have at least one sample")
  }
  if(nrow(otudata) < 1) {
    stop("`otudata` must have at least one OTU")
  }
  if(FALSE %in% sapply(otudata[,-1],function(x){any(class(x)==c("numeric","integer"))})) {
    stop("All columns of otudata except the first must be of class numeric or integer")
  }
  if(min(otudata[,-1]) < 0) {
    stop("`otudata` cannot have negative values")
  }
  colnames(otudata)[1] <- "otuid"
  if(length(unique(otudata$otuid)) < nrow(otudata)) {
    stop("`otudata` has duplicate OTU rows")
  }

  # Check for proportional otudata format
  if(all(otudata[,-1] %% 1 == 0)) {
    # Convert data to %
    otudata[,-1] <- 100 * as.data.frame(prop.table(as.matrix(otudata[,-1]),2))
  } else {
    if(max(otudata[,-1]) <= 1) {
      # values are proportions - convert to percentages
      if(!all(colSums(otudata[,-1])==1)) {
        stop("Column proportions in `otudata` don't add up to 1")
      }
      otudata[,-1] <- 100 * otudata[,-1]
    } else {
      if(max(otudata)[,-1] > 100) {
        stop("Percentage values in `otudata` are greater than 100")
      }
      if(!all(colSums(otudata[,-1])==100)) {
        stop("Column percentages in `otudata` don't add up to 100")
      }
    }
  }

  # Convert vector taxonomydata to data frame
  if(class(taxonomydata) %in% c("character","factor")) {
    if(length(taxonomydata) != nrow(otudata)) {
      stop("Length of `taxonomydata` does not match number of OTUs")
    }
    taxlist <- strsplit(as.character(taxonomydata),split=";")
    nTaxLevels <- max(sapply(taxlist,length))
    taxlist <- lapply(taxlist, function(x){c(x,rep(NA,nTaxLevels-length(x)))})
    taxonomydata <- taxlist %>% as.data.frame() %>% t() %>% as.data.frame()
    rownames(taxonomydata) <- c()
    taxonomydata$otuid <- otudata$otuid
    taxonomydata <- taxonomydata %>% dplyr::select(otuid,everything())
  }

  # Checks for single-OTU analysis
  if(!is.na(singleVarOTU)) {
    if(!(singleVarOTU %in% otudata$otuid)) {
      stop("singleVarOTU was not found in otudata")
    }
    if(length(singleVarOTU) > 1) {
      stop("singleVarOTU must specify only one OTU")
    }
    if(!(NA %in% subj) & length(subj)==1) {
      warning("Did you mean to perform single OTU analysis with just one subject selected? Set singleVarOTU=NA to display multiple OTUs.")
    }
    if(length(otulist) != 1 || !is.na(otulist)) {
      stop("Cannot perform single OTU analysis with OTUs manually selected - either set otulist=NA or singleVarOTU=NA")
    }
  }

  # Check filtering thresholds
  if(is.na(thresh_prevalence) | class(thresh_prevalence) != "numeric") {
    stop("thresh_prevalence must be of type numeric")
  }
  if(thresh_prevalence > 100 | thresh_prevalence < 0) {
    stop("thresh_prevalence must be a number between 0 and 100")
  }
  if(thresh_prevalence < 1 & thresh_prevalence > 0) {
    warning(paste("thresh_prevalence is less than 1%, and is supposed to be a number between 0 and 100. Did you mean thresh_prevalence=",100*thresh_prevalence,"?",sep=""))
  }
  if(is.na(thresh_abundance) | class(thresh_abundance) != "numeric") {
    stop("thresh_abundance must be of type numeric")
  }
  if(thresh_abundance > 100 | thresh_abundance < 0) {
    stop("thresh_abundance must be a number between 0 and 100")
  }
  if(!is.na(thresh_abundance_override)) {
    if(class(thresh_abundance_override) != "numeric") {
      stop("thresh_abundance_override must be of type numeric")
    }
    if(thresh_abundance_override<0 | thresh_abundance_override>100) {
      stop("thresh_abundance_override must be a number between 0 and 100 (or NA)")
    }
  }
  if(is.na(thresh_NA) | class(thresh_NA) != "numeric") {
    stop("thresh_NA must be of type numeric")
  }
  if(thresh_NA > 100 | thresh_NA < 0) {
    stop("thresh_NA must be a number between 0 and 100")
  }

  # Subj checks
  if(NA %in% subj) {
    subj <- NA
  }
  if(length(subj) > 1) {
    if(is.na(singleVarOTU)) {
      stop("Cannot select multiple subjects for filtering without performing single OTU analysis - select an OTU ID for `singleVarOTU`")
    }
    subjlist <- subj
    subj <- NA
  } else {
    subjlist <- NA
  }

  # Check format of metadata and filter by subject
  if(is.data.frame(metadata)) {
    if(!("sample" %in% colnames(metadata))) {
      stop("metadata must contain variable `sample`")
    }
    if(nrow(metadata) < 1) {
      stop("metadata must have at least one sample")
    }
    if(length(unique(metadata$sample)) != nrow(metadata)) {
      stop("metadata has duplicate sample rows")
    }
    if(!("subject" %in% colnames(metadata))) {
      if(!is.na(subj) | (length(subjlist)!=1 && !is.na(subjlist))) {
        stop("metadata must contain variable `subject` to filter by subject")
      }
      if(!is.na(singleVarOTU)) {
        stop("metadata must contain variable `subject` to perform single OTU analysis")
      }
    } else {
      if(is.na(subj)) {
        if(is.na(singleVarOTU) & length(unique(metadata$subject)) > 1) {
          stop("metadata has more than one subject; select a subject using the parameter `subj` to filter metadata")
        }
        if(!is.na(singleVarOTU) & length(unique(table(metadata$subject))) != 1) {
          stop("metadata must have the same number of samples for each subject to perform single OTU analysis")
        }
        metadata <- metadata %>% dplyr::arrange(subject)
      }
      # Filter subjects
      if(!is.na(subj)) {
        if(!(subj %in% metadata$subject)) {
          stop("`subj` was not found in subject column of metadata")
        }
        metadata <- metadata %>% dplyr::filter(subject == subj)
      }
      if((length(subjlist)!=1 && !is.na(subjlist))) {
        if(!all(subjlist %in% metadata$subject)) {
          stop("Some elements of `subj` were not found in subject column of metadata")
        }
        metadata <- metadata %>% dplyr::filter(subject %in% subjlist)
      }
    }
    if(!("collection_date" %in% colnames(metadata))) {
      if(!is.na(regularInterval) & regularInterval > 0) {
        stop("metadata must contain variable `collection_date` to regularize data")
      }
    } else {
      if(!any(class(metadata$collection_date) == c("Date","numeric","integer"))) {
        stop("variable `collection_date` in metadata must be of class Date, numeric or integer")
      }
      if(NA %in% metadata$collection_date) {
        stop("NAs found in metadata$collection_date")
      }
      if(!is.na(singleVarOTU)) {
        metadata <- metadata %>% dplyr::arrange(subject, as.numeric(collection_date))
        # Make sure sample collection dates are consistent between subjects
        datesBySubject <- lapply(metadata %>% dplyr::group_by(subject) %>% dplyr::group_split(), `[`,,"collection_date")
        if(!all(sapply(datesBySubject,function(x){x==datesBySubject[[1]]}))) {
          stop("Collection dates must be the same between subjects")
        }
      } else {
        metadata <- metadata %>% dplyr::arrange(as.numeric(collection_date))
      }
    }
  } else {
    if(!is.na(subj) | (length(subjlist)!=1 && !is.na(subjlist))) {
      stop("cannot filter otudata by subject without metadata on sample and subject")
    }
    if(!is.na(regularInterval) & regularInterval > 0) {
      stop("cannot regularize data without metadata containing collection dates")
    }
    if(!is.na(singleVarOTU)) {
      stop("cannot perform single OTU analysis without metadata on sample and subject")
    }
  }

  # Select and order samples in otudata, if metadata provided
  if(is.data.frame(metadata)) {
    if(sum(!(metadata$sample %in% colnames(otudata)[-1])) > 0) {
      stop("Some samples in metadata have no match in otudata")
    }
    if(is.na(singleVarOTU)) {
      otudata <- otudata %>% dplyr::select(otuid, as.character(metadata$sample))
    }
  }

  # Filter otudata
  if(is.na(singleVarOTU)) {
    if(length(otulist)==1 && is.na(otulist)) {
      otudata$nonzero <- rowSums(otudata[,-1] != 0, na.rm = TRUE)
      otudata$prevalence <- 100 * otudata$nonzero/(ncol(otudata)-2)
      otudata$abundance <- rowSums(otudata[,2:(ncol(otudata)-2)], na.rm=TRUE) / otudata$nonzero
      otudata$nacount <- rowSums(is.na(otudata[,2:(ncol(otudata)-3)]))
      otudata <- dplyr::filter(otudata, (100 * nacount / (ncol(otudata)-5) < thresh_NA) & ((prevalence > thresh_prevalence & abundance > thresh_abundance) | (!is.na(thresh_abundance_override) & abundance > thresh_abundance_override)))
      otu_stats <- otudata %>% dplyr::select(otuid, prevalence, abundance, nacount)
      otudata <- otudata %>% dplyr::select(-nonzero, -prevalence, -abundance, -nacount)
    } else {
      otudata <- otudata %>% dplyr::filter(otuid %in% otulist)
      if(nrow(otudata)==0) {
        stop("no values in otulist match OTU IDs in otudata; all rows filtered out")
      }
    }
  } else {
    otudata <- otudata %>% dplyr::filter(otuid == singleVarOTU)
  }

  # Covert otudata format for single variable analysis. Remove extra metadata columns
  if(!is.na(singleVarOTU)) {
    samplenames <- metadata %>% dplyr::select(subject, sample) %>% dplyr::group_by(subject) %>% dplyr::mutate(row_id=1:dplyr::n()) %>% dplyr::ungroup() %>% tidyr::spread(subject,sample) %>% dplyr::select(-row_id) %>% t() %>% as.data.frame() %>% tibble::rowid_to_column("otuid")
    ids <- samplenames$otuid
    samplenames <- samplenames %>% dplyr::select(-otuid)
    otudata <- matrix(otudata[1,][c(as.matrix(samplenames))],nrow(samplenames)) %>% as.data.frame() %>% dplyr::mutate(otuid=ids) %>% dplyr::select(otuid,everything())
    otudata$otuid <- unique(metadata$subject)
  }

  # Determine timestamps
  if(is.data.frame(metadata) & "collection_date" %in% colnames(metadata) & !is.na(regularInterval) & regularInterval>0) {
    if(is.na(singleVarOTU)) {
      dates <- as.numeric(metadata$collection_date)
    } else {
      dates <- as.numeric((metadata %>% dplyr::filter(subject==metadata$subject[1]))$collection_date)
    }
    timestamps <- 1 + dates - dates[1]
  } else {
    timestamps <- NA
  }

  # Check taxonomydata format
  if(is.data.frame(taxonomydata) & is.na(singleVarOTU)) {
    if(ncol(taxonomydata) > 9) {
      stop("taxonomydata has more than 9 columns; the function does not support classification past Subspecies")
    }
    taxonomynames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Subspecies")
    colnames(taxonomydata) <- c("otuid",taxonomynames[1:(ncol(taxonomydata)-1)])

    if(sum(!(otudata$otuid %in% taxonomydata$otuid)) > 0) {
      stop("Some OTU IDs in otudata have no match in taxonomydata")
    }
    taxonomydata <- taxonomydata %>% dplyr::filter(otuid %in% otudata$otuid) %>% dplyr::arrange(otudata$otuid)

    if(nrow(taxonomydata) != nrow(otudata)) {
      stop("number of OTU rows in taxonomydata does not match number of OTU rows in otudata")
    }
  } else {
    taxonomydata <- NA
  }

  # Set variables to correct values
  if(!is.na(regularInterval) & regularInterval<=0) {
    warning("regularInterval is <= 0, setting to NA")
    regularInterval <- NA
  }
  if(is.na(regularInterval)) {
    maxGap <- NA
  } else {
    if(!is.na(maxGap) & maxGap < regularInterval) {
      maxGap <- NA
      warning("maxGap is less than regularInterval; setting to default maxGap=NA")
    }
  }
  if(is.na(formatStep)) {
    formatStep <- FALSE
  } else {
    if(class(formatStep) != "logical") {
      stop("formatStep must be of type logical")
    }
  }
  if(is.na(interpolate_NA)) {
    interpolate_NA <- FALSE
  } else {
    if(class(interpolate_NA) != "logical") {
      stop("interpolate_NA must be of type logical")
    }
  }
  if(is.na(facetLabelsByTaxonomy)) {
    facetLabelsByTaxonomy <- FALSE
  } else {
    if(class(facetLabelsByTaxonomy) != "logical") {
      stop("facetLabelsByTaxonomy must be of type logical")
    }
    if(facetLabelsByTaxonomy & !is.data.frame(taxonomydata)) {
      stop("Need `taxonomydata` to label facets by taxonomy")
    }
  }
  if(length(customFacetLabels) > 1 && length(customFacetLabels) < nrow(otudata)) {
    warning("customFacetLabels has length less than the number of OTUs; using default labels.")
    customFacetLabels <- NA
  } else if(length(customFacetLabels) > nrow(otudata)) {
    warning("customFacetLabels has length greater than the number of OTUs, so the label vector was truncated.")
    customFacetLabels <- customFacetLabels[1:nrow(otudata)]
  }
  if(!(class(nbands) %in% c("numeric","integer"))) {
    stop("nbands must be of class numeric or integer")
  }
  if(nbands %% 1 != 0) {
    stop("nbands must be a whole number")
  }
  if(nbands < 3) {
    stop("nbands must be >=3")
  }

  # Set functions
  if(!is.function(origin) && is.na(origin)) {
    origin <- function(y) {
      median(y, na.rm=TRUE)
    }
  } else {
    if(is.numeric(origin)) {
      if(origin < 0 | origin > 100) {
        stop("Constant values of `origin` must be expressed as a percentage (0-100)")
      }
      originconstant <- origin
      origin <- function(y) {
        function(x) { y }
      }
      origin <- origin(originconstant)
    } else {
      if(!is.function(origin)) {
        warning("origin must be a function operating on a numeric vector. Using default function (median).")
        origin <- function(y) {
          median(y, na.rm=TRUE)
        }
      }
    }
  }

  if(!is.function(band.thickness) && is.na(band.thickness)) {
    band.thickness <- function(y) {
      max((abs(y - origin(y))), na.rm=TRUE) / nbands
    }
  } else {
    if(is.numeric(band.thickness)) {
      if(band.thickness < 0 | band.thickness > 100) {
        stop("Constant values of `band.thickness` must be expressed as a percentage (0-100)")
      }
      btconstant <- band.thickness
      band.thickness <- function(y) {
        function(x) { y }
      }
      band.thickness <- band.thickness(btconstant)
    } else {
      if(!is.function(band.thickness)) {
        warning("band.thickness must be a function operating on a numeric vector. Using default function.")
        band.thickness <- function(y) {
          max((abs(y - origin(y))), na.rm=TRUE) / nbands
        }
      } else {
        # If function is supplied, reset its environment
        environment(band.thickness) <- environment()
      }
    }
  }

  fill_NA <- ifelse(interpolate_NA, function(y){zoo::na.approx(y)}, function(y){zoo::na.fill(y,0)})

  if((length(subjlist)!=1 && !is.na(subjlist))) {
    subj <- subjlist
  }
  outputVars <- sapply(c("thresh_prevalence", "thresh_abundance", "thresh_abundance_override", "thresh_NA", "subj", "singleVarOTU"),function(x){get(x)})
  outputVars <- outputVars[!is.na(outputVars)]

  cat("Constructed an OTU table and other variables with the following settings:",
      paste(names(outputVars), ": ", outputVars, sep=""),
      sep="\n")

  if(length(otulist)==1 && is.na(otulist) & is.na(singleVarOTU)) {
    otu_stats <- dplyr::tibble("OTU_ID" = otu_stats$otuid,
                               "Average_abundance" = otu_stats$abundance,
                               "Prevalence" = otu_stats$prevalence,
                               "Num_missing_samples" = otu_stats$nacount)
    cat("\n", length(unique(otudata$otuid)), " OTUs met the filtering requirements, with the following stats:\n", sep="")
    print.data.frame(otu_stats)

    biomehorizonpkg_otu_stats <<- otu_stats
    cat("`biomehorizonpkg_otu_stats` was outputted to the environment.\n")
  }

  list(otudata,taxonomydata,timestamps,otulist,
       subj,
       regularInterval,maxGap,minSamplesPerFacet,
       band.thickness,origin,
       facetLabelsByTaxonomy,customFacetLabels,
       fill_NA,nbands,formatStep)
}

#' Construct a Microbiome Horizon Plot
#'
#' This is the main function that constructs and returns the microbiome horizon
#' plot.
#'
#' After data sets and other parameters have been properly formatted and checked
#' for errors in the \code{prepanel} function, they are entered into this
#' function, which constructs and returns the horizon plot. All customizations
#' of the graph should be specified in \code{prepanel()} and not here; no
#' alteratons should be made to the output list before it is entered into this
#' function.
#'
#' The refined version of \code{otudata} used in this function represents a
#' filtered OTU table, containing only the OTUs to be displayed on the graph,
#' and only the samples belonging to the subject selected. Sample values in this
#' refined table reflect difference in fractional abundance from the origin
#' value. Values are converted from raw sample reads to proportions of the
#' entire sample represented by a given OTU, and then the proportion values
#' within each OTU are centered to their respective \code{origin} values.
#'
#' The refined \code{taxonomydata} is filtered to just the OTUs in
#' \code{otudata}.
#'
#' @section Irregular Data:
#'
#' A common problem faced in visualizing time series data is plotting data
#'   spaced at irregular time intervals. A common solution for this problem is
#'   interpolating values at regular time intervals using nearby data. However,
#'   since microbiome data can change drastically in short periods of time, it
#'   doesn't make sense to interpolate through large timespans, and thus the
#'   function gives several options to plot irregular data as accurately as
#'   possible.
#'
#'   \enumerate{ \item Plot "real values" but with an inaccurate timescale. For
#'   this default option, samples will be plotted next to each other regardless
#'   of their timestamps. This is most accurate in that timepoints are plotted
#'   directly from sample values, but risks being misleading if the timescale is
#'   not clearly marked as inconsistent.  Additionally, this option removes the
#'   ability to visually compare temporal differences within the same plot.
#'   \item Plot artificial values but with a regularized timescale. New values
#'   can be interpolated using existing ones at a regular interval of time
#'   specified by \code{regularInterval}. The first sample is plotted as "day
#'   1," and a new value is interpolated using the closest previous and
#'   subsequent sample timepoints at a fixed interval throughout the rest of the
#'   data.  This "regularization" of the data allows for quick visual comparison
#'   of microbiome changes within the plot. The downside of this method is that
#'   "real" values are not plotted (except for rare cases where a sample
#'   timepoint happens to fall on the regular interval), and innaccuracies are
#'   created through interpolation. This is especially true given the
#'   continuous, rapid changes of bacterial abundances within the microbiome.
#'   \item Compromise between accuracy of values and a regular timescale:
#'   interpolation within clusters of closely-spaced data, which are separated
#'   by breaks in the time axis. This allows for temporal comparison within each
#'   cluster of timepoints and avoids interpolating across large timespans. This
#'   method is practical for datasets where samples are collected irregularly,
#'   arranged in periodic clusters of closely-spaced data separated by larger
#'   timespans with fewer samples. Clustering is done by specifying a value for
#'   \code{maxGap}, which defines the threshold of time without data to separate
#'   clusters.}
#'
#' @param parameterList The list of parameters for constructing the horizon
#'   plot. This should come directly from the output list of the
#'   \code{prepanel()} function, without alteration. The 15 parameters, in
#'   order, are \code{otudata}, \code{taxonomydata}, \code{timestamps},
#'   \code{otulist}, \code{subj}, \code{regularInterval}, \code{maxGap},
#'   \code{minSamplesPerFacet} \code{band.thickness}, \code{origin},
#'   \code{facetLabelsByTaxonomy}, \code{customFacetLabels}, \code{fill_NA},
#'   \code{nbands}, and \code{formatStep}.
#'
#'   Two of the parameters are not arguments of the \code{prepanel()} function,
#'   and are described below:
#'
#'   \code{timestamps} is an integer vector containing the time (days) each
#'   sample was collected, retrieved from the \code{collection_date} variable of
#'   \code{metadata}. The first element is designated as day 1.
#'
#'   \code{fill_NA} is the function that fills missing data in \code{otudata},
#'   based on the boolean value \code{interpolate_NA} supplied to
#'   \code{prepanel()}. This will be set to either assign missing data values to
#'   zero (\code{interpolate_NA == TRUE}) or interpolate them using adjacent
#'   data within the same OTU (\code{interpolate_NA == FALSE}). (Note: missing
#'   data in this case does not include entire timepoints without data, in which
#'   case data is either interpolated or a break in the time axis is created.)
#'
#' @param aesthetics A list of custom aesthetics to apply to the horizon plot.
#'   This should come directly from the output list of \code{horizonaes()},
#'   without alteration.
#'
#' @return Returns the horizon plot as a \code{ggplot} object.
#'
#' @examples
#' # Basic plot form. By default, samples are plotted next to each other.
#' plist <- prepanel(otudata = otusample_diet, metadata = metadatasample_diet,
#' taxonomysample = taxonomysample_diet, subj = "MCTs16")
#' horizonplot(plist)
#'
#' # For irregularly spaced time series, you can "regularize" the data to create
#' # an accurate timescale.
#'
#' # Adjust data to regular time intervals each 1 day. This will interpolate
#' # new data points for each OTU at day = 1, 2, 3 etc. based on values
#' # at previous and subsequent timepoints.
#' plist <- prepanel(otudata = otusample_diet, metadata = metadatasample_diet, 
#'                   subj = "MCTs16", regularInterval = 1)
#' horizonplot(plist)
#'
#' # If the data has large gaps of time without samples, interpolating data
#' # within these time intervals could be misleading. You can set a maximum
#' # amount of time without samples allowed to plot a timepoint. If a timepoint
#' # is eliminated, a break in the time axis will be created at that point, and
#' # data will be regularized separately on both sides of the break in two
#' # different facets.
#'
#' # Set maximum time without samples to 75 days
#' plist <- prepanel(otudata = otusample_baboon, metadata = metadatasample_baboon, 
#'                   subj = "Baboon_388", regularInterval = 25, maxGap = 75)
#' horizonplot(plist)
#'
#' # Remove facets with less than 5 samples
#' plist <- prepanel(otudata = otusample_baboon, metadata = metadatasample_baboon, 
#'                   subj = "Baboon_388", regularInterval = 25, maxGap = 75, 
#'                   minSamplesPerFacet = 5)
#' horizonplot(plist)
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr %>%
#'
#' @export
horizonplot <- function(parameterList, aesthetics=horizonaes()) {
  if(length(parameterList) != 15) {
    stop("`parameterList` must be of length 15; use output list from prepanel function")
  }

  otudata <- parameterList[[1]]
  taxonomydata <- parameterList[[2]]
  timestamps <- parameterList[[3]]
  otulist <- parameterList[[4]]
  subj <- parameterList[[5]]
  regularInterval <- parameterList[[6]]
  maxGap <- parameterList[[7]]
  minSamplesPerFacet <- parameterList[[8]]
  band.thickness <- parameterList[[9]]
  origin <- parameterList[[10]]
  facetLabelsByTaxonomy <- parameterList[[11]]
  customFacetLabels <- parameterList[[12]]
  fill_NA <- parameterList[[13]]
  nbands <- parameterList[[14]]
  formatStep <- parameterList[[15]]

  # Origin-center, save band-thicknesses into new column.
  otudata$bt <- numeric(nrow(otudata))
  for(i in 1:nrow(otudata)) {
    otudata[i,2:(ncol(otudata)-1)] <- unlist(otudata[i,2:(ncol(otudata)-1)]) - origin(unlist(otudata[i,2:(ncol(otudata)-1)]))
    otudata$bt[i] <- band.thickness(unlist(otudata[i,2:(ncol(otudata)-1)]))
  }

  # Change sample IDs to corresponding dates and melt data
  if(!is.na(regularInterval)) {
    newTimestamps <- seq(from=1, to=max(timestamps), by=regularInterval)
    indPrev <- numeric()
    for(i in 1:length(newTimestamps)) {
      indPrev[i] <- max(which(timestamps <= newTimestamps[i]))
    }

    templist <- list()
    for(i in 1:length(newTimestamps)) {
      y1 <- otudata[,1+indPrev[i]]
      y2 <- otudata[,2+indPrev[i]]
      t1 <- timestamps[indPrev[i]]
      t2 <- timestamps[indPrev[i]+1]
      Ti <- newTimestamps[i]

      templist[[i]] <- pmin(y1,y2) + (Ti-t1)*(abs(y2-y1)/(t2-t1))

      if(!is.na(maxGap) & !(Ti==t1 | Ti==t2) & (abs(Ti-t1) > maxGap | abs(Ti-t2) > maxGap)) {
        templist[[i]] <- NA
      }
    }
    newdata <- as.data.frame(templist)
    rm(templist)
    colnames(newdata) <- newTimestamps
    otudata <- cbind(otuid=otudata$otuid,newdata,bt=otudata$bt)
  } else {
    colnames(otudata)[2:(ncol(otudata)-1)] <- as.character(1:(ncol(otudata)-2))
  }

  # Fill NAs using selected method
  if(is.na(regularInterval) | is.na(maxGap)) {
    otudata[,2:(ncol(otudata)-1)] <- as.data.frame(t(fill_NA(t(otudata[,2:(ncol(otudata)-1)]))))
  } else {
    nums <- otudata %>% dplyr::select(-otuid, -bt)
    breakpts <- as.numeric(c(-1*regularInterval,colnames(nums)[colSums(is.na(nums))==nrow(nums)],regularInterval+as.numeric(colnames(nums)[ncol(nums)])))
    naCols <- as.character(breakpts[c(-1,-length(breakpts))])
    otudata <- otudata %>% dplyr::select(-naCols)
    rm(naCols)
  }

  # formatStep = FALSE: add extra points for interpolation
  extrapts <- data.frame(character(),numeric(),numeric(),numeric())
  colnames(extrapts) <- c("otuid","bt","day","value")
  dX <- ifelse(is.na(regularInterval), 1, regularInterval)

  if(!formatStep) {
    for(r in 1:nrow(otudata)) {
      bt <- otudata[r,ncol(otudata)]
      nums <- otudata[r,c(-1,-ncol(otudata))]

      bands <- 1 # Every point is at least a member of the lowest band
      for(n in (-1*nbands+1):(nbands-1)) {
        bands <- bands + (nums > n*bt)
      }
      bands <- as.numeric(bands-nbands)

      x <- numeric()
      y <- numeric()

      for(i in 1:(2*nbands-1)) {
        # UPCROSSINGS
        upcrossings <- as.numeric(which(bands-dplyr::lead(bands) == -i)) # The indices where the sample is about to increase band
        if(length(upcrossings>0)) {
          y1 <- as.numeric(nums[upcrossings])
          y2 <- as.numeric(nums[upcrossings+1])
          for(j in 0:(i-1)) {
            # Add a point.
            newY <- bt*(bands[upcrossings]+j)
            newX <- as.numeric(names(nums)[upcrossings]) + ((newY-y1)/(y2-y1)) * dX

            x <- c(x,newX)
            y <- c(y,newY)
          }
        }

        # DOWNCROSSINGS
        downcrossings <- as.numeric(which(bands-dplyr::lead(bands) == i)) # The indices where the sample is about to decrease band
        if(length(downcrossings>0)) {
          y1 <- as.numeric(nums[downcrossings])
          y2 <- as.numeric(nums[downcrossings+1])
          for(j in 0:(i-1)) {
            # Add a point.
            newY <- bt*(bands[downcrossings+1]+j)
            newX <- as.numeric(names(nums)[downcrossings])+dX - ((newY-y2)/(y1-y2)) * dX

            x <- c(x,newX)
            y <- c(y,newY)
          }
        }
      }

      if(length(x)>0) {
        extrapts <- rbind.data.frame(extrapts, data.frame(otuid=otudata[r,1],bt=bt,day=as.character(x),value=y, stringsAsFactors=FALSE))
      }
    }
  }

  # Melt data into a date column
  otudata <- otudata %>% tidyr::gather(key=day,value=value,-c("otuid","bt"))
  otudata$otuid <- as.character(otudata$otuid)
  otudata$day <- as.character(otudata$day)

  otudata <- rbind.data.frame(otudata,extrapts) # Merge regular points and intermediate points

  # Add band columns
  for (i in 1:nbands) {
    #positive
    otudata[,paste("ypos",i,sep="")] <- ifelse(otudata$value > 0,
                                               ifelse(abs(otudata$value) > otudata$bt * i,
                                                      otudata$bt,
                                                      ifelse(abs(otudata$value) - (otudata$bt * (i - 1)) > 0, abs(otudata$value) - (otudata$bt * (i - 1)), 0)), 0)
    #negative
    otudata[,paste("yneg",i,sep="")] <- ifelse(otudata$value < 0,
                                               ifelse(abs(otudata$value) > otudata$bt * i,
                                                      otudata$bt,
                                                      ifelse(abs(otudata$value) - (otudata$bt * (i - 1)) > 0, abs(otudata$value) - (otudata$bt * (i - 1)), 0)), 0)
  }

  # Melt data and graph it
  otudata <- otudata[,c(1,3,5:ncol(otudata))] %>% tidyr::gather(key=band,value=value,-(1:2))

  if(formatStep) {
    # Add additional steps to dataframe
    otudata <- otudata %>% dplyr::arrange(otuid)
    otudata_extraSteps <- otudata %>% dplyr::mutate(value = ifelse(day==1,NA,dplyr::lag(value)))
    otudata <- dplyr::bind_rows(old = otudata, new = otudata_extraSteps, .id = "source") %>%
      dplyr::arrange(otuid, band, day, source)
  }

  # Add breaks in timescale
  if(!is.na(maxGap) && length(breakpts)!=0) {
    if(is.na(minSamplesPerFacet) || length(minSamplesPerFacet) > 1) {
      minSamplesPerFacet <- 2
    }
    if(minSamplesPerFacet < 2) {
      warning("minSamplesPerFacet must be >=2; using default value")
      minSamplesPerFacet <- 2
    }
    if(!(class(minSamplesPerFacet) %in% c("numeric","integer")) || minSamplesPerFacet %% 1 != 0) {
      warning("minSamplesPerFacet must be an integer; using default value")
      minSamplesPerFacet <- 2
    }
    otudata <- addFacets(otudata,breakpts,regularInterval, minSamplesPerFacet)
  }


  # Arrange OTU rows according to order of otulist
  if(length(otulist)>1 || !is.na(otulist)) {
    otudata$otuid <- factor(otudata$otuid,levels=otulist)
    otudata <- otudata %>% dplyr::arrange(as.numeric(day),band,otuid)
  }
  # Arrange OTU rows according to order of subjlist
  if(length(subj) > 1) {
    otudata$otuid <- factor(otudata$otuid,levels=subj)
    otudata <- otudata %>% dplyr::arrange(as.numeric(day),band,otuid)
  }

  # Add facet labels
  if(length(customFacetLabels) > 1 || !is.na(customFacetLabels)) {
    facetLabels <- customFacetLabels
  } else {
    if(facetLabelsByTaxonomy & is.data.frame(taxonomydata)) {
      facetLabels <- taxonomydata %>%
        dplyr::arrange(factor(otuid,levels=unique(otudata$otuid))) %>%
        apply(1,function(x){dplyr::last(x[!is.na(x)])})
    } else {
      facetLabels <- as.character(unique(otudata$otuid))
      if(any(nchar(facetLabels) > 20)) {
        warning("One or more OTU IDs have long names; default facet labels will be truncated to 20 characters.")
        facetLabels <- unname(vapply(facetLabels, {function(s) substring(s,1,20)}, character(1)))
      }
    }
  }
  names(facetLabels) <- unique(otudata$otuid)

  # Get aesthetics
  col.outline <- aesthetics[[1]]
  col.border <- aesthetics[[2]]
  col.bands <- aesthetics[[3]]
  legendPosition <- aesthetics[[4]]
  aesthetics <- aesthetics[-(1:4)]

  if(is.null(col.bands) || length(col.bands) != 2*nbands) {
    if(!is.null(col.bands) && length(col.bands) != 2*nbands) {
      warning("Invalid length of `col.bands`, length must equal 2 * nbands\nUsing a default palette")
    }
    col.bands <- switch(as.character(nbands),
                       "3"=c("#B2182B", "#D6604D", "#F4A582", "#92C5DE", "#4393C3", "#2166AC"),
                       "4"=c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#92C5DE", "#4393C3", "#2166AC", "#053061"),
                       "5"=c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"),
                       NA)
    if(NA %in% col.bands) {
      stop("No default colorset for nbands > 5; provide a color palette for `col.bands` of length 2 * nbands")
    }
  }
  colorcodes <- numeric(nbands*2)
  for(i in 1:nbands) {
    colorcodes[i] <- paste("yneg",nbands-i+1,sep="")
  }
  for(i in (nbands+1):(nbands*2)) {
    colorcodes[i] <- paste("ypos",i-nbands,sep="")
  }
  names(col.bands) <- colorcodes

  p <- ggplot(data=otudata) +
    geom_area(aes(x = as.numeric(day), y = value, fill=band), position="identity", color=col.outline) +
    scale_fill_manual(values=col.bands,breaks=names(col.bands)[c((2*nbands):(1+nbands),nbands:1)],labels=c(paste("+",nbands:1,sep=""),(-1):(-1*nbands))) +
    theme_bw() +
    theme(axis.text.x=element_text(size=16), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid=element_blank(), panel.border=element_rect(color=col.border), strip.text.y.left=element_text(angle=0), panel.spacing.y=unit(0, units="cm"), legend.position=legendPosition) +
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) + # remove margins between plot and panel
    xlab(ifelse(is.na(timestamps), "Sample", "Day")) +
    ylab(element_blank())

  if(!is.na(maxGap) && maxGap!=0 && length(breakpts)!=0) {
    p <- p + facet_grid(otuid ~ type, scales="free", space="free_x", labeller=labeller(otuid=facetLabels), switch="y") + theme(strip.text.x=element_blank(), strip.background.x=element_blank()) # do new subplot for each otu
  } else {
    p <- p + facet_grid(otuid ~ ., scales="free", labeller=labeller(otuid=facetLabels), switch="y") # do new subplot for each otu, time cluster
  }
  p <- p + aesthetics

  outputNums <- sapply(c("regularInterval", "maxGap", "nbands"), function(x){get(x)})
  outputNums <- outputNums[!is.na(outputNums)]
  outputLgc <- c("facetLabelsByTaxonomy"=facetLabelsByTaxonomy, "formatStep"=formatStep)

  cat("Constructed a horizon plot with the following settings:",
      paste(names(outputNums), ": ", outputNums, sep=""),
      paste(names(outputLgc), ": ", outputLgc, sep=""),
      #paste("origin:", ifelse(is.na(originconstant), body(origin)[2], originconstant)),
      #paste("band.thickness:", ifelse(is.na(btconstant), body(band.thickness)[2], btconstant)),
      paste("fill_NA:", body(fill_NA)[2]),
      sep="\n")

  biomehorizonpkg_refined_otu <<- otudata
  cat("OTU table `biomehorizonpkg_refined_otu` was outputted to the environment.\n")

  if(length(taxonomydata) > 1 || !is.na(taxonomydata)) {
    biomehorizonpkg_taxonomy <<- taxonomydata
    cat("Taxonomy data table `biomehorizonpkg_taxonomy` was outputted to the environment.\n")
  }

  p
}

#' Add Custom Aesthetics to the Horizon Plot
#'
#' Add additional aesthetics to the horizon plot object, returning them in the
#' form of a list of ggplot aesthetics to supply to \code{horizonplot()}.
#'
#' Setting any aesthetic to \code{NA} will use the default value. For most
#' values, this means the aesthetic will be blank or not appear on the horizon
#' plot. Other values like xlabel will use their default text. If an aesthetic
#' is already in the horizon plot by default, and you want to remove it, you can
#' do so by setting the respective argument to \code{NULL}.
#'
#' This function provides an easy way to add the most common aesthetics to the
#' horizon plot, but if you want to add other aesthetics not included in this
#' function, you can do so by appending them to the horizon plot object using
#' the \code{+} operator. e.g. to add a gray background in the plotting area:\cr
#' \code{horizonplot(prepanel(otudata = otusample_diet, 
#' metadata = metadatasample_diet, taxonomydata = taxonomysample_diet, 
#' subj = "MCTs01")) + theme(panel.background = element_rect(fill="gray90"))}
#'
#' @param title character. The text for the title.
#' @param subtitle character. The text for the subtitle, displayed below the title.
#' @param xlabel character. The text for the x-axis label. If \code{NA}, uses a
#'   default label: "Sample" if data is not regularized, or "Day" if data is
#'   adjusted to regular time intervals.
#' @param ylabel character. The text for the y-axis label.
#' @param showColorLegend logical. If \code{TRUE} (the default) the color scale
#'   legend is drawn. Colors are displayed with the highest band on top and the
#'   lowest band on the bottom. For horizontal legends, higher bands will appear
#'   to the right and lower bands to the left.
#' @param showLegendLabels logical. If \code{TRUE} (the default), labels indicating
#'   band number are shown on the color scale legend.
#' @param legendPosition character. Where should the color scale legend be
#'   displayed? Possible values are "right" (the default), "left", "top", and
#'   "bottom". For the latter two options, the legend will appear horizontal.
#' @param legendTitle character. The text for the title of the color scale
#'   legend.
#' @param showPlotLabels logical. If \code{TRUE} (the default), labels for OTU
#'   subplots are shown.
#' @param col.bands character vector of hexadecimal color codes giving the color
#'   scale for horizon bands. Colors should be specified from the most negative
#'   band to the most positive band. Must be the same length as \code{2*nbands}.
#'   If \code{NA}, uses a default Red-Blue color gradient of a length determined
#'   by \code{nbands} in \code{prepanel()}. In this palette, darker shades of
#'   red indicate progressively more negative bands, while darker shades of blue
#'   indicate increasingly positive bands.
#' @param col.outline character string specifying the hexadecimal color code for
#'   the outline on top of the graph of each band. If \code{NULL}, the outline
#'   will be removed. Defaults to light gray (#CCCCCC).
#' @param col.border character string specifying the hexadecimal color code for
#'   panel borders. If \code{NULL}, panel borders will be removed. Defaults to
#'   light gray (#CCCCCC).
#'
#' @return A list containing custom ggplot aesthetics to override default values
#'   on the horizon plot. This list can then be supplied to \code{horizonplot()}
#'   to apply the aesthetics.
#'
#' @examples
#' plist <- prepanel(otudata = otusample_diet, metadata = metadatasample_diet, 
#' taxonomydata = taxonomysample_diet, subj = "MCTs01")
#'
#' # By default, the function is called with no arguments to use default aesthetics
#' horizonplot(plist, horizonaes())
#'
#' # Same plot as above
#' horizonplot(plist)
#'
#' # Add a custom title, ylabel
#' horizonplot(plist, horizonaes(title = "Microbiome Horizon Plot", ylabel = "OTU ID"))
#'
#' # Remove the default x-label
#' horizonplot(plist, xlabel = NULL)
#'
#' # Use a different colorscale
#' library(RColorBrewer)
#' horizonplot(plist, horizonaes(col.bands = brewer.pal(8, "PiYG")))
#'
#' # To add aesthetics not included in this function, append them to the
#' # horizon plot object. e.g. for a gray plotting area background:
#' horizonplot(plist) + theme(panel.background = element_rect(fill = "gray90"))
#'
#' @import ggplot2
#'
#' @export
horizonaes <- function(title=NA, subtitle=NA, xlabel=NA, ylabel=NA, showColorLegend=TRUE, showLegendLabels=TRUE, legendPosition="right", legendTitle=NA, showPlotLabels=TRUE, col.bands=NA, col.outline="#CCCCCC", col.border="#CCCCCC") {
  if(is.null(title)) {
    title <- NA
  }
  if(!is.na(title)) {
    if(!(class(title)=="character")) {
      warning("title must be of type character; removing title")
      title <- NA
    } else {
      if(!is.null(subtitle) && !is.na(subtitle)) {
        title <- ggtitle(title, subtitle)
      } else {
        title <- ggtitle(title)
      }
    }
  }

  if(is.null(xlabel)) {
    xlabel <- xlab(NULL)
  } else {
    if(!is.na(xlabel)) {
      if(!(class(xlabel)=="character")) {
        warning("xlabel must be of type character; using default xlab")
        xlabel <- NA
      } else {
        xlabel <- xlab(xlabel)
      }
    }
  }

  if(is.null(ylabel)) {
    ylabel <- NA
  }
  if(!is.na(ylabel)) {
    if(!(class(ylabel)=="character")) {
      warning("ylabel must be of type character; removing ylab")
      ylab <- NA
    } else {
      ylabel <- ylab(ylabel)
    }
  }

  if(is.null(showPlotLabels) || isFALSE(showPlotLabels)) {
    showPlotLabels <- theme(strip.background.y=element_blank(),strip.text.y=element_blank())
  } else {
    showPlotLabels <- NA
  }

  if(is.null(legendPosition)) {
    legendPosition <- "none"
  }
  if(is.na(legendPosition) || !(legendPosition %in% c("none","right","left","top","bottom"))) {
    legendPosition <- "right"
  }

  if(!is.null(showColorLegend) && !isFALSE(showColorLegend)) {
    if(!is.null(legendTitle) && (is.na(legendTitle) || !(class(legendTitle)=="character"))) {
      legendTitle <- ifelse(legendPosition %in% c("bottom","top"), "Band Colorscale", "Band\nColorscale")
    }
    if(legendPosition %in% c("bottom","top")) {
      showColorLegend <- ggplot2::guides(fill=ggplot2::guide_legend(title=legendTitle,label=ifelse(isFALSE(showLegendLabels),FALSE,TRUE),label.position="bottom",nrow=1,reverse=TRUE))
    } else {
      showColorLegend <- ggplot2::guides(fill=ggplot2::guide_legend(title=legendTitle,label=ifelse(isFALSE(showLegendLabels),FALSE,TRUE)))
    }
  } else {
    showColorLegend <- ggplot2::guides(fill=FALSE)
  }

  if(is.null(col.bands) || NA %in% col.bands) {
    col.bands <- NULL
  } else {
    if(class(col.bands) != "character") {
      warning("col.bands must be a hexadecimal color vector of type character. Using a default colorset.")
      col.bands <- NULL
    }
    if(length(col.bands) < 6) {
      warning("col.bands has <6 color values, not enough for at least 3 bands. Using a default colorset.")
      col.bands <- NULL
    }
    if(length(col.bands) %% 2 != 0) {
      warning("col.bands has an odd number of color values. Using a default colorset.")
      col.bands <- NULL
    }
    if(!is.null(col.bands)) {
      if(nchar(col.bands)!=7 || substr(col.bands,1,1)!="#") {
        warning("col.bands must be in hexadecimal format. Using a default colorset.")
        col.bands <- NULL
      }
    }
  }

  if(!is.null(col.outline)) {
    if(is.na(col.outline)) {
      col.outline <- "#CCCCCC"
    } else {
      if(class(col.outline) != "character") {
        stop("col.outline must be a hexadecimal color of type character")
      } else {
        if(nchar(col.outline) != 7) {
          stop("col.outline must be in hexadecimal color format (e.g. #CCCCCC)")
        }
      }
    }
  }

  if(!is.null(col.border)) {
    if(is.na(col.border)) {
      col.border <- "#CCCCCC"
    } else {
      if(class(col.outline) != "character") {
        stop("col.border must be a hexadecimal color of type character")
      } else {
        if(nchar(col.outline) != 7) {
          stop("col.border must be in hexadecimal color format (e.g. #CCCCCC)")
        }
      }
    }
  }

  aesthetics <- list(col.outline,col.border,col.bands,legendPosition,title,xlabel,ylabel,showColorLegend, showPlotLabels)
  aesthetics <- aesthetics[!is.na(aesthetics)]
  aesthetics
}

# Internal function used to create breaks in the time scale through facetting.
addFacets <- function(x, boundaries, regularInterval, minSamplesPerFacet) {
  if(!is.data.frame(x)) {
    stop("x must be of type data frame")
  }

  for(i in 1:nrow(x)) {
    for(j in 1:length(boundaries)) {
      if(as.numeric(x$day[i]) >= boundaries[j]+regularInterval & as.numeric(x$day[i]) <= boundaries[j+1]-regularInterval) {
        x$type[i] <- j
        break()
      }
      x$type[i] <- 0
    }
  }

  # filter out facets without enough samples (set type=0)
  x <- x %>% group_by(type) %>% do(mutate(.,type=ifelse(type==0 || sum(as.numeric(unique(day)) %% 1 == 0)<minSamplesPerFacet,0,type))) %>% as.data.frame()

  x %>% dplyr::filter(type!=0)
}

