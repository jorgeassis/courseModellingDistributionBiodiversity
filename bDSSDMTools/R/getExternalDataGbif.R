#' Download Occurrence Data from GBIF
#'
#' Fetches species occurrence records from the Global Biodiversity Information
#' Facility (GBIF). At most 100,000 georeferenced records are downloaded.
#' Smaller sets are downloaded with `dismo`. Larger sets are read from the GBIF
#' occurrence search API, which cannot return more than 100,000 records, and
#' the first 100,000 are kept. Optionally retrieves dataset citation
#' information.
#'
#' @param taxa Scientific name with at least a genus and a specific epithet.
#' @param getCitation If `TRUE`, add access rights and a bibliographic citation
#'   for each dataset.
#'
#' @return A data frame with scientific name, coordinates, date, access rights,
#'   and citation. Returns `NULL` when GBIF has no matching records.
#' @export
#' @importFrom stringi stri_trans_general
#' @importFrom dismo gbif
#' @importFrom utils download.file URLencode
#' @importFrom jsonlite fromJSON
#' @importFrom plyr rbind.fill
#' @importFrom rgbif gbif_citation
#'
getExternalDataGbif <- function(taxa, getCitation = TRUE) {

  if (missing(taxa) || !is.character(taxa) || length(taxa) != 1 || !nzchar(trimws(taxa))) {
    stop("Error :: taxa must be a scientific name with a genus and a species")
  }
  taxa <- stringi::stri_trans_general(str = taxa, id = "Latin-ASCII")
  nameParts <- strsplit(trimws(taxa), "\\s+")[[1]]
  if (length(nameParts) < 2) {
    stop("Error :: taxa must include at least a genus and a species")
  }
  genus <- nameParts[1]
  species <- nameParts[2]

  nRecords <- dismo::gbif(genus, species, geo = TRUE, removeZeros = TRUE, download = FALSE, ntries = 9)
  if (length(nRecords) != 1 || is.na(nRecords)) {
    stop("Error :: GBIF did not return a record count for ", taxa)
  }

  # The occurrence search API rejects any request whose offset + limit exceeds 100000.
  maxDownload <- 100000

  if (nRecords <= maxDownload) {
    my_occs_gbif <- dismo::gbif(genus, species, geo = TRUE, removeZeros = TRUE, download = TRUE, ntries = 999, nrecs = 100)
  } else {
    cat("GBIF reports ", nRecords, " georeferenced records for ", taxa, ". Downloading the first ", maxDownload, ".\n", sep = "")
    downloadTimeout <- options(timeout = 600)
    on.exit(options(downloadTimeout), add = TRUE)
    pageSize <- 300
    offsets <- seq(0, maxDownload - 1L, by = pageSize)
    tmpfile <- paste0(tempfile(), ".json")
    on.exit(unlink(tmpfile), add = TRUE)
    pages <- list()
    scientificName <- utils::URLencode(paste(genus, species), reserved = TRUE)

    for (offset in offsets) {
      limit <- min(pageSize, maxDownload - offset)
      url <- paste0(
        "https://api.gbif.org/v1/occurrence/search?scientificname=", scientificName,
        "&hasCoordinate=true&offset=", offset, "&limit=", limit
      )
      downloaded <- FALSE
      for (attempt in seq_len(3)) {
        status <- tryCatch(
          utils::download.file(url, tmpfile, quiet = TRUE, mode = "wb"),
          error = function(e) NA_integer_
        )
        if (is.finite(status) && status == 0 && file.exists(tmpfile) && file.info(tmpfile)$size > 0) {
          downloaded <- TRUE
          break
        }
        Sys.sleep(5)
      }
      if (!downloaded) {
        cat("GBIF download stopped at offset ", offset, ". Returning the records already retrieved.\n", sep = "")
        break
      }

      json <- scan(tmpfile, what = "character", quiet = TRUE, sep = "\n", encoding = "UTF-8")
      json <- chartr("\a\v", "  ", json)
      x <- tryCatch(jsonlite::fromJSON(json), error = function(e) NULL)
      if (is.null(x)) {
        cat("GBIF download stopped at offset ", offset, ". Returning the records already retrieved.\n", sep = "")
        break
      }
      r <- x$results
      if (is.null(r) || nrow(r) == 0) break
      r <- r[, !vapply(r, function(column) inherits(column, c("data.frame", "list")), logical(1)), drop = FALSE]
      rownames(r) <- NULL
      pages[[length(pages) + 1]] <- r
      if (isTRUE(x$endOfRecords)) break
    }

    my_occs_gbif <- if (length(pages) == 0) data.frame() else do.call(plyr::rbind.fill, pages)

    if (nrow(my_occs_gbif) > 0) {
      lon <- if ("decimalLongitude" %in% names(my_occs_gbif)) my_occs_gbif$decimalLongitude else my_occs_gbif$lon
      lat <- if ("decimalLatitude" %in% names(my_occs_gbif)) my_occs_gbif$decimalLatitude else my_occs_gbif$lat
      my_occs_gbif <- my_occs_gbif[!(lon == 0 | lat == 0), , drop = FALSE]
    }
  }

  if (is.null(my_occs_gbif) || nrow(my_occs_gbif) == 0) {
    cat("Data downloaded from GBIF\n")
    cat("Taxa: ", taxa, "\n")
    cat("Number of records: 0\n")
    return(NULL)
  }

  if (!"accessRights" %in% names(my_occs_gbif)) my_occs_gbif$accessRights <- NA_character_
  if (!"bibliographicCitation" %in% names(my_occs_gbif)) my_occs_gbif$bibliographicCitation <- NA_character_

  if (isTRUE(getCitation) && "datasetKey" %in% names(my_occs_gbif)) {
    for (datasetKey in unique(my_occs_gbif$datasetKey)) {
      z.Res <- tryCatch(rgbif::gbif_citation(x = datasetKey), error = function(e) NULL)
      if (is.null(z.Res)) next
      rows <- my_occs_gbif$datasetKey == datasetKey
      my_occs_gbif$accessRights[rows] <- if (!is.null(z.Res$rights)) z.Res$rights else ""
      citation <- z.Res$citation$citation
      if (!is.null(citation)) my_occs_gbif$bibliographicCitation[rows] <- citation
    }
  }

  my_occs_gbif <- .occurrenceTable(my_occs_gbif, c("lon", "decimalLongitude"), c("lat", "decimalLatitude"))
  cat("Data downloaded from GBIF\n")
  cat("Taxa: ", taxa, "\n")
  cat("Number of records: ", nrow(my_occs_gbif), "\n")
  my_occs_gbif

}
