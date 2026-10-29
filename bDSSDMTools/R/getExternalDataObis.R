#' Download Occurrence Data from OBIS
#'
#' Fetches species occurrence records from the Ocean Biodiversity Information
#' System (OBIS) using the `robis` package. Optionally retrieves dataset
#' citation information via the OBIS API.
#'
#' @param taxa Scientific name, preferably a WoRMS name.
#' @param getCitation If `TRUE`, add access rights and a bibliographic citation
#'   for each dataset. Records are still returned when a citation request fails.
#'
#' @return A data frame with scientific name, coordinates, date, access rights,
#'   and citation. Returns `NULL` when OBIS has no matching records.
#' @export
#' @importFrom stringi stri_trans_general
#' @importFrom robis occurrence
#'
getExternalDataObis <- function(taxa, getCitation = TRUE) {

  if (missing(taxa) || !is.character(taxa) || length(taxa) != 1 || !nzchar(trimws(taxa))) {
    stop("Error :: no taxa (WoRMS name) introduced")
  }
  taxa <- stringi::stri_trans_general(str = taxa, id = "Latin-ASCII")

  my_occs_obis <- tryCatch(
    robis::occurrence(scientificname = taxa),
    error = function(e) NULL
  )
  if (is.null(my_occs_obis)) my_occs_obis <- data.frame()
  my_occs_obis <- as.data.frame(my_occs_obis)

  if (nrow(my_occs_obis) > 0 && all(c("decimalLongitude", "decimalLatitude") %in% names(my_occs_obis))) {
    keep <- my_occs_obis$decimalLongitude != 0 & my_occs_obis$decimalLatitude != 0
    keep[is.na(keep)] <- FALSE
    my_occs_obis <- my_occs_obis[keep, , drop = FALSE]
  }

  if (nrow(my_occs_obis) == 0) {
    cat("Data downloaded from OBIS\n")
    cat("Taxa: ", taxa, "\n")
    cat("Number of records: 0\n")
    return(NULL)
  }

  my_occs_obis$accessRights <- NA_character_
  my_occs_obis$bibliographicCitation <- NA_character_

  if (isTRUE(getCitation) && "dataset_id" %in% names(my_occs_obis)) {
    for (datasetId in unique(my_occs_obis$dataset_id)) {
      z.Res <- NULL
      for (attempt in seq_len(10)) {
        z.Res <- tryCatch(
          RJSONIO::fromJSON(paste0("https://api.obis.org/v3/dataset/", datasetId)),
          error = function(e) NULL
        )
        if (!is.null(z.Res)) break
      }
      if (is.null(z.Res) || is.null(z.Res$results) || length(z.Res$results) < 1) next

      rows <- my_occs_obis$dataset_id == datasetId
      record <- z.Res$results[[1]]
      oneValue <- function(column) {
        if (!column %in% names(my_occs_obis)) return("")
        values <- unique(as.character(my_occs_obis[[column]][rows]))
        values <- values[!is.na(values) & nzchar(values)]
        if (length(values) == 0) "" else values[1]
      }
      my_occs_obis$accessRights[rows] <- if (!is.null(record$intellectualrights)) record$intellectualrights else NA_character_
      citation <- paste0(
        if (!is.null(record$citation)) record$citation else "",
        " ",
        oneValue("institutionCode"),
        " ",
        oneValue("collectionCode"),
        " (Available: Ocean Biogeographic Information System. Intergovernmental Oceanographic Commission of UNESCO. www.iobis.org. Accessed: ",
        Sys.time(),
        ")"
      )
      my_occs_obis$bibliographicCitation[rows] <- citation
    }
  }

  df <- .occurrenceTable(my_occs_obis, "decimalLongitude", "decimalLatitude")
  cat("Data downloaded from OBIS\n")
  cat("Taxa: ", taxa, "\n")
  cat("Number of records: ", nrow(df), "\n")
  df

}
