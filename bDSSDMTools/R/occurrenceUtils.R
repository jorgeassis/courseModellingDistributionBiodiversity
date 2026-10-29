# Shared column layout for GBIF and OBIS downloads. Not exported.

.occurrenceTable <- function(df, longitudeNames, latitudeNames) {

  take <- function(candidates) {
    for (name in candidates) {
      if (name %in% names(df)) return(df[[name]])
    }
    rep(NA, nrow(df))
  }

  data.frame(
    scientificName = take(c("scientificName", "species", "scientificname")),
    Lon = take(longitudeNames),
    Lat = take(latitudeNames),
    year = take("year"),
    month = take("month"),
    day = take("day"),
    accessRights = take("accessRights"),
    bibliographicCitation = take("bibliographicCitation"),
    stringsAsFactors = FALSE
  )

}
