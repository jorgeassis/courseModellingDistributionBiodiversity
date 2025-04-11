#' Download Occurrence Data from OBIS
#'
#' Fetches species occurrence records from the Ocean Biodiversity Information
#' System (OBIS) using the `robis` package. Optionally retrieves dataset
#' citation information via the OBIS API.
#'
#' @export
#' @importFrom stringi stri_trans_general
#' @importFrom robis occurrence
#' @importFrom jsonlite fromJSON
#'

getExternalDataObis <- function(taxa,getCitation) {

  df <- NULL
  taxa <- stri_trans_general(str = taxa, id = "Latin-ASCII")

  if( missing(taxa)) { errormessage("no taxa (Worms name) introduced.") }

  if( exists("my_occs_obis") ) { rm(my_occs_obis) }

  tryCatch( my_occs_obis <- robis::occurrence(scientificname = taxa  ) , error=function(e) { error <- TRUE })

  if( exists("my_occs_obis") ) { if( nrow(my_occs_obis) == 0 ) { my_occs_obis <- data.frame() } }

  if( ! exists("my_occs_obis") ) { my_occs_obis <- data.frame() }

  if( nrow(my_occs_obis) > 0) {

    my_occs_obis <- subset(my_occs_obis, my_occs_obis$decimalLongitude !=0 & my_occs_obis$decimalLatitude !=0)

  }

  my_occs_obis <- data.frame(my_occs_obis)

  if( getCitation ) {

    if( nrow(my_occs_obis) > 0) {

      my_occs_obisInfo <- my_occs_obis$dataset_id
      my_occs_obisInfo <- unique(my_occs_obis$dataset_id)

      for(z in 1:length(my_occs_obisInfo) ) {

        error <- TRUE
        errortrials <- 0

        while(error & errortrials < 10) {
          error <- FALSE
          errortrials <- errortrials + 1
          tryCatch(  z.Res <- RJSONIO::fromJSON(paste0("https://api.obis.org/v3/dataset/",my_occs_obisInfo[z])) , error=function(e) { error <- TRUE })
        }

        if(!error) {

          institutionCode <- my_occs_obis[my_occs_obis$dataset_id == my_occs_obisInfo[z],"institutionCode"]
          collectionCode <- my_occs_obis[my_occs_obis$dataset_id == my_occs_obisInfo[z],"collectionCode"]

          my_occs_obis[my_occs_obis$dataset_id == my_occs_obisInfo[z],"accessRights"] <- z.Res$results[[1]]$intellectualrights

          z.Res <- paste0( z.Res$results[[1]]$citation,
                           ifelse(!is.null(z.Res$results[[1]]$citation)," ",""),
                           ifelse(!is.na(institutionCode) | !is.null(institutionCode) , institutionCode , ""),
                           " ",
                           ifelse(!is.na(collectionCode) | !is.null(collectionCode) , collectionCode , ""),
                           " (Available: Ocean Biogeographic Information System. Intergovernmental Oceanographic Commission of UNESCO. www.iobis.org. Accessed: ", Sys.time())

          my_occs_obis[my_occs_obis$dataset_id == my_occs_obisInfo[z],"bibliographicCitation"] <- z.Res

        }

      }

      df <- my_occs_obis

    }

  }

  df <- data.frame(my_occs_obis)
  df <- df[,c("scientificName","decimalLongitude","decimalLatitude","year","month","day","accessRights","bibliographicCitation")]
  names(df) <- c("scientificName","Lon","Lat","year","month","day","accessRights","bibliographicCitation")
  if(nrow(df) == 0 ) { df <- NULL}

  cat("Data downloaded from OBIS\n")
  cat("Taxa: ", taxa, "\n")
  cat("Number of records: ", nrow(df), "\n")

  return(df)

}

