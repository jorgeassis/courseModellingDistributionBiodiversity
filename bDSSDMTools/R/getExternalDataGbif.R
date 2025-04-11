#' Download Occurrence Data from GBIF
#'
#' Fetches species occurrence records from the Global Biodiversity Information
#' Facility (GBIF). It handles potentially large datasets by downloading
#' records in chunks if the total count exceeds a threshold. Optionally retrieves
#' dataset citation information.
#'
#' @export
#' @importFrom stringi stri_trans_general
#' @importFrom dismo gbif
#' @importFrom utils download.file tempfile
#' @importFrom jsonlite fromJSON
#' @importFrom plyr rbind.fill
#' @importFrom rgbif gbif_citation
#'

getExternalDataGbif <- function(taxa, getCitation) {

  my_occs_gbif <- NULL
  if( exists("my_occs_gbif") ) { rm(my_occs_gbif) }

  taxa <- stri_trans_general(str = taxa, id = "Latin-ASCII")
  nRecords <- gbif(strsplit(as.character(taxa), " ")[[1]][1], strsplit(as.character(taxa), " ")[[1]][2], geo=T, removeZeros=T , download=FALSE, ntries=999)

  if( nRecords < 200000 ) { my_occs_gbif <- gbif(strsplit(as.character(taxa), " ")[[1]][1], strsplit(as.character(taxa), " ")[[1]][2], geo=T, removeZeros=T , download=TRUE, ntries=999, nrecs=100) }

  if( nRecords >= 200000 ) {

    seqListing <- seq(0,nRecords,by = 100)
    if(max(seqListing) < nRecords) { seqListing <- c(seqListing,nRecords) }
    parallelChunks <- data.frame(from = seqListing[-length(seqListing)], to = c(seqListing[-c(1,length(seqListing))] -1 , nRecords ) )

    tmpfile <- paste(tempfile(), ".json", sep = "")
    my_occs_gbif <- data.frame()

    for( ch in 1:nrow(parallelChunks)) {

      error <- TRUE
      chunck <- NULL

      while(error) {
        tryCatch( chunck <- download.file(paste0("https://api.gbif.org/v1/occurrence/search?scientificname=",strsplit(as.character(taxa), " ")[[1]][1],"+",strsplit(as.character(taxa), " ")[[1]][2],"&offset=",parallelChunks[ch,1],"&limit=100"), tmpfile, quiet = TRUE) , error=function(e) { error <- TRUE })
        if(!is.null("chunck")) { error <- FALSE }
      }

      json <- scan(tmpfile, what = "character", quiet = TRUE, sep = "\n", encoding = "UTF-8")
      json <- chartr("\a\v", "  ", json)
      x <- jsonlite::fromJSON(json)
      r <- x$results
      r <- r[, !sapply(r, class) %in% c("data.frame", "list")]
      rownames(r) <- NULL
      my_occs_gbif <- rbind.fill(my_occs_gbif,r)

    }
  }

  if( ! is.null(my_occs_gbif) & getCitation ) {

    my_occs_gbif_all <- unique(my_occs_gbif$datasetKey)

    for(z in 1:length(my_occs_gbif_all) ) {

      z.Res <- gbif_citation(x=my_occs_gbif_all[z])

      my_occs_gbif[my_occs_gbif$datasetKey == my_occs_gbif_all[z] ,"accessRights"] <- ifelse(!is.null(z.Res$rights),z.Res$rights,"")

      z.Res <- z.Res$citation$citation

      my_occs_gbif[my_occs_gbif$datasetKey == my_occs_gbif_all[z],"bibliographicCitation"] <- z.Res

    }

  }

  if( !is.null(my_occs_gbif)) { if(nrow(my_occs_gbif) == 0 ) { my_occs_gbif <- NULL}}

  my_occs_gbif <- data.frame(my_occs_gbif)
  my_occs_gbif <- my_occs_gbif[,c("scientificName","lon","lat","year","month","day","accessRights","bibliographicCitation")]
  names(my_occs_gbif) <- c("scientificName","Lon","Lat","year","month","day","accessRights","bibliographicCitation")
  if(nrow(my_occs_gbif) == 0 ) { my_occs_gbif <- NULL}

  cat("Data downloaded from GBIF\n")
  cat("Taxa: ", taxa, "\n")
  cat("Number of records: ", nrow(my_occs_gbif), "\n")

  return(my_occs_gbif)

}

