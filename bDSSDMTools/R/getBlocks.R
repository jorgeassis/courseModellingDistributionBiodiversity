#' Assign Spatial Block Membership using K-Means Clustering
#'
#' Clusters presence coordinates with k-means and assigns every absence to the
#' nearest cluster. Longitude-latitude coordinates are converted to approximate
#' kilometres before distances are calculated. The returned coordinates remain
#' longitude and latitude.
#'
#' The clustering uses a fixed seed internally and restores the session random
#' seed afterwards.
#'
#' @param records A terra `SpatVector` of points with a numeric `PA` column.
#'   Presences are `1` and absences are `0`.
#' @param kFolds Number of clusters. Cannot exceed the number of presences.
#'
#' @return A list with one element, `blocks`: a data frame with columns `PA`,
#'   `Lon`, `Lat`, and `membership`.
#' @export
#' @importFrom terra crds is.lonlat
#' @importFrom stats kmeans
#'
getBlocks <- function(records, kFolds) {

  if (!inherits(records, "SpatVector")) {
    stop("Error :: records must be a SpatVector object")
  }
  if (!"PA" %in% names(records)) {
    stop("Error :: records must contain a column named PA")
  }
  if (!is.numeric(kFolds) || length(kFolds) != 1 || kFolds < 1) {
    stop("Error :: kFolds must be a positive number")
  }
  kFolds <- as.integer(kFolds)

  presences <- records[records$PA == 1, ]
  absences <- records[records$PA == 0, ]
  if (nrow(presences) < kFolds) {
    stop("Error :: kFolds (", kFolds, ") cannot exceed the number of presences (", nrow(presences), ")")
  }

  presenceCoords <- terra::crds(presences)
  absenceCoords <- if (nrow(absences) > 0) terra::crds(absences) else matrix(numeric(), ncol = 2)
  lonlat <- isTRUE(terra::is.lonlat(records))
  lat0 <- mean(c(presenceCoords[, 2], absenceCoords[, 2]), na.rm = TRUE) * pi / 180
  toDistance <- function(xy) {
    if (nrow(xy) == 0) return(xy)
    if (!lonlat) return(xy)
    cbind(xy[, 1] * 111.32 * cos(lat0), xy[, 2] * 110.574)
  }

  hasSeed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (hasSeed) oldSeed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (hasSeed) {
      assign(".Random.seed", oldSeed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(123)

  kmeansResult <- stats::kmeans(toDistance(presenceCoords), centers = kFolds, nstart = 25)
  centres <- kmeansResult$centers

  presenceMembership <- kmeansResult$cluster
  absenceMembership <- integer(0)
  if (nrow(absences) > 0) {
    absenceKm <- toDistance(absenceCoords)
    absenceMembership <- apply(absenceKm, 1, function(point) {
      which.min(rowSums(sweep(centres, 2, point, "-")^2))
    })
  }

  combCoords <- rbind(
    data.frame(PA = rep(1, nrow(presences)), Lon = presenceCoords[, 1], Lat = presenceCoords[, 2], membership = presenceMembership),
    data.frame(PA = rep(0, nrow(absences)), Lon = absenceCoords[, 1], Lat = absenceCoords[, 2], membership = absenceMembership)
  )
  rownames(combCoords) <- NULL

  list(blocks = combCoords)

}
