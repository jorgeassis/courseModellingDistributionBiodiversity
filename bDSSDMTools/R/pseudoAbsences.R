#' Generate Pseudo-Absence Points Excluding Buffer Around Presences
#'
#' Creates a SpatVector of pseudo-absence points by randomly sampling locations
#' from the non-NA cells of a template raster layer, while excluding areas
#' within a specified buffer distance around known presence points.
#'
#' @param rasterLayers A terra `SpatRaster`. The first layer supplies the cells
#'   that can be sampled.
#' @param presences A terra `SpatVector` of presence points. Cells within 50 km
#'   of these points are not sampled.
#' @param n Number of pseudo-absence points to return.
#' @export
#' @importFrom terra values xyFromCell vect

pseudoAbsences <- function(rasterLayers,presences,n=1000) {

  shape <- rasterLayers[[1]]
  terra::crs(presences) <- terra::crs(shape)

  # Make a buffer on presences using distanceThreshold in meters
  presencesBuffer <- terra::buffer(presences, width=50000)
  presencesBufferRaster <- terra::rasterize(presencesBuffer, shape, field=1)
  shape[which(!is.na(presencesBufferRaster[][,1]))] <- NA

  nonNACells <- which(!is.na( values(shape)[,1] ))
  absences <- xyFromCell(shape, nonNACells)

  absences.n <- sample( 1:nrow(absences) , min(n,nrow(absences)) , replace=FALSE)
  absences <- absences[absences.n,]
  colnames(absences) <- c("Lon","Lat")

  absences <- vect(as.data.frame(absences), geom=c("Lon","Lat"))

  return(absences)

}

