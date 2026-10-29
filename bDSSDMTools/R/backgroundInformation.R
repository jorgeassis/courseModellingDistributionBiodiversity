#' Generate Random Background Points from Raster Extent
#'
#' Samples a specified number of random background points from the non-NA cells
#' within the spatial extent defined by the input raster layers. The first layer
#' of the input SpatRaster is used as the template for identifying valid sampling locations.
#'
#' @param rasterLayers A terra `SpatRaster`. The first layer supplies the cells
#'   that can be sampled.
#' @param n Number of background points to return. Fewer are returned when the
#'   raster has fewer non-NA cells.
#' @export
#' @importFrom terra values xyFromCell vect
#'

backgroundInformation <- function(rasterLayers,n=10000) {

  shape <- rasterLayers[[1]]

  nonNACells <- which(!is.na( values(shape)[,1] ))
  absences <- xyFromCell(shape, nonNACells)

  absences.n <- sample( 1:nrow(absences) , min(n,nrow(absences)) , replace=FALSE)
  absences <- absences[absences.n,]
  colnames(absences) <- c("Lon","Lat")

  absences <- vect(as.data.frame(absences), geom=c("Lon","Lat"))

  return(absences)

}

