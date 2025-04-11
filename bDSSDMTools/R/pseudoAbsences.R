#' Generate Pseudo-Absence Points Excluding Buffer Around Presences
#'
#' Creates a SpatVector of pseudo-absence points by randomly sampling locations
#' from the non-NA cells of a template raster layer, while excluding areas
#' within a specified buffer distance around known presence points.
#'
#' @export
#' @importFrom terra SpatRaster SpatVector crs buffer rasterize values `values<-` xyFromCell vect is.projected ncell `[` subset `names<-` units `units<-` wrap unwrap deepcopy `ext<-` ext maths Compare Logic Arith `Arith,SpatRaster,SpatRaster-method` `Arith,SpatRaster,numeric-method` `Arith,numeric,SpatRaster-method` `Compare,SpatRaster,numeric-method` `Compare,SpatRaster,SpatRaster-method` `Compare,numeric,SpatRaster-method` `Logic,SpatRaster,logical-method` `Logic,SpatRaster,SpatRaster-method` `Logic,!SpatRaster-method` subset extract `subset,SpatRaster,numeric-method` `subset,SpatRaster,missing-method` `extract,SpatRaster,SpatVector-method` `extract,SpatRaster,matrix-method` `is.na,SpatRaster-method` `anyNA,SpatRaster-method`
#' @importFrom methods is

pseudoAbsences <- function(rasterLayers,presences,n=1000) {

  shape <- rasterLayers[[1]]
  crs(presences) <- crs(shape)

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

