#' Generate Random Background Points from Raster Extent
#'
#' Samples a specified number of random background points from the non-NA cells
#' within the spatial extent defined by the input raster layers. The first layer
#' of the input SpatRaster is used as the template for identifying valid sampling locations.
#'
#' @export
#' @importFrom terra SpatRaster SpatVector crs values xyFromCell vect ncell `[` subset `names<-` units `units<-` wrap unwrap deepcopy `ext<-` ext maths Compare Logic Arith `Arith,SpatRaster,SpatRaster-method` `Arith,SpatRaster,numeric-method` `Arith,numeric,SpatRaster-method` `Compare,SpatRaster,numeric-method` `Compare,SpatRaster,SpatRaster-method` `Compare,numeric,SpatRaster-method` `Logic,SpatRaster,logical-method` `Logic,SpatRaster,SpatRaster-method` `Logic,!SpatRaster-method` subset extract `subset,SpatRaster,numeric-method` `subset,SpatRaster,missing-method` `extract,SpatRaster,SpatVector-method` `extract,SpatRaster,matrix-method` `is.na,SpatRaster-method` `anyNA,SpatRaster-method`
#' @importFrom methods is
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

