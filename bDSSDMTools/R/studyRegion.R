#' Define Study Region Mask Based on Multiple Criteria
#'
#' Creates a spatial mask (SpatRaster with values 1 or NA) defining a study
#' region. The region is determined by the intersection (AND logic) of several
#' optional criteria: inclusion within a polygon mask, inclusion based on depth
#' preferences, inclusion based on an intertidal layer, and/or inclusion within
#' a specified distance threshold around presence points.
#'
#' @param rasterLayers A terra `SpatRaster` object. The first layer is used as a
#'   template for the output geometry (extent, resolution, CRS) and initial
#'   non-NA area. Must have a defined CRS.
#' @param presences A terra `SpatVector` object with point geometry. Required only
#'   if `distanceThreshold` is specified. Must have the same CRS as `rasterLayers`.
#' @param distanceThreshold Numeric or NULL (Default: 200000). If numeric, defines
#'   the maximum distance (in CRS units, typically meters) from presence points
#'   to include in the study region. Cells outside this buffer will be masked out.
#'   Setting to NULL disables this criterion. Requires projected CRS.
#' @param bathymetryLayer Optional. A `terra SpatRaster` or a character string
#'   path to a bathymetry raster file. Used with `depthPref`. Must have the
#'   same CRS as `rasterLayers`.
#' @param depthPref Optional. A numeric vector of length 2: `c(min_depth, max_depth)`.
#'   Defines the acceptable depth range. Values should typically be non-positive
#'   (<= 0). Cells outside this depth range will be masked out. Requires
#'   `bathymetryLayer` to be provided.
#' @param intertidalLayer Optional. A `terra SpatRaster` or a character string
#'   path to an intertidal zone raster file. Cells that are NA in this layer
#'   (after cropping) will be masked out from the study region. Must have the
#'   same CRS as `rasterLayers`.
#' @param maskPolygon Optional. A `terra SpatVector` or `sf` polygon object.
#'   Defines the outer boundary. Cells outside this polygon will be masked out.
#'   Must have the same CRS as `rasterLayers`.
#' @param trimToRegion Logical (Default: TRUE). If TRUE, the final output mask
#'   will be trimmed to the smallest extent containing non-NA cells using `terra::trim()`.
#'
#' @return A `terra SpatRaster` mask with value 1 for cells within the finally
#'   defined study region and NA elsewhere. It retains the CRS of the input `rasterLayers`.
#'   Returns NULL or stops if essential inputs are missing or invalid.
#' @export
#' @importFrom terra SpatRaster SpatVector crs buffer rasterize values `values<-` xyFromCell vect is.projected ncell `[` subset `names<-` units `units<-` wrap unwrap deepcopy `ext<-` ext maths Compare Logic Arith `Arith,SpatRaster,SpatRaster-method` `Arith,SpatRaster,numeric-method` `Arith,numeric,SpatRaster-method` `Compare,SpatRaster,numeric-method` `Compare,SpatRaster,SpatRaster-method` `Compare,numeric,SpatRaster-method` `Logic,SpatRaster,logical-method` `Logic,SpatRaster,SpatRaster-method` `Logic,!SpatRaster-method` subset extract mask classify crop trim rast `rast,character-method` project `project,SpatVector,character-method` `ifel` `NAflag<-` `NAflag` `crs` `is.na,SpatRaster-method` `anyNA,SpatRaster-method`
#' @importFrom methods is
#'

studyRegion <- function(rasterLayers, presences, distanceThreshold=200000, bathymetryLayer=NULL, depthPref=NULL, intertidalLayer=NULL, maskPolygon=NULL, trimToRegion=TRUE) {

  forceStudyRegion <- numeric(0)
  nonStudyRegion <- numeric(0)

  shape <- rasterLayers[[1]]
  shape[] <- NA
  names(shape) <- "studyRegion"

  if( ! is.null(maskPolygon)) {

    if( class(maskPolygon) == "sf" ) { maskPolygon <- vect(maskPolygon) }
    studyRegionLayer.i <- terra::mask(shape,maskPolygon)
    forceStudyRegion <- c(forceStudyRegion,which(!is.na(studyRegionLayer.i[][,1])))

  }

  if( ! is.null(intertidalLayer) ) {

    if( class(intertidalLayer) == "character" ) { intertidalLayer <- rast(intertidalLayer) }

    intertidalLayer <- terra::crop(intertidalLayer, shape)
    forceStudyRegion <- c(forceStudyRegion,which(!is.na(intertidalLayer[][,1])))

  }

  if( ! is.null(bathymetryLayer) & ! is.null(depthPref) ) {

    if( max(depthPref) > 0 ) { stop("Error :: depthPref must be lower or equal to 0") }
    if( length(depthPref) != 2 ) { stop("Error :: depthPref must be a vector of length 2") }

    if( class(bathymetryLayer) == "character" ) { bathymetryLayer <- rast(bathymetryLayer) }

    bathymetryLayer <- terra::crop(bathymetryLayer, shape)

    m_specific <- data.frame(from=min(depthPref), to=max(depthPref), becomes=1)
    bathymetryLayer.i <- terra::classify(bathymetryLayer, rcl = m_specific, others = NA)
    forceStudyRegion <- c(forceStudyRegion,which(!is.na(bathymetryLayer.i[][,1])))

  }

  if( ! is.null(distanceThreshold) ) {

    crs(presences, warn=FALSE) <- shape

    # Make a buffer on presences using distanceThreshold in meters
    presencesBuffer <- terra::buffer(presences, width=distanceThreshold)
    presencesBufferRaster <- terra::rasterize(presencesBuffer, shape, field=1)
    nonStudyRegion <- c(nonStudyRegion,which(is.na(presencesBufferRaster[][,1])))

  }

  shape[unique(forceStudyRegion)] <- 1
  shape[unique(nonStudyRegion)] <- NA

  if( trimToRegion ) { shape <- trim(shape) }

  return(shape)
}

