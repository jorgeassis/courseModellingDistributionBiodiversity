#' Prepare Environmental Data for Niche Overlap
#'
#' Builds the data table expected by [meanNicheOverlap()], [plotNicheOverlap()]
#' and [nicheSize()]. Presence cells (value `1`) are read from two binary
#' rasters, and the environmental layers are extracted at those cells.
#'
#' The two presence rasters must each have one layer and different layer names.
#' Those names become the species labels. Every cell with value `1` is kept, so
#' large rasters can produce a large table. Use [prepareDataNicheOverlap2()]
#' when the next step is an ecospat grid, which works from a sample of cells.
#'
#' @param environment A terra `SpatRaster` of environmental predictors.
#' @param rasterSPA A single-layer terra `SpatRaster`. Cells with value `1` are
#'   presences of the first species. The layer name is used as the species name.
#' @param rasterSPB A single-layer terra `SpatRaster` for the second species.
#'   Its layer name must differ from `rasterSPA`.
#'
#' @return A data frame. Column `species` identifies the presence raster, and
#'   the remaining columns are the environmental values at those cells. The
#'   mean environment of each species is printed.
#' @export
#' @importFrom terra as.data.frame extract nlyr
#' @importFrom stats aggregate complete.cases
#' @seealso [prepareDataNicheOverlap2()], [meanNicheOverlap()], [plotNicheOverlap()], [nicheSize()]
#'
#' @examples
#' \dontrun{
#' overlapData <- prepareDataNicheOverlap(environment, rasterSPA, rasterSPB)
#' meanNicheOverlap(overlapData, nsamples = 1000, colors.sp = c("red", "blue"))
#' }
prepareDataNicheOverlap <- function(environment, rasterSPA, rasterSPB) {

  if (terra::nlyr(rasterSPA) != 1 || terra::nlyr(rasterSPB) != 1) {
    stop("Error :: rasterSPA and rasterSPB must be single-layer presence rasters")
  }

  nameA <- names(rasterSPA)
  nameB <- names(rasterSPB)
  if (identical(nameA, nameB)) {
    stop("Error :: rasterSPA and rasterSPB must have different layer names, because those names are used as species names")
  }

  environment.used.by.species.1 <- .environmentAtPresences(environment, rasterSPA)
  environment.used.by.species.2 <- .environmentAtPresences(environment, rasterSPB)

  data <- cbind(
    data.frame(
      species = c(
        rep(nameA, nrow(environment.used.by.species.1)),
        rep(nameB, nrow(environment.used.by.species.2))
      )
    ),
    rbind(environment.used.by.species.1, environment.used.by.species.2)
  )
  rownames(data) <- NULL

  print(stats::aggregate(data[, -1, drop = FALSE], by = data["species"], FUN = mean))
  data

}
