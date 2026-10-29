#' Generate Correlation Plot and Matrix for Raster Layers
#'
#' Calculates and plots the correlation matrix for a sample of cell values
#' from a multi-layer SpatRaster. Cells are sampled directly from the raster,
#' up to 5000 complete cells, so the full raster is not copied into memory.
#'
#' @param rasterLayers A terra `SpatRaster` object with at least two layers.
#'
#' @return A list containing:
#'   \item{plot}{A recorded base-graphics plot. Draw it again with `print()`.}
#'   \item{correlationPairs}{The correlation matrix calculated with `stats::cor`.}
#' @export
#' @importFrom terra spatSample nlyr
#' @importFrom stats cor
#' @importFrom corrplot corrplot
#'
#' @examples
#' \dontrun{
#' r1 <- terra::rast(nrows = 10, ncols = 10, vals = rnorm(100))
#' r2 <- r1 * 1.5 + rnorm(100, 0, 0.5)
#' my_rasters <- c(r1, r2)
#' names(my_rasters) <- c("Var1", "Var2")
#' corr_info <- getCorrPlot(my_rasters)
#' print(corr_info$plot)
#' print(corr_info$correlationPairs)
#' }
getCorrPlot <- function(rasterLayers) {

  if (!inherits(rasterLayers, "SpatRaster")) {
    stop("Error :: rasterLayers must be a SpatRaster")
  }
  if (terra::nlyr(rasterLayers) < 2) {
    stop("Error :: at least two layers are required to calculate correlations")
  }

  sampled <- terra::spatSample(
    rasterLayers,
    size = min(5000, terra::ncell(rasterLayers)),
    method = "random",
    na.rm = TRUE,
    as.df = TRUE,
    warn = FALSE
  )
  if (nrow(sampled) < 3) {
    stop("Error :: not enough complete cells to calculate correlations")
  }

  cor_matrix <- stats::cor(sampled)
  if (any(!is.finite(cor_matrix))) {
    stop("Error :: correlations could not be calculated. A sampled layer may have no variance.")
  }

  openedDevice <- grDevices::dev.cur() == 1
  plotFile <- tempfile(fileext = ".pdf")
  if (openedDevice) grDevices::pdf(plotFile)
  on.exit({
    if (openedDevice && grDevices::dev.cur() > 1) grDevices::dev.off()
    if (openedDevice) unlink(plotFile)
  }, add = TRUE)

  corrplot::corrplot(cor_matrix, type = "upper", diag = FALSE)
  recorded <- grDevices::recordPlot()

  list(plot = recorded, correlationPairs = cor_matrix)

}
