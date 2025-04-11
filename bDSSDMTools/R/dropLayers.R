#' Drop Raster Layers Based on Variance Threshold
#'
#' Removes layers from a SpatRaster object that have a variance of cell values
#' below (or equal to) a specified threshold. This is useful for removing constant
#' or near-constant predictor layers before modeling. Variance is calculated
#' across all non-NA cells for each layer.
#'
#' @param rasters A terra `SpatRaster` object with one or more layers.
#' @param threshold Numeric. The minimum variance (exclusive) required for a layer
#'   to be kept. Layers with variance less than or equal to this threshold,
#'   or with undefined variance (e.g., fewer than 2 non-NA cells, constant values),
#'   will be dropped.
#'
#' @return A `terra SpatRaster` containing only the layers whose variance was
#'   strictly greater than the threshold. If no layers meet the criteria,
#'   an error might occur depending on subsequent usage (consider checking
#'   `nlyr` of the result).
#' @export
#' @importFrom terra global subset SpatRaster nlyr values
#' @importFrom stats var # Although used via fun="var", importing is safer practice
#'
#' @examples
#' \dontrun{
#' library(terra)
#' # Create example rasters
#' r1 <- rast(nrows=10, ncols=10, vals=rnorm(100)) # Variable layer
#' r2 <- rast(nrows=10, ncols=10, vals=1)          # Constant layer
#' r3 <- rast(nrows=10, ncols=10, vals=runif(100)/100) # Low variance layer
#' r4 <- rast(nrows=10, ncols=10, vals=c(NA, 1:99)) # Layer with NAs
#' names(r1) <- "Variable"
#' names(r2) <- "Constant"
#' names(r3) <- "LowVariance"
#' names(r4) <- "WithNAs"
#'
#' rasters_combined <- c(r1, r2, r3, r4)
#' print(rasters_combined)
#'
#' # Calculate variances (using terra::global)
#' print(global(rasters_combined, fun="var", na.rm=TRUE))
#'
#' # Drop layers with variance <= 0.01
#' rasters_filtered <- dropLayers(rasters_combined, threshold = 0.01)
#' print(rasters_filtered) # Should keep Variable and WithNAs
#'
#' # Drop layers with variance <= 0
#' rasters_filtered_zero <- dropLayers(rasters_combined, threshold = 0)
#' print(rasters_filtered_zero) # Should keep Variable, LowVariance, WithNAs
#' }
dropLayers <- function(rasters, threshold) {

  # --- Input Checks ---
  if (!methods::is(rasters, "SpatRaster")) {
    stop("Error: 'rasters' must be a SpatRaster object.")
  }
  if (!is.numeric(threshold) || length(threshold) != 1) {
    stop("Error: 'threshold' must be a single numeric value.")
  }
  if (terra::nlyr(rasters) == 0) {
    warning("Input SpatRaster has no layers.")
    return(rasters)
  }

  # --- Calculate Variance per Layer ---
  # Use terra::global for efficiency, avoids creating large data frame
  # Ask for variance, removing NAs during calculation
  message("Calculating variance for each layer...")
  layer_stats <- terra::global(rasters, fun = "var", na.rm = TRUE)

  # The result is a data frame, the column name is 'var'
  variances <- layer_stats$var

  # --- Identify Layers to Keep ---
  # Keep layers where variance is not NA and is strictly greater than the threshold
  indices_to_keep <- which(!is.na(variances) & variances > threshold)

  if (length(indices_to_keep) == 0) {
    warning("No layers found with variance > ", threshold, ". Returning an empty SpatRaster (or original if only 1 layer).")
    # Return an empty raster with same structure or handle as needed
    # For safety, let subset handle the empty case if nlyr > 1
    if (terra::nlyr(rasters) == 1 && length(indices_to_keep)==0) return(rasters[[0L]]) # Explicitly empty
    if (terra::nlyr(rasters) == 0) return(rasters) # Already empty
  } else if (length(indices_to_keep) == terra::nlyr(rasters)) {
    message("All layers have variance > ", threshold, ". Returning original raster.")
    return(rasters)
  }

  message("Keeping ", length(indices_to_keep), " out of ", terra::nlyr(rasters), " layers.")

  # --- Subset the Raster ---
  # Use the numeric indices to select layers
  rasters_subset <- terra::subset(rasters, indices_to_keep)

  return(rasters_subset)

}
