# Internal helpers for the niche-overlap functions. Not exported.

.speciesLabels <- function(data) {

  if (!is.data.frame(data) || ncol(data) < 1 || names(data)[1] != "species") {
    stop("Error :: the first column of data must be named species, followed by environmental variables")
  }

  species <- data$species
  if (anyNA(species)) {
    stop("Error :: species contains missing values")
  }
  if (is.factor(species)) {
    droplevels(species)
  } else {
    as.character(species)
  }

}

.nichePosterior <- function(data, nsamples) {

  if (!is.data.frame(data)) {
    stop("Error :: data must be a data.frame")
  }
  if (ncol(data) < 2) {
    stop("Error :: the first column of data must be named species, followed by environmental variables")
  }
  if (!is.numeric(nsamples) || length(nsamples) != 1 || !is.finite(nsamples) || nsamples < 1) {
    stop("Error :: nsamples must be a positive number")
  }

  species <- .speciesLabels(data)

  predictors <- data[, -1, drop = FALSE]
  if (!all(vapply(predictors, is.numeric, logical(1)))) {
    stop("Error :: environmental columns must be numeric")
  }
  if (anyNA(predictors)) {
    stop("Error :: environmental columns contain missing values")
  }

  counts <- table(species)
  if (any(counts <= ncol(predictors))) {
    stop("Error :: each species needs more observations than environmental variables")
  }

  tapply(seq_len(nrow(data)), species, function(ii) {
    nicheROVER::niw.post(nsamples = nsamples, X = predictors[ii, , drop = FALSE])
  })

}

.checkSpeciesColors <- function(colors, species.par) {

  speciesNames <- names(species.par)
  if (length(colors) != length(speciesNames)) {
    stop(
      "Error :: provide one colour per species, in this order: ",
      paste(speciesNames, collapse = ", ")
    )
  }
  colors

}

.environmentAtPresences <- function(environment, presenceRaster) {

  if (!inherits(environment, "SpatRaster")) {
    stop("Error :: environment must be a terra SpatRaster")
  }
  if (!inherits(presenceRaster, "SpatRaster") || terra::nlyr(presenceRaster) != 1) {
    stop("Error :: each presence raster must be a single-layer SpatRaster")
  }

  cells <- terra::as.data.frame(presenceRaster, xy = TRUE)
  layerName <- names(presenceRaster)
  present <- which(cells[[layerName]] == 1)
  if (length(present) == 0) {
    stop("Error :: no presence cells (value 1) in layer ", layerName)
  }

  values <- terra::extract(environment, cells[present, c("x", "y"), drop = FALSE], ID = FALSE)
  values <- values[stats::complete.cases(values), , drop = FALSE]
  if (nrow(values) == 0) {
    stop("Error :: no complete environmental values at presence cells of ", layerName)
  }
  if (nrow(values) <= ncol(values)) {
    stop("Error :: layer ", layerName, " needs more presence cells with complete data than environmental variables")
  }

  values

}

.sampleRows <- function(x, maxPoints) {

  n <- nrow(x)
  nKeep <- min(n, maxPoints)
  if (nKeep >= n) {
    return(x)
  }
  x[sample.int(n, nKeep), , drop = FALSE]

}
