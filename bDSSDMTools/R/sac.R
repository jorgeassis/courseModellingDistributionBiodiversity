#' Estimate Spatial Autocorrelation Range and Suggest Block Size
#'
#' Estimates the distance at which each environmental predictor is no longer
#' spatially autocorrelated, and suggests a block size for spatial
#' cross-validation. A random sample of cells is drawn from longitude-latitude
#' layers. Within geographic distance classes, the correlation between
#' geographic distance and environmental difference is tested. The range of a
#' predictor is the first class, in kilometres, whose correlation is no longer
#' significant.
#'
#' The sampling uses a fixed seed internally and restores the session random
#' seed afterwards.
#'
#' @param environmentalLayers A terra `SpatRaster` in longitude-latitude
#'   coordinates, with one or more predictor layers.
#' @param num_sample Number of cells to sample. Defaults to 2500.
#' @param autocorrelationClassDistance Width of each distance class, in
#'   kilometres. Defaults to 10.
#' @param autocorrelationMaxDistance Largest distance tested, in kilometres.
#'   A predictor whose correlation never becomes non-significant is given this
#'   range. Defaults to 500.
#' @param autocorrelationSignif Significance level used to decide that the
#'   correlation within a distance class is no longer significant. Defaults
#'   to 0.05.
#'
#' @return A list with:
#'   \item{plot}{A ggplot of the range of each predictor and the mean block size.}
#'   \item{distance}{Suggested block size in kilometres, the mean of the predictor ranges.}
#' @export
#' @importFrom terra spatSample is.lonlat
#' @importFrom sp spDists
#' @importFrom stats lm pf cor reorder
#' @importFrom ggplot2 ggplot aes geom_segment geom_point geom_vline annotate labs theme_minimal
#'
sac <- function(environmentalLayers, num_sample = 2500, autocorrelationClassDistance = 10, autocorrelationMaxDistance = 500, autocorrelationSignif = 0.05) {

  if (!inherits(environmentalLayers, "SpatRaster")) {
    stop("Error :: environmentalLayers must be a SpatRaster")
  }
  if (!isTRUE(terra::is.lonlat(environmentalLayers))) {
    stop("Error :: sac expects environmental layers in longitude-latitude coordinates so ranges are in kilometres")
  }
  if (terra::nlyr(environmentalLayers) < 1) {
    stop("Error :: environmentalLayers must contain at least one layer")
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
  set.seed(42)

  environmentalLayers.r <- terra::spatSample(
    environmentalLayers,
    size = num_sample,
    method = "random",
    na.rm = TRUE,
    xy = TRUE,
    as.df = TRUE,
    warn = FALSE
  )
  if (nrow(environmentalLayers.r) < 3) {
    stop("Error :: not enough complete cells to estimate spatial autocorrelation")
  }

  space <- as.vector(sp::spDists(
    as.matrix(environmentalLayers.r[, c("x", "y")]),
    as.matrix(environmentalLayers.r[, c("x", "y")]),
    longlat = TRUE
  ))
  predictors <- environmentalLayers.r[, !names(environmentalLayers.r) %in% c("x", "y"), drop = FALSE]
  classFrom <- seq(0, autocorrelationMaxDistance - autocorrelationClassDistance, by = autocorrelationClassDistance)
  classTo <- seq(autocorrelationClassDistance, autocorrelationMaxDistance, by = autocorrelationClassDistance)

  range <- numeric(ncol(predictors))

  for (variableIndex in seq_len(ncol(predictors))) {

    envDistance <- as.vector(sp::spDists(
      as.matrix(predictors[, variableIndex, drop = FALSE]),
      as.matrix(predictors[, variableIndex, drop = FALSE]),
      longlat = FALSE
    ))
    pValue <- rep(NA_real_, length(classFrom))

    for (bin in seq_along(classFrom)) {

      keep <- which(space >= classFrom[bin] & space <= classTo[bin])
      if (length(keep) < 3) next

      fit <- stats::lm(space[keep] ~ envDistance[keep])
      fstat <- summary(fit)$fstatistic
      if (!is.null(fstat) && length(fstat) == 3) {
        pValue[bin] <- stats::pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
      }

    }

    nonsignificant <- which(pValue >= autocorrelationSignif)
    range[variableIndex] <- if (length(nonsignificant) == 0) {
      autocorrelationMaxDistance
    } else {
      round(classTo[nonsignificant[1]])
    }

  }

  df <- data.frame(
    variable = colnames(predictors),
    range = as.numeric(range),
    stringsAsFactors = FALSE
  )

  plot1 <- ggplot2::ggplot(df, ggplot2::aes(x = range, y = stats::reorder(variable, range))) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = range, yend = variable), color = "#5aa0d6", linewidth = 1.2) +
    ggplot2::geom_point(size = 3, color = "#5aa0d6") +
    ggplot2::geom_vline(xintercept = mean(range), linetype = "dashed", color = "red") +
    ggplot2::annotate(
      "text", x = mean(range), y = 1, label = "Block size",
      angle = 90, vjust = 1.5, color = "red", size = 4
    ) +
    ggplot2::labs(
      title = "Spatial autocorrelation range",
      subtitle = "First distance class whose correlation is no longer significant",
      x = "Range (km)",
      y = "Variable"
    ) +
    ggplot2::theme_minimal()

  list(plot = plot1, distance = mean(range))

}
