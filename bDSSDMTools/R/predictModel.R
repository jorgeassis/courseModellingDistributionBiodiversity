#' Predict Ensemble SDM Results onto Raster Layers
#'
#' Generates spatial predictions from an ensemble of SDM models (the list
#' returned by [trainModel()]). Each model is applied with `terra::predict`
#' in blocks, so the raster is not loaded into a single data frame. The
#' result is the mean probability across models and the standard deviation
#' of those probabilities. A single model has no standard deviation.
#'
#' @param model The list returned by [trainModel()].
#' @param newData A terra `SpatRaster` of predictor layers. Layer names must
#'   include the predictors used to fit each model.
#'
#' @return A list with `rasterLayer` (mean probability) and `rasterLayerSD`.
#' @export
#' @importFrom terra predict rast mean stdev
#'
predictModel <- function(model, newData) {

  if (!inherits(newData, "SpatRaster")) {
    stop("Error :: newData must be a SpatRaster")
  }
  if (!is.list(model) || is.null(model$model) || length(model$model) == 0) {
    stop("Error :: model must be the list returned by trainModel, with at least one fitted model")
  }

  predictorNames <- function(model.j) {
    cls <- class(model.j)
    if ("glm" %in% cls) return(attr(stats::terms(model.j), "term.labels"))
    if ("gbm" %in% cls && !is.null(model.j$var.names)) return(model.j$var.names)
    if ("maxnet" %in% cls) {
      if (!is.null(model.j$varnames)) return(model.j$varnames)
      if (!is.null(model.j$samplemeans)) return(names(model.j$samplemeans))
    }
    if (any(c("xgb.Booster", "xgboost") %in% cls) && !is.null(model.j$feature_names)) return(model.j$feature_names)
    NULL
  }

  predictLayer <- function(model.j) {
    needed <- predictorNames(model.j)
    predictors <- newData
    if (!is.null(needed)) {
      missingPredictors <- setdiff(needed, names(newData))
      if (length(missingPredictors) > 0) {
        stop(
          "Error :: newData is missing predictors required by a ",
          paste(class(model.j), collapse = "/"),
          " model: ",
          paste(missingPredictors, collapse = ", ")
        )
      }
      predictors <- newData[[needed]]
    }

    cls <- class(model.j)
    if (any(cls %in% c("gbm", "brt"))) {
      terra::predict(predictors, model.j, na.rm = TRUE, fun = function(mod, data) {
        as.numeric(stats::predict(mod, data, n.trees = mod$n.trees, type = "response"))
      })
    } else if ("maxnet" %in% cls) {
      terra::predict(predictors, model.j, na.rm = TRUE, fun = function(mod, data) {
        as.numeric(stats::predict(mod, data, type = "cloglog"))
      })
    } else if ("glm" %in% cls) {
      terra::predict(predictors, model.j, na.rm = TRUE, fun = function(mod, data) {
        as.numeric(stats::predict(mod, data, type = "response"))
      })
    } else if (any(c("xgb.Booster", "xgboost") %in% cls)) {
      terra::predict(predictors, model.j, na.rm = TRUE, fun = function(mod, data) {
        as.numeric(stats::predict(mod, as.matrix(data)))
      })
    } else {
      stop("Error :: predictions are not implemented for models of class ", paste(cls, collapse = ", "))
    }
  }

  stacked <- terra::rast(lapply(model$model, predictLayer))
  prediction <- terra::mean(stacked)
  if (terra::nlyr(stacked) == 1) {
    predictionSD <- prediction
    terra::values(predictionSD) <- NA_real_
  } else {
    predictionSD <- terra::stdev(stacked)
  }
  names(prediction) <- "prediction"
  names(predictionSD) <- "sd"
  terra::units(prediction) <- "mean probability"
  terra::units(predictionSD) <- "sd of the mean probability"

  list(rasterLayer = prediction, rasterLayerSD = predictionSD)

}
