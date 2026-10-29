#' Download Multiple Bio-Oracle Layers
#'
#' Downloads Bio-ORACLE layers for one experiment and one decade. Each
#' requested layer is a variable, a statistic, and a depth realm.
#'
#' `variables` can be any of:
#' - a list of length-3 character vectors, `c(variable, predictor, realm)`
#' - one such vector
#' - a table with three columns, in that order
#' - a table with three rows (variable, predictor, realm) and one column per layer
#'
#' Predictors are `ltmax`, `ltmin`, `mean`, `max`, `min`, and `range`.
#' Realms are `depthsurf`, `depthmax`, `depthmean`, and `depthmin`.
#' Experiments are `baseline`, `ssp119`, `ssp126`, `ssp245`, `ssp370`,
#' `ssp460`, and `ssp585`.
#'
#' @param variables Layer requests, as described above.
#' @param experiment One Bio-ORACLE experiment name.
#' @param decade One decade, such as `2010` or `2050`.
#' @param latitude Optional numeric vector of length 2, minimum and maximum.
#' @param longitude Optional numeric vector of length 2, minimum and maximum.
#' @param exportRaster If `TRUE` (the default), return the downloaded layers
#'   as one SpatRaster. If `FALSE`, return the paths of the GeoTIFF files.
#' @param directory Optional folder in which to also save the GeoTIFF files.
#'
#' @return A SpatRaster when `exportRaster` is `TRUE`, or a character vector
#'   of file paths when it is `FALSE`.
#' @export
#' @importFrom terra rast writeRaster
#'
download_multiple_layers <- function(variables, experiment, decade, latitude = NULL, longitude = NULL, exportRaster = TRUE, directory = NULL) {

  if (is.null(longitude)) longitude <- c(-179.975, 179.975)
  if (is.null(latitude)) latitude <- c(-89.975, 89.975)
  if (length(longitude) != 2) stop("Error :: longitude must have two elements, minimum and maximum")
  if (length(latitude) != 2) stop("Error :: latitude must have two elements, minimum and maximum")
  if (length(experiment) != 1) stop("Error :: experiment must have one element")
  if (length(decade) != 1) stop("Error :: decade must have one element")
  if (!is.null(directory) && !dir.exists(directory)) {
    stop("Error :: directory does not exist: ", directory)
  }

  requests <- .biooracleRequests(variables)
  requests$experiment <- as.character(experiment)
  requests$decade <- decade

  realms <- c("depthsurf", "depthmax", "depthmean", "depthmin")
  experiments <- c("baseline", "ssp119", "ssp126", "ssp245", "ssp370", "ssp460", "ssp585")
  predictors <- c("ltmax", "ltmin", "mean", "max", "min", "range")
  badRealm <- unique(requests$realm[!requests$realm %in% realms])
  badExperiment <- unique(requests$experiment[!requests$experiment %in% experiments])
  badPredictor <- unique(requests$predictor[!requests$predictor %in% predictors])
  if (length(badRealm) > 0) stop("Error :: unrecognized realm: ", paste(badRealm, collapse = ", "), ". Use ", paste(realms, collapse = ", "))
  if (length(badExperiment) > 0) stop("Error :: unrecognized experiment: ", paste(badExperiment, collapse = ", "), ". Use ", paste(experiments, collapse = ", "))
  if (length(badPredictor) > 0) stop("Error :: unrecognized predictor: ", paste(badPredictor, collapse = ", "), ". Use ", paste(predictors, collapse = ", "))

  datasetLayers <- biooracler::list_layers()
  combLayers <- NULL
  written <- character()
  skipped <- character()

  for (i in seq_len(nrow(requests))) {

    variable <- requests$variable[i]
    predictor <- requests$predictor[i]
    realm <- requests$realm[i]
    experiment.i <- requests$experiment[i]
    decade.i <- requests$decade[i]
    requestLabel <- paste(variable, predictor, realm)

    if (variable %in% c("kdpar_mean", "par_mean")) {
      variable_string <- paste0(variable, "_")
    } else {
      variable_string <- paste0(variable, "_")
    }

    datasetLayers.i <- datasetLayers[grepl(variable_string, datasetLayers$dataset_id), ]
    datasetLayers.i <- datasetLayers.i[grepl(experiment.i, datasetLayers.i$dataset_id), ]
    datasetLayers.i <- datasetLayers.i[grepl(realm, datasetLayers.i$dataset_id), ]

    if (nrow(datasetLayers.i) == 0) {
      skipped <- c(skipped, paste0(requestLabel, " (no matching dataset)"))
      next
    }

    if (variable == "par_mean") {
      datasetLayers.i <- datasetLayers.i[!grepl("kdpar_mean", datasetLayers.i$dataset_id), ]
    }
    if (variable == "kdpar_mean") {
      datasetLayers.i <- datasetLayers.i[grepl("kdpar_mean", datasetLayers.i$dataset_id), ]
    }

    if (nrow(datasetLayers.i) != 1) {
      stop(
        "Error :: expected one Bio-ORACLE dataset for ", requestLabel,
        " and experiment ", experiment.i, ", but found ", nrow(datasetLayers.i)
      )
    }

    var_name <- datasetLayers.i$dataset_id
    var_name <- gsub("kdpar_mean", "kdparmean", var_name)
    var_name <- gsub("par_mean", "parmean", var_name)
    time <- unlist(strsplit(var_name, split = "_"))
    time <- as.numeric(time[3:4])
    time <- floor(time / 10 + 0.5) * 10
    time <- seq(time[1], time[2], by = 10)
    time <- time[-length(time)]

    if (!decade.i %in% time) {
      skipped <- c(skipped, paste0(requestLabel, " (decade ", decade.i, " is outside ", min(time), "-", max(time), ")"))
      next
    }

    timeWindow <- rep(paste0(decade.i, "-01-01T00:00:00Z"), 2)
    constraints <- list(time = timeWindow, latitude = latitude, longitude = longitude)
    myDir <- file.path(tempdir(), paste0("biooracle-", sample.int(1e9, 1)))
    dir.create(myDir, recursive = TRUE)

    predictor_var <- paste0(variable, "_", predictor)

    biooracler::download_layers(
      dataset_id = datasetLayers.i$dataset_id,
      variables = predictor_var,
      constraints = constraints,
      directory = myDir
    )

    fileDownloaded <- list.files(myDir, pattern = "\\.nc$", full.names = TRUE)
    if (length(fileDownloaded) == 0) {
      stop("Error :: Bio-ORACLE did not write a NetCDF file for ", requestLabel)
    }
    fileDownloaded <- fileDownloaded[which.max(file.info(fileDownloaded)$mtime)]

    filename <- datasetLayers.i$title
    filename <- gsub("Bio-Oracle ", "", filename)
    filename <- gsub("\\.", "", filename)
    filename <- sub(paste0(sub(".* ", "", filename), "$"), "", filename)
    filename <- paste0(filename, decade.i, " ", predictor)
    filename <- gsub("\\[|\\]", "", filename)
    filename <- gsub(" ", "_", filename)

    raster <- terra::rast(fileDownloaded)
    names(raster) <- filename
    tifPath <- file.path(myDir, paste0(filename, ".tif"))
    terra::writeRaster(raster, filename = tifPath, overwrite = TRUE)

    if (!is.null(directory)) {
      savedPath <- file.path(directory, paste0(filename, ".tif"))
      terra::writeRaster(raster, filename = savedPath, overwrite = TRUE)
      written <- c(written, savedPath)
    } else {
      written <- c(written, tifPath)
    }

    raster <- terra::rast(tifPath)
    names(raster) <- filename
    if (is.null(combLayers)) combLayers <- raster else combLayers <- c(combLayers, raster)
    file.remove(fileDownloaded)

  }

  if (is.null(combLayers)) {
    stop("Error :: no Bio-ORACLE layers were downloaded. ", paste(skipped, collapse = "; "))
  }
  if (length(skipped) > 0) {
    message("Skipped layers: ", paste(skipped, collapse = "; "))
  }
  if (isTRUE(exportRaster)) combLayers else written

}

.biooracleRequests <- function(variables) {

  if (is.character(variables) && is.null(dim(variables)) && length(variables) == 3) {
    requests <- data.frame(
      variable = variables[1],
      predictor = variables[2],
      realm = variables[3],
      stringsAsFactors = FALSE
    )
    return(requests)
  }

  if (is.list(variables) && !is.data.frame(variables)) {
    lengths <- vapply(variables, length, integer(1))
    if (length(variables) == 0 || any(lengths != 3)) {
      stop("Error :: each layer in variables must be c(variable, predictor, realm)")
    }
    requests <- as.data.frame(do.call(rbind, lapply(variables, as.character)), stringsAsFactors = FALSE)
    names(requests) <- c("variable", "predictor", "realm")
    return(requests)
  }

  if (is.null(variables) || (!is.data.frame(variables) && !is.matrix(variables))) {
    stop("Error :: variables must be a list of c(variable, predictor, realm), or a table with those three fields")
  }

  vars <- as.data.frame(variables, stringsAsFactors = FALSE)
  if (nrow(vars) == 3) {
    requests <- as.data.frame(t(vars), stringsAsFactors = FALSE)
  } else if (ncol(vars) >= 3) {
    requests <- vars[, 1:3, drop = FALSE]
  } else {
    stop("Error :: variables must have three fields: variable, predictor, and realm")
  }
  names(requests) <- c("variable", "predictor", "realm")
  requests[] <- lapply(requests, as.character)
  requests

}
