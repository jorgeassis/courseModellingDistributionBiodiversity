#' Download Multiple Bio-Oracle Layers
#'
#' Downloads specified Bio-Oracle layers based on combinations of variables,
#' predictors, experiments, realms, and decades. Allows spatial subsetting
#' and optionally returns a combined multi-layer SpatRaster.
#'
#' @export
#'
#' @importFrom terra rast writeRaster c SpatRaster ext rowColFromCell xyFromCell rasterize buffer values `values<-` `names<-` crop trim mask classify ncell nlyr crs `crs<-` project units `units<-` wrap unwrap deepcopy init `res<-` relate is.related centroids values as.data.frame vect crds geomtype subset maths Compare Logic Arith `Arith,SpatRaster,SpatRaster-method` `Arith,SpatRaster,numeric-method` `Arith,numeric,SpatRaster-method` `Compare,SpatRaster,numeric-method` `Compare,SpatRaster,SpatRaster-method` `Compare,numeric,SpatRaster-method` `Logic,SpatRaster,logical-method` `Logic,SpatRaster,SpatRaster-method` `Logic,!SpatRaster-method` subset extract `subset,SpatRaster,numeric-method` `subset,SpatRaster,missing-method` `extract,SpatRaster,SpatVector-method` `extract,SpatRaster,matrix-method` `is.na,SpatRaster-method` `anyNA,SpatRaster-method`
#' @importFrom utils file.info tempdir

download_multiple_layers <- function(variables, experiment, decade, latitude=NULL, longitude=NULL, exportRaster=TRUE, directory=NULL) {

  if( is.null(longitude) ) {
    longitude <- c(-179.975, 179.975)
  }
  if( is.null(latitude) ) {
    latitude <- c(-89.975, 89.975)
  }

  if( length(longitude) != 2 ) {
    stop("Error :: Longitude must have two elements, minimum and maximum values")
  }
  if( length(latitude) != 2 ) {
    stop("Error :: Latitude must have two elements, minimum and maximum values")
  }
  if( length(experiment) != 1 ) {
    stop("Error :: Experiment must have one element")
  }
  if( length(decade) != 1 ) {
    stop("Error :: Decade must have one element")
  }

  # -----------------------

  combinations <- data.frame(t(data.frame(variables)),experiment, decade, stringsAsFactors = FALSE)
  combLayers <- NULL

  if( sum(! combinations[,3] %in% c("depthsurf","depthmax","depthmean","depthmin")) > 0) { stop("Unrecognized realm in list") }
  if( sum(! combinations[,4] %in% c("baseline","ssp119","ssp126","ssp245","ssp370","ssp460","ssp585")) > 0 ) { stop("Unrecognized experiment in list") }
  if( sum(! combinations[,2] %in% c("ltmax","ltmin","mean","max","min","range")) > 0 ) { stop("Unrecognized predictor in list") }

  # -----------------------

  for(i in 1:nrow(combinations)) {

    variable <- combinations[i,1]
    predictor <- combinations[i,2]
    realm <- combinations[i,3]
    experiment <- combinations[i,4]
    decade <- combinations[i,5]

    if( variable == "kdpar_mean" ) {
      variable_string <- paste0("kdpar_mean_")
    }
    if( variable == "par_mean" ) {
      variable_string <- paste0("par_mean_")
    }
    if( variable != "kdpar_mean" & variable != "par_mean" ) {
      variable_string <- paste0(variable,"_")
    }

    datasetLayers <- biooracler::list_layers()
    datasetLayers.i <- datasetLayers[grepl(variable_string, datasetLayers$dataset_id),]
    datasetLayers.i <- datasetLayers.i[grepl(experiment, datasetLayers.i$dataset_id),]
    datasetLayers.i <- datasetLayers.i[grepl(realm, datasetLayers.i$dataset_id),]

    if( nrow(datasetLayers.i) == 0 ) { next }

    if( nrow(datasetLayers.i) != 1 ) {

      if( variable == "par_mean" ) {
        datasetLayers.i <- datasetLayers.i[!grepl("kdpar_mean", datasetLayers.i$dataset_id),]
      }
      if( variable == "kdpar_mean" ) {
        datasetLayers.i <- datasetLayers.i[grepl("kdpar_mean", datasetLayers.i$dataset_id),]
      }

    }

    if( nrow(datasetLayers.i) != 1 ) { stop("Error :: 001") }

    var_name <- datasetLayers.i$dataset_id
    var_name <- gsub("kdpar_mean","kdparmean",var_name)
    var_name <- gsub("par_mean","parmean",var_name)

    time <- unlist(strsplit(var_name, split = "_"))
    time <- as.numeric(time[3:4])
    time <- floor(time / 10 + 0.5) * 10
    time <- seq(time[1],time[2],by=10)
    time <- time[-length(time)]

    if( ! decade %in% time ) { next }

    timeWindow <- rep(paste0(decade,"-01-01T00:00:00Z") ,2)
    constraints <- list(time = timeWindow , latitude = latitude,longitude = longitude)

    myDir <- paste0(tempdir(),"/",sample(1:100000000000,1))
    if( ! dir.exists(myDir) ) { dir.create(myDir, recursive = TRUE) }

    predictor_var <- paste0(sub("_.*", "", datasetLayers.i$dataset_id),"_",predictor)

    if( grepl("kdpar", predictor_var) ) {
      predictor_var <- gsub("kdpar","kdpar_mean", predictor_var)
    }
    if( grepl("parmean", predictor_var) ) {
      predictor_var <- gsub("par","par_mean", predictor_var)
    }

    biooracler::download_layers(dataset_id=datasetLayers.i$dataset_id,
                                variables = predictor_var,
                                constraints = constraints,
                                directory= myDir )

    fileDownloaded <- list.files(myDir, pattern = ".nc", full.names = TRUE)
    fileDownloadedDetails <- file.info(fileDownloaded)
    fileDownloaded <- fileDownloaded[which.max(fileDownloadedDetails$mtime)]

    filename <- datasetLayers.i$title
    filename <- gsub("Bio-Oracle ","",filename)
    filename <- gsub("\\.","",filename)

    filename <- sub(paste0(sub(".* ", "", filename), "$"), "", filename, fixed = FALSE)
    filename <- paste0(filename,decade," ",predictor)

    pattern_to_remove <- "\\[|\\]"
    filename <- gsub(pattern_to_remove, "", filename)

    filename <- gsub(" ","_",filename)

    raster <- terra::rast(fileDownloaded)
    names(raster) <- filename

    writeRaster(raster, filename = paste0(myDir,"/",filename,".tif"), overwrite=TRUE)

    if( ! is.null(directory) ) {
      writeRaster(raster, filename = paste0(directory,"/",filename,".tif"), overwrite=TRUE)
    }

    if(exportRaster) {
      raster <- rast(paste0(myDir,"/",filename,".tif"))
      names(raster) <- gsub(" ","_",names(raster))
      if( is.null(combLayers) ) { combLayers <- raster } else { combLayers <- c(combLayers, raster) }
    }

    file.remove(fileDownloaded)

  }

  return(combLayers)

}

