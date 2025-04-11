#' Predict Ensemble SDM Results onto Raster Layers
#'
#' Generates spatial predictions from an ensemble of SDM models (typically the
#' output list from `trainModel`). It applies each individual model to the new
#' environmental data (`newData` SpatRaster) using `terra::predict`. It then
#' calculates the mean and standard deviation across all model predictions for
#' each grid cell, providing an ensemble prediction and an uncertainty estimate.
#' This approach avoids loading large rasters into memory.
#'
#' @export
#' @importFrom terra SpatRaster predict rast c mean stdev nlyr `names<-` `units<-` crs ext res `is.na<-` project maths Compare Logic Arith `Arith,SpatRaster,SpatRaster-method` `Arith,SpatRaster,numeric-method` `Arith,numeric,SpatRaster-method` `Compare,SpatRaster,numeric-method` `Compare,SpatRaster,SpatRaster-method` `Compare,numeric,SpatRaster-method` `Logic,SpatRaster,logical-method` `Logic,SpatRaster,SpatRaster-method` `Logic,!SpatRaster-method` subset extract mask classify crop trim values `values<-` xyFromCell vect is.projected ncell `[` `anyNA,SpatRaster-method`
#' @importFrom stats predict sd # Generic predict and sd
#' @importFrom gbm predict.gbm # Ensure predict method is available for terra::predict
#' @importFrom maxnet predict.maxnet # Ensure predict method is available for terra::predict
#' @importFrom xgboost predict.xgb.Booster xgb.DMatrix # Ensure predict method is available
#' @importFrom methods is
#'
#'
predictModel <- function(model,newData) {

  shape <- newData[[1]]
  shape[!is.na(shape)] <- 1
  cells <- which(!is.na( values(shape)[,1] ))

  newData <- as.data.frame(newData)

  for(j in 1:length(model$model) ) {

    model.j <- model$model[j][[1]]
    algorithmType <- class(model.j)[1]
    prediction <- NULL

    if( algorithmType == "brt" | algorithmType == "gbm" ) {
      n.trees <- model.j$n.trees
      prediction <- predict(model.j,newData,n.trees=n.trees, type="response")
    }

    if( algorithmType == "maxnet" ) {
      prediction <- predict(model.j,newData,type="cloglog")
    }

    if( algorithmType == "glm" ) {
      prediction <- predict(model.j,newData,type="response")
    }

    if(algorithmType == "mboost") {
      prediction <- predict(model.j,newData, type="response")[,1]
    }

    if(algorithmType == "xgb.Booster") {
      newData <- newData[model.j$feature_names]
      newData.x = xgb.DMatrix(data = data.matrix(newData[,setdiff(names(newData),c("PA"))]),  nthread = 1)
      prediction <- predict(model.j,newData.x, type="response")
    }

    if( algorithmType != "maxnet" & algorithmType != "glm" & algorithmType != "gbm" & algorithmType != "mboost" & algorithmType != "xgb.Booster"  ) {
      stop("Error :: 032")
    }

    if( j == 1 ) { predictedDistribution <- prediction } else { predictedDistribution <- cbind(predictedDistribution,prediction) }

  }

  predictedDistributionSD <- apply(predictedDistribution,1,sd,na.rm=T)
  predictedDistribution <- apply(predictedDistribution,1,mean,na.rm=T)

  shape[cells] <- predictedDistribution
  prediction <- shape
  names(prediction) <- "prediction"
  units(prediction) <- "mean probability"

  shape[cells] <- predictedDistributionSD
  predictionSD <- shape
  names(predictionSD) <- "sd"
  units(predictionSD) <- "sd of the mean probability"


  return(list(rasterLayer=prediction,rasterLayerSD=predictionSD))

}
