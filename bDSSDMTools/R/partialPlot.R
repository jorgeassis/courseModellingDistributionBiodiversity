#' Generate Partial Dependence Plot (Response Curve)
#'
#' Calculates and plots the partial dependence (response curve) for a specified
#' predictor variable from an ensemble of models. It shows how the predicted
#' probability, averaged across models, changes across the observed range of the
#' chosen variable, while holding other predictor variables constant at their mean
#' values. Predictions that are already probabilities are not transformed again.
#' The curve is not forced to be monotonic. Includes a rug plot showing the
#' distribution of the variable for presence points.
#'
#' @param model A model, or the list returned by [trainModel()].
#' @param modelData A data frame of records, or the list returned by
#'   [generateModelData()]. It must contain `PA` and `variablePlot`.
#' @param variablePlot Name of the predictor to vary. Other predictors are held
#'   at their mean.
#' @export
#' @importFrom stats predict quantile
#' @importFrom ggplot2 ggplot aes geom_point geom_line theme labs element_text element_rect element_line element_blank margin xlab ylab
#' @importFrom gbm predict.gbm
#' @importFrom maxnet maxnet
#' @importFrom xgboost xgboost
#'

partialPlot <- function(model,modelData,variablePlot) {

  themePlot <- theme(
    text = element_text(size=12) ,
    panel.background = element_rect(fill = "#F0F0F0", colour = "#F0F0F0",linewidth = 0, linetype = "solid"),
    panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid',colour = "#EFEFEF"),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(margin = margin(t = 0, r = 12, b = 0, l = 0)) ,
    axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)) ,
    axis.text.x = element_text(size=11, margin = margin(t = 8, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(size=11, margin = margin(t = 0, r = 8, b = 0, l = 0)))

  model.i <- model
  modelData.i <- modelData

  if( "model" %in% names(model.i) ) { model.i <- model.i$model }
  if( "modelData" %in% names(modelData.i) ) { modelData.i <- modelData.i$modelData }

  modelData.i.PA <- modelData.i[,"PA"]
  modelData.i <- modelData.i[, -which(colnames(modelData.i) %in% c("PA","Lon","Lat","bootstrap","Set","rep"))]

  if( variablePlot %in% colnames(modelData.i) == FALSE) { stop("Error :: variablePlot not in modelData") }

  predictProbability <- function(fitted, newdata) {
    cls <- class(fitted)
    if (any(cls %in% c("xgboost", "xgb.Booster"))) {
      as.numeric(predict(fitted, newdata))
    } else if ("maxnet" %in% cls) {
      as.numeric(predict(fitted, newdata, type = "cloglog"))
    } else if ("gbm" %in% cls) {
      as.numeric(predict(fitted, newdata, n.trees = fitted$n.trees, type = "response"))
    } else {
      as.numeric(predict(fitted, newdata, type = "response"))
    }
  }

  ## -----------------------

  data.to.plot.i <- modelData.i[,variablePlot]
  data.to.plot <- modelData.i
  data.to.plot <- apply(data.to.plot,2,mean)
  data.to.plot <- do.call("rbind", replicate(100, data.to.plot, simplify = FALSE))
  data.to.plot[,variablePlot] <- seq(min(data.to.plot.i,na.rm=T),max(data.to.plot.i,na.rm=T),length.out=100)
  data.to.plot <- as.data.frame(data.to.plot)

  matrixEffect <- NULL

  for( m in 1:length(model.i)) {

    model.i.i <- model.i[[m]]

    matrixEffect.m <- data.frame(Effect = predictProbability(model.i.i, data.to.plot))

    if( ! is.null(matrixEffect) ) { matrixEffect <- cbind(matrixEffect,data.frame(m=matrixEffect.m)) }
    if( is.null(matrixEffect) ) { matrixEffect <- data.frame(m=matrixEffect.m) }

  }

  data.to.plot <- data.frame(Variable=data.to.plot[,variablePlot],value=apply(matrixEffect,1,mean,na.rm=T))
  data.to.plot <- data.to.plot[sort(data.to.plot$Variable, index.return =T)$ix,]

  ## -----------------------

  speciesDataUse <- modelData.i[ modelData.i.PA == 1,variablePlot]
  speciesDataUse <- data.frame(x=speciesDataUse,y=min(data.to.plot$value))
  speciesDataUse <- speciesDataUse[sort(speciesDataUse[,1], index.return =T)$ix,]

  ## -----------------------

  partialPlot <- ggplot() +
    geom_point( data=speciesDataUse,aes(x=x, y=y), size=2,shape=15, fill="gray", color="gray") +
    geom_line( data=data.to.plot,aes(x=Variable, y=value),color="Black", linewidth=0.5) +
    themePlot +
    ggplot2::xlab(gsub("\\."," ",variablePlot)) + ggplot2::ylab("Effect on response (probability)")

  ## -----------------------

  span.vals <- data.to.plot$value
  span.vals.var <- numeric(length(span.vals)-1)
  for(i in 1:(length(span.vals)-1)) { span.vals.var[i] <- (abs(span.vals[i+1]) - abs(span.vals[i])) / (max(span.vals)-min(span.vals)) }
  span.vals.var <- abs(span.vals.var)
  rateChangeThreshold <- which.max(span.vals.var)
  rateChangeThreshold <- data.to.plot[rateChangeThreshold,1]

  ## --

  if( data.to.plot[1,2] > data.to.plot[nrow(data.to.plot),2]  ) { tail <- "Max" } else { tail <- "Min" }

  if( tail == "Min") {

    usedThreshold <- min(speciesDataUse[,1])
    usedThresholdQ95 <- quantile(speciesDataUse[,1],0.05,na.rm=T)
    nearMinimum <- data.to.plot$value <= min(data.to.plot$value) + quantile(data.to.plot$value, 0.05)
    extremeThreshold <- if (any(nearMinimum, na.rm = TRUE)) max(data.to.plot$Variable[nearMinimum], na.rm = TRUE) else max(data.to.plot$Variable[data.to.plot$value <= min(data.to.plot$value)], na.rm = TRUE)

  }

  if( tail == "Max") {

    usedThreshold <- max(speciesDataUse[,1])
    usedThresholdQ95 <- quantile(speciesDataUse[,1],0.95,na.rm=T)
    belowMinimum <- data.to.plot$value <= min(data.to.plot$value) - quantile(data.to.plot$value, 0.05)
    atMinimum <- data.to.plot$value <= min(data.to.plot$value)
    extremeThreshold <- if (any(belowMinimum, na.rm = TRUE)) min(data.to.plot$Variable[belowMinimum], na.rm = TRUE) else min(data.to.plot$Variable[atMinimum], na.rm = TRUE)

  }

  if( data.to.plot[1,2] == data.to.plot[nrow(data.to.plot),2]  ) { rateChangeThreshold <- NA }

  if(length(rateChangeThreshold) == 0) { rateChangeThreshold <- NA }
  if(length(extremeThreshold) == 0) { extremeThreshold <- NA }
  if(length(usedThreshold) == 0) { usedThreshold <- NA }
  if(length(usedThresholdQ95) == 0) { usedThresholdQ95 <- NA }

  tippingPoints <- data.frame( maxRateChangeThreshold=rateChangeThreshold,extremeThreshold=extremeThreshold,usedThreshold=usedThreshold,usedThresholdQ95=usedThresholdQ95 )

  return(list(partialPlot=partialPlot,tippingPoints=tippingPoints, partialPlotData=data.to.plot))

}

