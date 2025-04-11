#' Generate Partial Dependence Plot (Response Curve)
#'
#' Calculates and plots the partial dependence (response curve) for a specified
#' predictor variable from an ensemble of models. It shows how the predicted
#' suitability (averaged across models) changes across the observed range of the
#' chosen variable, while holding other predictor variables constant at their mean
#' values. Optionally calculates potential threshold/tipping points. Includes a
#' rug plot showing the distribution of the variable for presence points.
#'
#' @export
#' @importFrom stats predict na.omit quantile coefficients lm sd glm binomial kmeans dist as.matrix aggregate family formula terms reformulate update delete.response get_all_vars weighted.mean t.test
#' @importFrom ggplot2 ggplot aes geom_point geom_line theme labs element_text element_rect element_line element_blank margin annotate ggplotGrob is.ggplot rel
#' @importFrom gbm predict.gbm
#' @importFrom maxnet predict.maxnet
#' @importFrom xgboost xgb.DMatrix predict.xgb.Booster
#' @importFrom methods is
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

  if( "model" %in% names(model) ) { model <- model$model }
  if( "modelData" %in% names(modelData) ) { modelData <- modelData$modelData }

  modelData.PA <- modelData[,"PA"]
  modelData <- modelData[, -which(colnames(modelData) %in% c("PA","Lon","Lat","cvFold","Set","rep"))]

  if( variablePlot %in% colnames(modelData) == FALSE) { stop("Error :: variablePlot not in modelData") }

  yLimits <- c(0,1)

  logit2prob <- function(logit){
    odds <- exp(logit)
    prob <- odds / (1 + odds)
    return(prob)
  }

  ## -----------------------

  data.to.plot.i <- modelData[,variablePlot]
  data.to.plot <- modelData
  data.to.plot <- apply(data.to.plot,2,mean)
  data.to.plot <- do.call("rbind", replicate(100, data.to.plot, simplify = FALSE))
  data.to.plot[,variablePlot] <- seq(min(data.to.plot.i,na.rm=T),max(data.to.plot.i,na.rm=T),length.out=100)
  data.to.plot <- as.data.frame(data.to.plot)

  matrixEffect <- NULL

  for( m in 1:length(model)) {

    model.i <- model[[m]]

    if ( "xgb.Booster" %in% class(model.i) ) {

      matrixEffect.m <- data.to.plot
      matrixEffect.m <- matrixEffect.m[model.i$feature_names]
      matrixEffect.m <- xgb.DMatrix(data = data.matrix(matrixEffect.m), label = rep(0,nrow(matrixEffect.m)))
      matrixEffect.m <- data.frame(Effect= predict(model.i,matrixEffect.m, type = "response") )

    }

    if ( ! "xgb.Booster" %in% class(model.i)  ) {

      matrixEffect.m <- data.frame(Effect= predict(model.i,data.to.plot))

    }

    matrixEffect.m <- logit2prob(matrixEffect.m)

    if( ! is.null(matrixEffect) ) { matrixEffect <- cbind(matrixEffect,data.frame(m=matrixEffect.m)) }
    if( is.null(matrixEffect) ) { matrixEffect <- data.frame(m=matrixEffect.m) }

  }

  data.to.plot <- data.frame(Variable=data.to.plot[,variablePlot],value=apply(matrixEffect,1,mean,na.rm=T))
  data.to.plot <- data.to.plot[sort(data.to.plot$Variable, index.return =T)$ix,]

  if(data.to.plot[nrow(data.to.plot),2] > data.to.plot[1,2]) {
    for( x in 2:nrow(data.to.plot) ) { if( data.to.plot[x,2] < data.to.plot[x-1,2] ) { data.to.plot[x,2] <- data.to.plot[x-1,2] } }
  }
  if(data.to.plot[nrow(data.to.plot),2] < data.to.plot[1,2]) {
    for( x in 2:nrow(data.to.plot) ) { if( data.to.plot[x,2] > data.to.plot[x-1,2] ) { data.to.plot[x,2] <- data.to.plot[x-1,2] } }
  }

  ## -----------------------

  speciesDataUse <- modelData[ modelData.PA == 1,variablePlot]
  speciesDataUse <- data.frame(x=speciesDataUse,y=min(data.to.plot$value))
  speciesDataUse <- speciesDataUse[sort(speciesDataUse[,1], index.return =T)$ix,]

  ## -----------------------

  partialPlot <- ggplot() +
    geom_point( data=speciesDataUse,aes(x=x, y=y), size=2,shape=15, fill="gray", color="gray") +
    geom_line( data=data.to.plot,aes(x=Variable, y=value),color="Black", size=0.5) +
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
    extremeThreshold <- max(data.to.plot[data.to.plot$value <= min(data.to.plot$value) + quantile(data.to.plot$value,0.05),"Variable"])
    if(extremeThreshold == Inf | extremeThreshold == -Inf) { extremeThreshold <- max(data.to.plot[data.to.plot$value <= min(data.to.plot$value),"Variable"]) }

  }

  if( tail == "Max") {

    usedThreshold <- max(speciesDataUse[,1])
    usedThresholdQ95 <- quantile(speciesDataUse[,1],0.95,na.rm=T)
    extremeThreshold <- min(data.to.plot[data.to.plot$value <= min(data.to.plot$value) - quantile(data.to.plot$value,0.05) ,"Variable"])
    if(extremeThreshold == Inf | extremeThreshold == -Inf) { extremeThreshold <- min(data.to.plot[data.to.plot$value <= min(data.to.plot$value) ,"Variable"]) }

  }

  if( data.to.plot[1,2] == data.to.plot[nrow(data.to.plot),2]  ) { rateChangeThreshold <- NA }

  if(length(rateChangeThreshold) == 0) { rateChangeThreshold <- NA }
  if(length(extremeThreshold) == 0) { extremeThreshold <- NA }
  if(length(usedThreshold) == 0) { usedThreshold <- NA }
  if(length(usedThresholdQ95) == 0) { usedThresholdQ95 <- NA }

  tippingPoints <- data.frame( maxRateChangeThreshold=rateChangeThreshold,extremeThreshold=extremeThreshold,usedThreshold=usedThreshold,usedThresholdQ95=usedThresholdQ95 )

  return(list(partialPlot=partialPlot,tippingPoints=tippingPoints, partialPlotData=data.to.plot))

}

