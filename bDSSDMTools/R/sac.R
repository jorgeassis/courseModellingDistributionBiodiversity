#' Estimate Spatial Autocorrelation Range and Suggest Block Size
#'
#' Computes spatial autocorrelation ranges for a set of environmental predictors
#' using empirical variograms and theoretical spherical model fitting. The function 
#' summarizes the estimated autocorrelation range (in km) for each variable, and produces 
#' a diagnostic plot showing the range per variable alongside the selected block size.
#'
#' @param environmentalLayers A \code{SpatRaster} object (from \code{terra})
#' containing environmental predictor layers.
#' @param num_sample An integer specifying the number of random spatial points to 
#' sample for calculating the empirical variogram. Default is 10000.
#'
#' @return A list with:
#' \describe{
#'   \item{plot}{A \code{ggplot} object showing spatial autocorrelation ranges per variable.}
#'   \item{distance}{A numeric value representing the suggested block size (in km).}
#' }
#'
#' @details
#' The function:
#' \itemize{
#'   \item Generates 5 random 10x10 degree geographic subsets (crops) across the overall spatial extent.
#'   \item Iteratively computes spatial autocorrelation ranges (in km) for each predictor layer within each crop using \code{gstat}.
#'   \item Extracts the minimum estimated range across the 5 random iterations for each variable to represent a conservative spatial scale.
#'   \item Caps extreme autocorrelation ranges at 500 km.
#'   \item Defines the final recommended block size as the overall mean of these capped ranges.
#' }
#'
#' @export
#' @importFrom terra SpatRaster ext crop subset spatSample
#' @importFrom gstat variogram vgm fit.variogram
#' @importFrom sf st_as_sf
#' @importFrom stats var runif apply
#' @importFrom ggplot2 ggplot aes geom_segment geom_point geom_vline annotate labs theme_minimal
#' 
#' 
#' 

sac <- function(environmentalLayers,num_sample=2500,autocorrelationClassDistance=10,autocorrelationMaxDistance=500,autocorrelationSignif=0.05) {
  
  set.seed(42)

  environmentalLayers.r <- as.data.frame(environmentalLayers, xy=TRUE, na.rm=TRUE)
  environmentalLayers.r <- environmentalLayers.r[!duplicated(environmentalLayers.r),]
  environmentalLayers.r <- environmentalLayers.r[sample(1:nrow(environmentalLayers.r), min(c(num_sample,nrow(environmentalLayers.r)))),]

  # --------

  space <- spDists(as.matrix(environmentalLayers.r[,c("x","y")]),as.matrix(environmentalLayers.r[,c("x","y")]),longlat=TRUE)
  environmentalLayers.r <- environmentalLayers.r[,! names(environmentalLayers.r) %in% c("x","y")]
  n.class <- round(autocorrelationMaxDistance / autocorrelationClassDistance)

  range <- numeric(0)
    
  for( i in 1:ncol(environmentalLayers.r)) {

      #data <- ecodist::distance( environmentalLayers.r[,1] , method = "euclidean")
      data <- spDists(as.matrix(environmentalLayers.r[,i]), as.matrix(environmentalLayers.r[,i]),longlat=FALSE)
  
      resultsMatrix <- data.frame(classdistanceFrom=seq(0,autocorrelationMaxDistance-autocorrelationClassDistance,autocorrelationClassDistance),
                                  classdistanceTo=seq(autocorrelationClassDistance,autocorrelationMaxDistance,autocorrelationClassDistance),
                                  R=NA,
                                  pVal=NA,
                                  P=NA)

      for( i in 1:nrow(resultsMatrix)) {
        
        d1 = resultsMatrix[i,1]
        d2 = resultsMatrix[i,2]
        
        data.d <- as.vector(data)
        space.d <- as.vector(space)
        
        remove <- which(space.d < d1 | space.d > d2)
        data.d <- data.d[-remove]
        space.d <- space.d[-remove]

        modelobject <- lm(space.d~data.d)
        
        f <- summary(modelobject)$fstatistic
        p <-0
        
        tryCatch( p <- pf(f[1],f[2],f[3],lower.tail=FALSE) , error=function(e) { Error <<- TRUE })
        
        resultsMatrix[i,3] <- summary(modelobject)$adj.r.squared
        resultsMatrix[i,4] <- p
        resultsMatrix[i,5] <- cor(data.d,space.d)

      }
    
      distance <- round( resultsMatrix[ which(resultsMatrix[,4] >= autocorrelationSignif)  , 2 ][1] )
      if( is.na(distance)) { distance <- autocorrelationMaxDistance }
      range <- c(range,distance)
  
  }
    
  df <- data.frame(
    variable = colnames(environmentalLayers.r),
    range = as.numeric(range)
  )

  plot1 <- ggplot(df, aes(x = range, y = reorder(variable, range))) +
      geom_segment(aes(x = 0, xend = range, yend = variable),color = "#5aa0d6", size = 1.2) +
      geom_point(size = 3, color = "#5aa0d6") +
      geom_vline(xintercept = mean(range), linetype = "dashed", color = "red") +
      annotate("text", x = mean(range), y = 1, label = "Block size",
              angle = 90, vjust = 1.5, color = "red", size = 4) +
      
      labs(
        title = "spatial autocorrelation range",
        subtitle = "Based on sample points",
        x = "Range (km)",
        y = "Variable"
      ) +
        
      theme_minimal()

  return(list(plot=plot1,distance=mean(range)))

  
  
}
