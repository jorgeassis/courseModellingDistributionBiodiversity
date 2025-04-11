#' Generate Training, Validation, and Test Datasets for SDM
#'
#' Prepares data for model training and evaluation by cleaning records,
#' splitting data into train/validation/test sets using either random folds or
#' spatial blocks, optionally creating multiple data replicates by re-sampling
#' absences (useful when presence data is scarce), extracting environmental
#' predictor values, and generating plots of the data splits.
#'
#' @export
#' @importFrom terra SpatRaster SpatVector extract nlyr values `values<-` nrow ncol ncell subset crs `crs<-` ext rast project vect as.data.frame classify mask rasterize buffer trim complete.cases init xyFromCell relate is.related centroids geomtype `geomtype<-` `names<-` units `units<-` wrap unwrap deepcopy `ext<-` maths Compare Logic Arith `Arith,SpatRaster,SpatRaster-method` `Arith,SpatRaster,numeric-method` `Arith,numeric,SpatRaster-method` `Compare,SpatRaster,numeric-method` `Compare,SpatRaster,SpatRaster-method` `Compare,numeric,SpatRaster-method` `Logic,SpatRaster,logical-method` `Logic,SpatRaster,SpatRaster-method` `Logic,!SpatRaster-method` `subset,SpatRaster,numeric-method` `subset,SpatRaster,missing-method` `extract,SpatRaster,SpatVector-method` `extract,SpatRaster,matrix-method` `is.na,SpatRaster-method` `anyNA,SpatRaster-method`
#' @importFrom stats complete.cases kmeans dist sd na.omit quantile glm binomial lm AIC var
#' @importFrom blockCV cv_spatial
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs theme_minimal coord_equal theme element_rect element_blank ggplotGrob is.ggplot
#' @importFrom cowplot plot_grid
#' @importFrom methods is
#' @importFrom grid is.grob zeroGrob # Needed by internal get_plot_legend
#'

generateModelData <- function(records, envConditions, method="random", kFolds=NULL, proportionTrain=0.7, proportionTest=0.15, proportionVal=0.15, paMinimum=100, paRounds = 10, paRatio=1 ) {

  if(! "PA" %in% names(records)) { stop("Error :: Records must contain a column named 'PA'") }
  if( class(envConditions)[1] != "SpatRaster") { stop("Error :: envConditions must be a SpatRaster object") }
  if( class(records)[1] != "SpatVector") { stop("Error :: records must be a SpatVector object") }
  if (method != "blocks" & method != "random") { stop("Error :: method must be 'blocks' or 'random'") }

  initRecords <- nrow(records)
  recordsData <- which(complete.cases(terra::extract(envConditions,records)))
  records <- records[recordsData,]
  difference <- initRecords - nrow(records)

  if( difference > 0) { cat("Warning :: ", difference, " records were removed because they have missing values in the environmental data\n") }
  if( nrow(records) == 0 ) { stop("Error :: No records available") }

  initRecords <- nrow(records)
  DuplicatedRecordsData <- which(duplicated(terra::extract(envConditions,records)))

  if( length(DuplicatedRecordsData) > 0 ) {
    records <- records[DuplicatedRecordsData,]
    difference <- initRecords - nrow(records)
    if( difference > 0) { cat("Warning :: ", difference, " records were removed because they have duplicated entries\n") }
  }

  if( nrow(records) == 0 ) { stop("Error :: No records available") }


  if( proportionTrain + proportionTest + proportionVal != 1 ) { stop("Error :: The sum of the proportions must be equal to 1") }

  if( paMinimum >= sum(records$PA == 1) ) {
    rep <- 1:paRounds
  }

  if( paMinimum < sum(records$PA == 1) ) {
    rep <- 1 # 1:paRounds  # Force repetition even when n > paMinimum
  }

  if(method == "blocks") {

    if( is.null(kFolds)){ kFolds <- min(c(sum(records$PA == 1),10)) }

    records.i <- records
    crs(records.i) <- crs(envConditions[[1]])

    blaockSize <- min( c( (((diff( c(ext(envConditions[[1]])[1],ext(envConditions[[1]])[2])) * 111 ) / 10 ) * 1000 ),
                          (((diff( c(ext(envConditions[[1]])[3],ext(envConditions[[1]])[4])) * 111 ) / 10 ) * 1000 ) ))

    blocks <- cv_spatial(x = records.i,
                         r = envConditions[[1]], # optionally add a raster layer
                         k = kFolds,
                         size = blaockSize,
                         hexagon = FALSE, # use square blocks
                         selection = "random",
                         progress = FALSE, # turn off progress bar for vignette
                         biomod2 = TRUE, report=FALSE, plot=FALSE)

    blocks <- vect(blocks$blocks)

    proportionTrain <- ( sum(records$PA == 1) - floor(sum(records$PA == 1) / kFolds) ) / ( sum(records$PA == 1) )
    proportionVal <- 1 - proportionTrain
    recordsFolds <- getBlocks(records, kFolds=kFolds )
    recordsFolds <- recordsFolds$blocks

    membership <- extract(blocks, vect(recordsFolds, geom=c("Lon","Lat")))$folds
    recordsFolds$membership <- membership
    names(recordsFolds) <- c("PA","Lon","Lat","membership")

  }

  if(method == "random") {

    if( is.null(kFolds)){ kFolds <- floor(1/proportionTest) }

    recordsFolds <- as.data.frame(records, geom="XY")
    recordsFolds$membership <- sample( rep(1:kFolds,floor((nrow(recordsFolds)+kFolds)/(kFolds)))[1:nrow(recordsFolds)] )
    recordsFolds <- recordsFolds[, c("PA","x","y","membership")]
    names(recordsFolds) <- c("PA","Lon","Lat","membership")

  }

  modelData <- data.frame()

  for (c in 1:kFolds) {

    numberRecords.test.P <- sum(recordsFolds$membership == c & recordsFolds$PA == 1)
    numberRecord.tests.A <- sum(recordsFolds$membership == c & recordsFolds$PA == 0)

    numberRecords.train.P <- sum(recordsFolds$membership != c & recordsFolds$PA == 1)
    numberRecord.train.A <- sum(recordsFolds$membership != c & recordsFolds$PA == 0)

    if( numberRecords.test.P == 0 | numberRecord.tests.A == 0 | numberRecords.train.P == 0 | numberRecord.train.A == 0 ) { next }

    modelDataTest.P <- which(recordsFolds$membership == c & recordsFolds$PA == 1)
    modelDataTrain.P <- which(recordsFolds$membership != c & recordsFolds$PA == 1)

    modelDataTrain.P <- sample(modelDataTrain.P, size = floor((1-proportionVal) * length(modelDataTrain.P) ), replace = FALSE)
    modelDataValidation.P <- setdiff( which(recordsFolds$membership != c & recordsFolds$PA == 1) , c(modelDataTest.P,modelDataTrain.P) )

    paNTrain <- numberRecords.train.P
    paNTest <- numberRecords.test.P
    paNVal <- length(modelDataValidation.P)

    if( length(modelDataTest.P) == 0 ) { next }
    if( length(modelDataTrain.P) == 0 ) { next }
    if( length(modelDataValidation.P) == 0 ) { next }

    if( sum(modelDataTest.P %in% modelDataTrain.P) != 0 ) { stop("Wrong data generation") }
    if( sum(modelDataTrain.P %in% modelDataValidation.P) != 0 ) { stop("Wrong data generation") }
    if( sum(modelDataTest.P %in% modelDataValidation.P) != 0 ) { stop("Wrong data generation") }

    for( paRep in rep ) {

      modelDataTest.A <- which(recordsFolds$membership == c & recordsFolds$PA == 0)
      modelDataTest.A <- sample( modelDataTest.A , size = min(c(paNTest,length(modelDataTest.A))), replace = FALSE)

      modelDataTrain.A <- which(recordsFolds$membership != c & recordsFolds$PA == 0)
      modelDataTrain.A <- sample(modelDataTrain.A, size = floor((1-proportionVal) * length(modelDataTrain.A) ), replace = FALSE)
      modelDataTrain.A <- sample( modelDataTrain.A , size = min(c(paNTrain,length(modelDataTrain.A))), replace = FALSE)

      modelDataValidation.A <- setdiff( which(recordsFolds$membership != c & recordsFolds$PA == 0) , c(modelDataTest.A,modelDataTrain.A) )
      modelDataValidation.A <- sample( modelDataValidation.A , size = min(c(paNVal,length(modelDataValidation.A))), replace = FALSE)

      if( length(modelDataTest.A) == 0 ) { next }
      if( length(modelDataTrain.A) == 0 ) { next }
      if( length(modelDataValidation.A) == 0 ) { next }

      if( sum(modelDataTest.A %in% modelDataTrain.A) != 0 ) { stop("Wrong data generation") }
      if( sum(modelDataTrain.A %in% modelDataValidation.A) != 0 ) { stop("Wrong data generation") }
      if( sum(modelDataTest.A %in% modelDataValidation.A) != 0 ) { stop("Wrong data generation") }

      modelData <- rbind(modelData,
                         data.frame( cvFold=c, rep=paRep, Set="train", recordsFolds[modelDataTrain.P, ]),
                         data.frame( cvFold=c, rep=paRep, Set="train", recordsFolds[modelDataTrain.A, ]),
                         data.frame( cvFold=c, rep=paRep, Set="validation", recordsFolds[modelDataValidation.P, ]),
                         data.frame( cvFold=c, rep=paRep, Set="validation", recordsFolds[modelDataValidation.A, ]),
                         data.frame( cvFold=c, rep=paRep, Set="test", recordsFolds[modelDataTest.P, ]),
                         data.frame( cvFold=c, rep=paRep, Set="test", recordsFolds[modelDataTest.A, ]))

      modelData <- unique(modelData)

    }

  }

  # Make plot

  plots <- list()

  for( k in 1:kFolds) {

    combCoords <- modelData[modelData$cvFold == k,]
    combCoords[combCoords$Set == "validation","Set"] <- "train"

    if(nrow(combCoords) == 0 | length(unique(combCoords$PA)) != 2 ) { next }

    num_clusters <- length(unique(combCoords$Set))
    random_colors <- c("#FF8C00","#4682B4")# rgb(runif(num_clusters), runif(num_clusters), runif(num_clusters))

    if( ! exists("legendPlot")) {
      plot <- ggplot(combCoords, aes(x = Lon, y = Lat, color = Set)) +
        geom_point(alpha = 0.8, size = 1) +
        scale_color_manual(values = random_colors ) +
        theme(
          legend.position = "bottom",      # Position legend at the bottom (or "top")
          legend.direction = "horizontal", # Arrange items horizontally
          legend.background = element_rect(fill = "transparent",  colour = NA),
          legend.key = element_rect(fill = "transparent",  colour = NA)
        )

      legendPlot <- get_plot_legend(plot)

    }

    plot <- ggplot(combCoords, aes(x = Lon, y = Lat, color = Set)) +
      geom_point(alpha = 0.8, size = 0.5) +
      scale_color_manual(values = random_colors ) +
      labs(x = "", y = "", color = "Fold") +
      theme_minimal() +
      theme(legend.position = "none") +
      coord_equal()

    plots <- c(plots, list(plot))

  }

  plot <- cowplot::plot_grid(plotlist = plots, ncol = 3)
  plot2 <- cowplot::plot_grid(plot, legendPlot, ncol = 1, rel_heights = c(1, 0.1))

  modelData <- vect(modelData, geom=c("Lon","Lat"))
  modelData <- data.frame(modelData, terra::extract(envConditions, modelData, ID=FALSE))
  modelData <- modelData[complete.cases(modelData),]
  modelData <- modelData[, -which(names(modelData) %in% c("membership"))]

  names(modelData) <- c("cvFold","rep","Set","PA",names(envConditions))

  return(list(modelData=modelData,plotDatasets=plot2))

}
