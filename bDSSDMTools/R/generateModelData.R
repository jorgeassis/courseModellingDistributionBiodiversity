#' Generate Training, Validation, and Test Datasets for SDM
#'
#' splitting data into train/validation/test sets using either random folds or
#' spatial blocks, optionally creating multiple data replicates by re-sampling
#' absences (useful when presence data is scarce), extracting environmental
#' predictor values, and generating plots of the data splits.
#'
#' @export
#' @importFrom terra SpatRaster SpatVector extract nlyr values `values<-` nrow ncol ncell subset crs `crs<-` ext rast project vect as.data.frame classify mask rasterize buffer trim complete.cases init xyFromCell relate is.related centroids geomtype `geomtype<-` `names<-` units `units<-` wrap unwrap deepcopy `ext<-` maths Compare Logic Arith `Arith,SpatRaster,SpatRaster-method` `Arith,SpatRaster,numeric-method` `Arith,numeric,SpatRaster-method` `Compare,SpatRaster,numeric-method` `Compare,SpatRaster,SpatRaster-method` `Compare,numeric,SpatRaster-method` `Logic,SpatRaster,logical-method` `Logic,SpatRaster,SpatRaster-method` `Logic,!SpatRaster-method` `subset,SpatRaster,numeric-method` `subset,SpatRaster,missing-method` `extract,SpatRaster,SpatVector-method` `extract,SpatRaster,matrix-method` `is.na,SpatRaster-method` `anyNA,SpatRaster-method`
#' @importFrom stats complete.cases kmeans dist sd na.omit quantile glm binomial lm AIC var
#' @importFrom blockCV cv_spatial
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs theme_minimal coord_equal theme element_rect element_blank ggplotGrob is_ggplot
#' @importFrom cowplot plot_grid
#' @importFrom methods is
#' @importFrom grid is.grob zeroGrob # Needed by internal get_plot_legend
#'

generateModelData <- function(records, environmentalLayers, method="random", proportionTrain=0.7, proportionTest=0.15, proportionVal=0.15, bootstrapRounds = 6, paMinimum=100, paRatio=1, sacDist=NULL ) {

  if(! "PA" %in% names(records)) { stop("Error :: Records must contain a column named 'PA'") }
  if( class(environmentalLayers)[1] != "SpatRaster") { stop("Error :: environmentalLayers must be a SpatRaster object") }
  if( class(records)[1] != "SpatVector") { stop("Error :: records must be a SpatVector object") }
  if (method != "blocks" & method != "random") { stop("Error :: method must be 'blocks' or 'random'") }

  if(is.null(sacDist)) { sacDist <- 0 }

  records <- data.frame(PA=records$PA,terra::extract(environmentalLayers, records, xy=TRUE)[,c("x","y")])
  records <- records[!duplicated(records),]
  records <- vect(records, geom=c("x","y"))

  initRecords <- nrow(records)
  recordsData <- which(complete.cases(terra::extract(environmentalLayers,records)))
  records <- records[recordsData,]
  difference <- initRecords - nrow(records)

  if( difference > 0) { cat("Warning :: ", difference, " records were removed because they have NA values in the environmental data\n") }
  if( nrow(records) == 0 ) { stop("Error :: No records available") }

  if( proportionTrain + proportionTest + proportionVal != 1 ) { stop("Error :: The sum of the proportions must be equal to 1") }

  initRecords <- nrow(records)
  rep <- 1:bootstrapRounds

  if(method == "blocks") {

    records.i <- records
    crs(records.i) <- crs(environmentalLayers[[1]])

    blocks <- blockCV::cv_spatial(x = records.i,
                                  r = environmentalLayers[[1]], # optionally add a raster layer
                                  size = sacDist * 1000,
                                  hexagon = FALSE, # use square blocks
                                  selection = "random",
                                  progress = FALSE, # turn off progress bar for vignette
                                  biomod2 = TRUE, report=FALSE, plot=FALSE)

    blocks <- vect(blocks$blocks)

    recordsFolds <- data.frame( PA= as.data.frame(records.i)[,"PA"],
                                Lon = geom(records.i)[,c("x")] , Lat = geom(records.i)[,c("y")] ,
                                membership = terra::extract(blocks, records.i)[,"block_id"][which(!duplicated(terra::extract(blocks, records.i)[,"id.y"]))])

  }

  if(method == "random") {

    recordsFolds <- as.data.frame(records, geom="XY")
    recordsFolds$membership <- 1:nrow(recordsFolds)
    recordsFolds <- recordsFolds[, c("PA","x","y","membership")]
    names(recordsFolds) <- c("PA","Lon","Lat","membership")

  }

  modelData <- data.frame()
  paRep <- 0

  while( paRep <= max(rep) ) {

    recordsFolds.p <- recordsFolds[recordsFolds$PA == 1,]
    recordsFolds.a <- recordsFolds[recordsFolds$PA == 0,]
    if( nrow(recordsFolds.p) < paMinimum ) {
      recordsFolds.a <- split_f(recordsFolds.a, target=paMinimum )
    } else { recordsFolds.a <- split_f(recordsFolds.a, target=nrow(recordsFolds.p) ) }
    recordsFolds.r <- rbind(recordsFolds.p, recordsFolds.a)
    
    recordsFolds.r <- split(recordsFolds.r, recordsFolds.r$membership)
    recordsFolds.r <- lapply(recordsFolds.r, function(group) {
      shuffled_group <- group[sample(nrow(group)), , drop = FALSE]
      return(shuffled_group)
    })
    recordsFolds.r <- recordsFolds.r[sample(length(recordsFolds.r))]
    recordsFolds.r <- do.call(rbind, recordsFolds.r)

    v <- recordsFolds.r$membership
    target_idx <- floor(length(v) * proportionTrain)
    boundary_id <- v[target_idx]
    actual_cut_idx <- max(which(v == boundary_id))
    train.id <- v[1:actual_cut_idx]

    v <- recordsFolds.r$membership[(actual_cut_idx+1):length(recordsFolds.r$membership)]
    target_idx <- floor(length(recordsFolds.r$membership) * proportionTest)
    boundary_id <- v[min(c(length(v),target_idx))]
    actual_cut_idx <- max(which(v == boundary_id))
    test.id <- v[1:actual_cut_idx]
  
    val.id <- recordsFolds.r$membership[ ! recordsFolds.r$membership %in% c(train.id,test.id) ]

    if( sum(train.id %in% test.id) > 0) { stop("Error :: Train and test sets are not mutually exclusive") }
    if( sum(train.id %in% val.id) > 0) { stop("Error :: Train and validation sets are not mutually exclusive") }
    if( sum(test.id %in% val.id) > 0) { stop("Error :: Test and validation sets are not mutually exclusive") }

    if( length(train.id) == 0 | length(test.id) == 0 | length(val.id) == 0 ) { next }

    paRep <- paRep + 1

    modelData.train <- data.frame( bootstrap=paRep, Set="train", recordsFolds.r[which(recordsFolds.r$membership %in% train.id),] )
    modelData.val <- data.frame( bootstrap=paRep, Set="validation", recordsFolds.r[which(recordsFolds.r$membership %in% val.id),] )
    modelData.test <- data.frame( bootstrap=paRep, Set="test", recordsFolds.r[which(recordsFolds.r$membership %in% test.id),] )

    if( sum(modelData.train$membership %in% modelData.test$membership) > 0) { stop("Error :: Train and test sets are not mutually exclusive") }

    modelData <- rbind(modelData, modelData.train, modelData.val, modelData.test)

  }


  summary <- data.frame( dataset = c("train","validation","test"),
                         presence = c( sum(modelData$Set == "train" & modelData$PA == 1 & modelData$bootstrap %in% 1),
                                       sum(modelData$Set == "validation" & modelData$PA == 1 & modelData$bootstrap %in% 1),
                                       sum(modelData$Set == "test" & modelData$PA == 1 & modelData$bootstrap %in% 1) ),
                         absence = c( sum(modelData$Set == "train" & modelData$PA == 0 & modelData$bootstrap %in% 1),
                                      sum(modelData$Set == "validation" & modelData$PA == 0& modelData$bootstrap %in% 1),
                                      sum(modelData$Set == "test" & modelData$PA == 0 & modelData$bootstrap %in% 1) ) )

  # Make plot
  plots <- list()

  for( k in rep) {

    combCoords <- modelData[modelData$bootstrap == k,]
    if(nrow(combCoords) == 0 | length(unique(combCoords$PA)) != 2 ) { next }

    random_colors <- c("#FF8C00","#4682B4","#54ac32ff") # rgb(runif(num_clusters), runif(num_clusters), runif(num_clusters))

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

    plot <- ggplot() +  geom_spatvector(data = blocks, aes(fill = as.factor(folds)), fill=NA, color = "black", size = 0.2) +
                        geom_point(data = combCoords,aes(x = Lon, y = Lat, color = Set),alpha = 0.8, size = 0.5) +
                        scale_color_manual(values = random_colors) +
                        labs(x = "", y = "", color = "Fold") +
                        theme_minimal() +
                        theme(legend.position = "none") +
                        coord_sf()

    plots <- c(plots, list(plot))

  }

  plot <- cowplot::plot_grid(plotlist = plots, ncol = 3)
  plot2 <- cowplot::plot_grid(plot, legendPlot, ncol = 1, rel_heights = c(1, 0.1))

  modelData <- vect(modelData, geom=c("Lon","Lat"))
  modelData <- data.frame(modelData, terra::extract(environmentalLayers, modelData, ID=FALSE))
  modelData <- modelData[complete.cases(modelData),]
  modelData <- modelData[, -which(names(modelData) %in% c("membership"))]

  names(modelData) <- c("bootstrap","Set","PA",names(environmentalLayers))

  return(list(modelData=modelData,plotDatasets=plot2,summary=summary))

  }
