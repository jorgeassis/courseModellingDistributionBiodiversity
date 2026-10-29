#' Generate Training, Validation, and Test Datasets for SDM
#'
#' Presences are split into training, validation, and test sets, either at
#' random or by spatial blocks. Absences are then drawn inside each set.
#' Validation and test keep `paRatio` absences per presence. Training keeps
#' that ratio as well, and at least `paMinimum` absences, so a large MaxEnt
#' background does not unbalance the performance sets. Each bootstrap round
#' draws a new absence sample. Environmental values are extracted for the
#' records that remain, and the splits are plotted.
#'
#' @param records A terra `SpatVector` of points with a numeric `PA` column.
#'   Presences are `1` and absences are `0`.
#' @param environmentalLayers A terra `SpatRaster` of predictor layers.
#' @param method `"random"` gives each record its own group. `"blocks"` uses
#'   `blockCV` squares whose side is `sacDist` kilometres.
#' @param proportionTrain Share of presences placed in the training set.
#'   Absences are drawn afterwards from within each set.
#' @param proportionTest Share of presences placed in the test set.
#' @param proportionVal Share of presences placed in the validation set. The three
#'   proportions must sum to 1.
#' @param bootstrapRounds Number of absence resamples. Each round becomes one
#'   `bootstrap` in the returned table.
#' @param paMinimum Minimum number of absences in the training set. Training
#'   keeps at least this many, and at least `paRatio` absences for each
#'   training presence. Use a large value, such as `10000`, for MaxEnt.
#'   Validation and test do not use this minimum.
#' @param paRatio Absences kept for each presence in the validation and test
#'   sets. Training uses the same ratio unless `paMinimum` asks for more. If a
#'   set does not contain enough absences, every absence in that set is kept.
#' @param sacDist Block side in kilometres, used when `method` is `"blocks"`.
#' @export
#' @importFrom stats complete.cases reshape
#' @importFrom blockCV cv_spatial
#' @importFrom tidyterra geom_spatvector
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual labs theme_minimal theme element_rect element_blank coord_sf
#' @importFrom cowplot plot_grid
#'

generateModelData <- function(records, environmentalLayers, method="random", proportionTrain=0.7, proportionTest=0.15, proportionVal=0.15, bootstrapRounds = 6, paMinimum=100, paRatio=1, sacDist=NULL ) {

  split_f <- function(data, target, setName) {
    if (nrow(data) == 0) {
      warning("No absences are available for the ", setName, " set; ", target, " were requested.")
      return(data)
    }
    split_data <- split(data, data$membership)
    ranked_list <- lapply(split_data, function(group) {
      shuffled_group <- group[sample(nrow(group)), , drop = FALSE]
      shuffled_group$draw_rank <- seq_len(nrow(shuffled_group))
      return(shuffled_group)
    })

    ranked_df <- do.call(rbind, ranked_list)
    rownames(ranked_df) <- NULL

    ranked_df <- ranked_df[sample(nrow(ranked_df)), , drop = FALSE]
    sorted_df <- ranked_df[order(ranked_df$draw_rank), , drop = FALSE]
    if (target > nrow(sorted_df)) {
      warning("Only ", nrow(sorted_df), " absences are available for the ", setName, " set; ", target, " were requested. All absences were kept.")
      target <- nrow(sorted_df)
    }
    final_sample <- sorted_df[seq_len(target), , drop = FALSE]
    final_sample <- final_sample[, -which(names(final_sample) %in% c("draw_rank")), drop = FALSE]
    return(final_sample)
  }

  cutMemberships <- function(points) {
    if (nrow(points) < 3) return(NULL)
    grouped <- split(points, points$membership)
    grouped <- lapply(grouped, function(group) group[sample(nrow(group)), , drop = FALSE])
    grouped <- grouped[sample(length(grouped))]
    ordered <- do.call(rbind, grouped)
    v <- ordered$membership

    target_idx <- floor(length(v) * proportionTrain)
    if (target_idx < 1) return(NULL)
    boundary_id <- v[target_idx]
    actual_cut_idx <- max(which(v == boundary_id))
    if (!is.finite(actual_cut_idx) || actual_cut_idx < 1 || actual_cut_idx >= length(v)) return(NULL)
    train.id <- unique(v[seq_len(actual_cut_idx)])

    rest <- v[(actual_cut_idx + 1):length(v)]
    target_idx <- floor(length(v) * proportionTest)
    if (length(rest) == 0 || target_idx < 1) return(NULL)
    boundary_id <- rest[min(length(rest), target_idx)]
    actual_cut_idx <- max(which(rest == boundary_id))
    if (!is.finite(actual_cut_idx) || actual_cut_idx < 1 || actual_cut_idx >= length(rest)) return(NULL)
    test.id <- unique(rest[seq_len(actual_cut_idx)])
    val.id <- unique(rest[!rest %in% test.id])
    if (length(train.id) == 0 || length(test.id) == 0 || length(val.id) == 0) return(NULL)
    if (length(intersect(train.id, test.id)) > 0 || length(intersect(train.id, val.id)) > 0 || length(intersect(test.id, val.id)) > 0) {
      return(NULL)
    }
    list(train = train.id, test = test.id, validation = val.id)
  }

  if(! "PA" %in% names(records)) { stop("Error :: Records must contain a column named 'PA'") }
  if( class(environmentalLayers)[1] != "SpatRaster") { stop("Error :: environmentalLayers must be a SpatRaster object") }
  layerNames <- names(environmentalLayers)
  spacedNames <- layerNames[grepl(" ", layerNames)]
  if (length(spacedNames) > 0) {
    stop("Error :: raster layer names cannot have spaces: ", paste(spacedNames, collapse = ", "))
  }
  if( class(records)[1] != "SpatVector") { stop("Error :: records must be a SpatVector object") }
  if (method != "blocks" & method != "random") { stop("Error :: method must be 'blocks' or 'random'") }

  if(is.null(sacDist)) { sacDist <- 0 }

  records <- data.frame(PA=records$PA,terra::extract(environmentalLayers, records, xy=TRUE)[,c("x","y")])
  records <- records[!duplicated(records),]
  records <- terra::vect(records, geom=c("x","y"), crs=terra::crs(environmentalLayers[[1]]))

  initRecords <- nrow(records)
  recordsData <- which(complete.cases(terra::extract(environmentalLayers,records)))
  records <- records[recordsData,]
  difference <- initRecords - nrow(records)

  if( difference > 0) { cat("Warning :: ", difference, " records were removed because they have NA values in the environmental data\n") }
  if( nrow(records) == 0 ) { stop("Error :: No records available") }

  if( proportionTrain + proportionTest + proportionVal != 1 ) { stop("Error :: The sum of the proportions must be equal to 1") }
  if( !is.numeric(paRatio) || length(paRatio) != 1 || !is.finite(paRatio) || paRatio <= 0 ) {
    stop("Error :: paRatio must be a positive number")
  }
  if( !is.numeric(paMinimum) || length(paMinimum) != 1 || !is.finite(paMinimum) || paMinimum < 1 ) {
    stop("Error :: paMinimum must be at least 1")
  }

  initRecords <- nrow(records)
  rep <- 1:bootstrapRounds

  if(method == "blocks") {

    records.i <- records
    terra::crs(records.i) <- terra::crs(environmentalLayers[[1]])

    blocks <- blockCV::cv_spatial(x = records.i,
                                  r = environmentalLayers[[1]], # optionally add a raster layer
                                  size = sacDist * 1000,
                                  hexagon = FALSE, # use square blocks
                                  selection = "random",
                                  progress = FALSE, # turn off progress bar for vignette
                                  biomod2 = TRUE, report=FALSE, plot=FALSE)

    blocks <- terra::vect(blocks$blocks)

    extracted <- terra::extract(blocks, records.i)
    idCol <- if ("id.y" %in% names(extracted)) "id.y" else if ("ID" %in% names(extracted)) "ID" else NA_character_
    if (is.na(idCol) || !"block_id" %in% names(extracted)) {
      stop("Error :: spatial blocks did not return a block id for each record")
    }
    coords <- terra::crds(records.i)
    membership <- rep(NA_integer_, nrow(records.i))
    if (nrow(extracted) > 0) {
      first <- !duplicated(extracted[[idCol]])
      idx <- extracted[[idCol]][first]
      ok <- !is.na(idx) & idx >= 1 & idx <= length(membership)
      membership[idx[ok]] <- extracted$block_id[first][ok]
    }
    outside <- sum(is.na(membership))
    if (outside > 0) {
      cat("Warning :: ", outside, " records were removed because they fall outside spatial blocks\n")
    }
    keep <- which(!is.na(membership))
    if (length(keep) == 0) stop("Error :: no records fall inside spatial blocks")

    recordsFolds <- data.frame(
      PA = as.data.frame(records.i)[keep, "PA"],
      Lon = coords[keep, 1],
      Lat = coords[keep, 2],
      membership = membership[keep]
    )

  }

  if(method == "random") {

    recordsFolds <- as.data.frame(records, geom="XY")
    recordsFolds$membership <- 1:nrow(recordsFolds)
    recordsFolds <- recordsFolds[, c("PA","x","y","membership")]
    names(recordsFolds) <- c("PA","Lon","Lat","membership")

  }

  modelData <- data.frame()
  paRep <- 1
  splitAttempts <- 0

  while( paRep <= max(rep) ) {

    splitAttempts <- splitAttempts + 1
    if (splitAttempts > bootstrapRounds * 25) {
      stop("Error :: could not build mutually exclusive train, validation, and test sets")
    }

    recordsFolds.p <- recordsFolds[recordsFolds$PA == 1, , drop = FALSE]
    recordsFolds.a <- recordsFolds[recordsFolds$PA == 0, , drop = FALSE]
    if (nrow(recordsFolds.p) == 0) stop("Error :: no presence records")
    if (nrow(recordsFolds.a) == 0) stop("Error :: no absence records")

    sets <- cutMemberships(recordsFolds.p)
    if (is.null(sets)) next

    if (method == "blocks") {
      absenceOnly <- recordsFolds.a[!recordsFolds.a$membership %in% unique(recordsFolds.p$membership), , drop = FALSE]
      if (nrow(absenceOnly) > 0) {
        extra <- cutMemberships(absenceOnly)
        if (is.null(extra)) {
          sets$train <- c(sets$train, unique(absenceOnly$membership))
        } else {
          sets$train <- c(sets$train, extra$train)
          sets$test <- c(sets$test, extra$test)
          sets$validation <- c(sets$validation, extra$validation)
        }
      }
    }

    presTrain <- recordsFolds.p[recordsFolds.p$membership %in% sets$train, , drop = FALSE]
    presVal <- recordsFolds.p[recordsFolds.p$membership %in% sets$validation, , drop = FALSE]
    presTest <- recordsFolds.p[recordsFolds.p$membership %in% sets$test, , drop = FALSE]
    if (nrow(presTrain) == 0 || nrow(presVal) == 0 || nrow(presTest) == 0) next

    if (method == "random") {
      pool <- recordsFolds.a
      absTest <- split_f(pool, max(1, round(nrow(presTest) * paRatio)), "test")
      pool <- pool[!pool$membership %in% absTest$membership, , drop = FALSE]
      absVal <- split_f(pool, max(1, round(nrow(presVal) * paRatio)), "validation")
      pool <- pool[!pool$membership %in% absVal$membership, , drop = FALSE]
      absTrain <- split_f(pool, max(1, round(nrow(presTrain) * paRatio), paMinimum), "train")
    } else {
      absencesIn <- function(ids) recordsFolds.a[recordsFolds.a$membership %in% ids, , drop = FALSE]
      absTest <- split_f(absencesIn(sets$test), max(1, round(nrow(presTest) * paRatio)), "test")
      absVal <- split_f(absencesIn(sets$validation), max(1, round(nrow(presVal) * paRatio)), "validation")
      absTrain <- split_f(absencesIn(sets$train), max(1, round(nrow(presTrain) * paRatio), paMinimum), "train")
    }
    if (nrow(absTrain) == 0 || nrow(absVal) == 0 || nrow(absTest) == 0) next

    modelData.train <- data.frame(bootstrap = paRep, Set = "train", rbind(presTrain, absTrain))
    modelData.val <- data.frame(bootstrap = paRep, Set = "validation", rbind(presVal, absVal))
    modelData.test <- data.frame(bootstrap = paRep, Set = "test", rbind(presTest, absTest))

    if (sum(modelData.train$membership %in% modelData.test$membership) > 0) stop("Error :: Train and test sets are not mutually exclusive")
    if (sum(modelData.train$membership %in% modelData.val$membership) > 0) stop("Error :: Train and validation sets are not mutually exclusive")
    if (sum(modelData.test$membership %in% modelData.val$membership) > 0) stop("Error :: Test and validation sets are not mutually exclusive")

    modelData <- rbind(modelData, modelData.train, modelData.val, modelData.test)

    paRep <- paRep + 1

  }

  counts_long <- as.data.frame(table(modelData$bootstrap, modelData$Set, modelData$PA))
  colnames(counts_long) <- c("bootstrap", "Set", "PA", "Count")
  counts_wide <- reshape(counts_long,
                         idvar = c("bootstrap", "Set"),
                         timevar = "PA",
                         direction = "wide")
  rownames(counts_wide) <- NULL
  summary <- counts_wide
  summary <- summary[order(summary$bootstrap),]

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

    plot <- ggplot()

    if(method == "blocks") {

      plot <- plot + geom_spatvector(data = blocks, aes(fill = as.factor(folds)), fill=NA, color = "black", size = 0.2)

    }

    plot <- plot + geom_point(data = combCoords,aes(x = Lon, y = Lat, color = Set),alpha = 0.8, size = 0.5) +
      scale_color_manual(values = random_colors) +
      labs(x = "", y = "", color = "Fold") +
      theme_minimal() +
      theme(legend.position = "none") +
      coord_sf()

    plots <- c(plots, list(plot))

  }

  plot <- cowplot::plot_grid(plotlist = plots, ncol = 3)
  plot2 <- cowplot::plot_grid(plot, legendPlot, ncol = 1, rel_heights = c(1, 0.1))

  modelData <- terra::vect(modelData, geom=c("Lon","Lat"))
  modelData <- data.frame(modelData, terra::extract(environmentalLayers, modelData, ID=FALSE))
  modelData <- modelData[complete.cases(modelData),]
  modelData <- modelData[, -which(names(modelData) %in% c("membership"))]

  names(modelData) <- c("bootstrap","Set","PA",names(environmentalLayers))

  return(list(modelData=modelData,plotDatasets=plot2,summary=summary))

}
