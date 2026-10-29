#' Train and Evaluate SDM Algorithm with Hyperparameter Tuning via CV
#'
#' Trains a specified Species Distribution Model algorithm (GLM, BRT, MaxEnt, XGBoost)
#' using data prepared with cross-validation folds and repetitions (e.g., by
#' `generateModelData`). It performs hyperparameter tuning via grid search,
#' selecting the best parameters for each fold/repetition based on AUC performance
#' on the validation set. It then retrains a final model for each fold/repetition
#' on combined training and validation data using the best parameters. Test performance
#' is evaluated using an ensemble (average) prediction from all repetition models within a fold.
#' Finally, it aggregates performance metrics across folds.
#'
#' @param modelData A data frame with `bootstrap`, `Set`, and `PA` columns, or
#'   the list returned by [generateModelData()].
#' @param algorithm One of `"glm"`, `"brt"`, `"xgboost"`, or `"maxent"`.
#' @param hyperparameters A named list of values to try. `NULL` uses one
#'   default setting for the chosen algorithm.
#' @param monotonicity Numeric constraints, one per predictor, for `"brt"` and
#'   `"xgboost"`. `NULL` leaves every predictor unconstrained.
#' @param reclassCriteria Threshold rule passed to [modelPerformance()] when
#'   the test set is scored.
#' @export
#' @importFrom stats glm binomial predict aggregate sd
#' @importFrom gbm gbm predict.gbm
#' @importFrom maxnet maxnet
#' @importFrom xgboost xgb.DMatrix xgboost
#'

trainModel <- function(modelData, algorithm ="brt", hyperparameters=NULL, monotonicity=NULL,reclassCriteria="minimumTrain95") {

  warnOption <- options(warn = -1)
  on.exit(options(warnOption), add = TRUE)

  if (algorithm == "adaboost") {
    stop("Error :: adaboost is not implemented. Use brt, xgboost, glm, or maxent.")
  }
  if( ! algorithm %in% c("brt","xgboost","glm","maxent") ) { stop("Error :: Unknown algorithm. Use brt, xgboost, glm, or maxent.") }

  if( "modelData" %in% names(modelData) ) { modelData <- modelData$modelData }

  performanceIndex <- "auc.prg"

  if( is.null(hyperparameters) ) {

    if( algorithm == "brt" ) {
      hyperparameters.grid <- expand.grid(learning.rate=0.01 , interaction.depth=3, n.trees=1000)
    }

    if( algorithm == "xgboost" ) {
      hyperparameters.grid <- expand.grid(shrinkage=0.1 , min_split_loss=1, interaction.depth=3, rounds=50)
    }

    if( algorithm == "maxent" ) {
      hyperparameters.grid <- expand.grid(betamultiplier=1 , feature="ht")
    }

    if( algorithm == "glm" ) { hyperparameters.grid <- data.frame(1) }

  }

  if( ! is.null(hyperparameters) ) {

    hyperparameters.grid <- expand.grid(hyperparameters)

    if( algorithm == "brt" ) {

      if( ! "learning.rate" %in% names(hyperparameters.grid) ) {
        hyperparameters.grid[,"learning.rate"] <- 0.01
      }
      if( ! "interaction.depth" %in% names(hyperparameters.grid) ) {
        hyperparameters.grid[,"interaction.depth"] <- 3
      }
      if( ! "n.trees" %in% names(hyperparameters.grid) ) {
        hyperparameters.grid[,"n.trees"] <- 1000
      }

    }

    if( algorithm == "xgboost" ) {
      if( ! "shrinkage" %in% names(hyperparameters.grid) ) {
        hyperparameters.grid[,"shrinkage"] <- 0.1
      }
      if( ! "min_split_loss" %in% names(hyperparameters.grid) ) {
        hyperparameters.grid[,"min_split_loss"] <- 1
      }
      if( ! "interaction.depth" %in% names(hyperparameters.grid) ) {
        hyperparameters.grid[,"interaction.depth"] <- 3
      }
      if( ! "rounds" %in% names(hyperparameters.grid) ) {
        hyperparameters.grid[,"rounds"] <- 50
      }
    }

    if( algorithm == "maxent" ) {
      if( ! "betamultiplier" %in% names(hyperparameters.grid) ) {
        hyperparameters.grid[,"betamultiplier"] <- 1
      }
      if( ! "feature" %in% names(hyperparameters.grid) ) {
        hyperparameters.grid[,"feature"] <- "ht"
      }
    }

    if( algorithm == "glm" ) { hyperparameters.grid <- data.frame("none") }

  }

  varsInModelData <- colnames(modelData)[! colnames(modelData) %in% c("bootstrap","Set","PA","Lon","Lat")]

  if( algorithm == "maxent" & sum(grepl(" ",varsInModelData)) > 0 ) { stop("In maxent algorithm, variable names cannot have spaces.")  }

  if( ! is.null(monotonicity) ) {

    if( ! algorithm %in% c("brt","xgboost") ) { stop("Error :: The algorithm chosen does not have monotonic constrains implemented") }

    if( length(monotonicity) != length(varsInModelData) ) { stop("Error :: The length of monotonicity must be equal to the number of variables in the model") }
    monotonicity.m <- monotonicity

  }

  if( is.null(monotonicity) ) { monotonicity.m <- rep(0,length(varsInModelData)); names(monotonicity.m) <- varsInModelData }

  paRep <- unique(modelData$bootstrap)

  modelsTestPerformance <- data.frame()
  modelsValidationPerformance <- data.frame()
  modelsTrainPerformance <- data.frame()
  modelsHyperparameters <- data.frame()

  allModels <- list()
  allModelsLoc <- 0

  for( i in paRep) {

    trainData <- modelData[modelData$bootstrap == i & modelData$Set == "train",]
    validationData <- modelData[modelData$bootstrap == i & modelData$Set == "validation",]
    testData <- modelData[modelData$bootstrap == i & modelData$Set == "test" ,]

    trainData <- trainData[,which(colnames(trainData) %in% c("PA",varsInModelData))]
    validationData <- validationData[,which(colnames(validationData) %in% c("PA",varsInModelData))]
    testData <- testData[,which(colnames(testData) %in% c("PA",varsInModelData))]

    if( nrow(trainData) == 0 ) { next }
    if( nrow(validationData) == 0 ) { next }
    if( nrow(testData) == 0 ) { next }

    if( length(unique(testData$PA)) != 2 ) { next }
    if( length(unique(trainData$PA)) != 2 ) { next }
    if( length(unique(validationData$PA)) != 2 ) { next }

    hyperparameters.validation.m <- data.frame()
    hyperparameters.train.m <- data.frame()
    gridRow <- integer()

    for( m in 1:nrow(hyperparameters.grid)) {

      if(exists("model.m")) { rm("model.m") }

      if( algorithm == "glm" ) {

        model.m <- glm(PA ~ ., data=trainData, family=binomial)
        pred.train <- predict(model.m, newdata=trainData, type="response")
        pred.validation <- predict(model.m, newdata=validationData, type="response")

      }

      if( algorithm == "brt" ) {

        model.m <- gbm::gbm(PA ~ .,
                            data=trainData,
                            distribution="bernoulli",
                            n.trees=hyperparameters.grid[m,"n.trees"],
                            shrinkage=hyperparameters.grid[m,"learning.rate"],
                            interaction.depth=hyperparameters.grid[m,"interaction.depth"],
                            bag.fraction=0.5,
                            var.monotone=monotonicity.m,
                            n.cores=1)

        pred.train <- predict(model.m, newdata=trainData, n.trees=model.m$n.trees, type="response")
        pred.validation <- predict(model.m, newdata=validationData, n.trees=model.m$n.trees, type="response")

      }

      if( algorithm == "maxent" ) {

        model.m <- tryCatch(
          maxnet::maxnet(p=trainData$PA,
                         data=trainData[,which(names(trainData) != "PA")],
                         regmult=hyperparameters.grid[m,"betamultiplier"],
                         classes=hyperparameters.grid[m,"feature"]),
          error = function(e) e
        )
        if (inherits(model.m, "error")) next

        pred.train <- predict(model.m, newdata=trainData, type="cloglog")
        pred.validation <- predict(model.m, newdata=validationData, type="cloglog")

      }

      if(algorithm == "xgboost") {

        trainData.x <- xgboost::xgb.DMatrix(data = data.matrix(trainData[,-which(names(trainData) == "PA")]), label = trainData[,which(names(trainData) == "PA")], nthread = 1)
        model.m <- xgboost::xgboost(x = data.matrix(trainData[,-which(names(trainData) == "PA")]),
                                    y = as.factor(trainData$PA),
                                    monotone_constraints=monotonicity.m[names(trainData[,which(colnames(trainData) != "PA")])],
                                    max_depth = hyperparameters.grid[m,"interaction.depth"],
                                    min_split_loss=hyperparameters.grid[m,"min_split_loss"],
                                    learning_rate=hyperparameters.grid[m,"shrinkage"],
                                    nrounds=hyperparameters.grid[m,"rounds"],
                                    nthread = 1, objective="binary:logistic")

        newData.x <- trainData[names(trainData)[!names(trainData) %in% "PA"]]
        newData.x <- data.matrix(newData.x)
        newData.x =  xgboost::xgb.DMatrix(data = newData.x,  nthread = 1)
        pred.train <- predict(model.m,trainData, type="response")
        newData.x =  xgboost::xgb.DMatrix(data = data.matrix(validationData[,setdiff(names(validationData),c("PA"))]),  nthread = 1)
        pred.validation <- predict(model.m,validationData, type="response")

      }

      hyperparameters.validation <- modelPerformance(validationData$PA,pred.validation,index=performanceIndex)
      hyperparameters.validation <- hyperparameters.validation[which.max(hyperparameters.validation[,performanceIndex]),,drop=FALSE]
      hyperparameters.train <- modelPerformance(trainData$PA,pred.train,index=performanceIndex)
      hyperparameters.train <- hyperparameters.train[which.max(hyperparameters.train[,performanceIndex]),,drop=FALSE]
      hyperparameters.validation.m <- rbind(hyperparameters.validation.m, hyperparameters.validation)
      hyperparameters.train.m <- rbind(hyperparameters.train.m, hyperparameters.train)
      gridRow <- c(gridRow, m)

    }

    if (nrow(hyperparameters.validation.m) == 0) next

    if( class(hyperparameters.validation.m)[1] != "data.frame") { hyperparameters.validation.m <- t(hyperparameters.validation.m)}
    if( class(hyperparameters.train.m)[1] != "data.frame") { hyperparameters.train.m <- t(hyperparameters.train.m)}

    hyperparameters.validation.m <- data.frame(fold=i,hyperparameters.validation.m)
    hyperparameters.train.m <- data.frame(fold=i,hyperparameters.train.m)

    hyperparameters.validation.m.index <- hyperparameters.validation.m[, performanceIndex]
    chosen <- which( hyperparameters.validation.m.index == hyperparameters.validation.m.index[which.max(hyperparameters.validation.m.index)])[1]
    m <- gridRow[chosen]
    hyperparameters.m <- hyperparameters.grid[m,]
    performance.validation.m <- hyperparameters.validation.m[chosen,]
    performance.train.m <- hyperparameters.train.m[chosen,]

    modelsHyperparameters <- rbind(modelsHyperparameters,data.frame(fold=i,hyperparameters.m))
    modelsValidationPerformance <- rbind(modelsValidationPerformance,performance.validation.m)
    modelsTrainPerformance <- rbind(modelsTrainPerformance,performance.train.m)

    # Final model

    if( algorithm == "glm" ) {

      model.m <- glm(PA ~ ., data=rbind(trainData,validationData), family=binomial)
      pred.test <- predict(model.m, newdata=testData, type="response")

    }

    if(algorithm == "xgboost") {

      trainData.x = rbind(trainData,validationData)
      #trainData.x <- xgb.DMatrix(data = data.matrix(trainData.x[,-which(names(trainData.x) == "PA")]), label = trainData.x[,which(names(trainData.x) == "PA")], nthread = 1)
      model.m <- xgboost::xgboost(x = data.matrix(trainData.x[,-which(names(trainData.x) == "PA")]),
                                    y = as.factor(trainData.x$PA),
        monotone_constraints=monotonicity.m[names(trainData.x[,which(colnames(trainData.x) != "PA")])],
                          max_depth = hyperparameters.grid[m,"interaction.depth"],
                          min_split_loss=hyperparameters.grid[m,"min_split_loss"],
                          learning_rate=hyperparameters.grid[m,"shrinkage"],
                          nrounds=hyperparameters.grid[m,"rounds"],
                          nthread = 1, objective="binary:logistic")
      #newData.x <- testData[model.m$feature_names]
      #newData.x = xgb.DMatrix(data = data.matrix(newData.x),  nthread = 1)
      pred.test <- predict(model.m,testData, type="response")

    }

    if( algorithm == "brt" ) {

      model.m <- gbm::gbm(PA ~ .,
                          data=rbind(trainData,validationData),
                          distribution="bernoulli",
                          n.trees=hyperparameters.grid[m,"n.trees"],
                          shrinkage=hyperparameters.grid[m,"learning.rate"],
                          interaction.depth=hyperparameters.grid[m,"interaction.depth"],
                          bag.fraction=0.5,
                          var.monotone=monotonicity.m,
                          n.cores=1)

      pred.test <- predict(model.m, newdata=testData, n.trees=model.m$n.trees, type="response")

    }

    if( algorithm == "maxent" ) {

      model.m <- maxnet::maxnet(p=rbind(trainData,validationData)$PA,
                        data=rbind(trainData,validationData)[,which(names(trainData) != "PA")],
                        regmult=hyperparameters.grid[m,"betamultiplier"],
                        classes=hyperparameters.grid[m,"feature"])
      pred.test <- predict(model.m, newdata=testData, type="cloglog")

    }

    allModelsLoc <- allModelsLoc + 1
    allModels[[allModelsLoc]] <- model.m

    modelsTestPerformance.e <- modelPerformance(testData$PA,pred.test,index=reclassCriteria)
    modelsTestPerformance.e <- modelsTestPerformance.e[which.max(modelsTestPerformance.e[,performanceIndex]),,drop=FALSE]
    modelsTestPerformance <- rbind(modelsTestPerformance,data.frame(fold=i,modelsTestPerformance.e))

  }

  if (nrow(modelsTestPerformance) == 0) {
    stop("Error :: no bootstrap replicate produced a model. Each replicate needs presences and absences in the train, validation, and test sets.")
  }

  maskNegative <- function(df) {
    numericColumns <- vapply(df, is.numeric, logical(1))
    numericColumns[names(df) %in% c("fold", "bootstrap")] <- FALSE
    df[numericColumns] <- lapply(df[numericColumns], function(values) {
      values[!is.na(values) & values < 0] <- NA
      values
    })
    df
  }
  modelsTestPerformance <- maskNegative(modelsTestPerformance)
  modelsValidationPerformance <- maskNegative(modelsValidationPerformance)
  modelsTrainPerformance <- maskNegative(modelsTrainPerformance)

  res <- rbind(data.frame(dataset="test",index="Average",t( apply(modelsTestPerformance[,-which(colnames(modelsTestPerformance) %in% c("fold","criteria"))],2,mean,na.rm=T) )),
               data.frame(dataset="test",index="sd",t( apply(modelsTestPerformance[,-which(colnames(modelsTestPerformance) %in% c("fold","criteria"))],2,sd,na.rm=T) )),
               data.frame(dataset="validation",index="Average",t( apply(modelsValidationPerformance[,-which(colnames(modelsValidationPerformance) %in% c("fold","criteria"))],2,mean,na.rm=T) )),
               data.frame(dataset="validation",index="sd",t( apply(modelsValidationPerformance[,-which(colnames(modelsValidationPerformance) %in% c("fold","criteria"))],2,sd,na.rm=T) )),
               data.frame(dataset="train",index="Average",t( apply(modelsTrainPerformance[,-which(colnames(modelsTrainPerformance) %in% c("fold","criteria"))],2,mean,na.rm=T) )),
               data.frame(dataset="train",index="sd",t( apply(modelsTrainPerformance[,-which(colnames(modelsTrainPerformance) %in% c("fold","criteria"))],2,sd,na.rm=T) ))
  )

  names(modelsHyperparameters)[1] <- "bootstrap"
  res <- res[,-which(colnames(res) %in% c("omission.rate","prop.correct","aicc","prevalence","Kappa","deviance"))]

  return(list(performance=res, model=allModels, hyperparameters=modelsHyperparameters))

}
