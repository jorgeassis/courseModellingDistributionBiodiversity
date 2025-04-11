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
#' @export
#' @importFrom stats glm binomial predict predict.glm aggregate sd family formula coefficients lm na.omit var cor complete.cases kmeans dist as.matrix quantile model.matrix t.test weighted.mean formula terms reformulate update delete.response get_all_vars
#' @importFrom gbm gbm predict.gbm # Changed from gbm::gbm because of import tag convention
#' @importFrom maxnet maxnet predict.maxnet
#' @importFrom xgboost xgb.DMatrix xgboost predict.xgb.Booster
#' @importFrom methods is
#'

trainModel <- function(modelData, algorithm ="brt", hyperparameters=NULL, monotonicity=NULL) {

  if( ! algorithm %in% c("adaboost","brt","xgboost","glm","maxent") ) { stop("Error :: Unknown algorithm") }

  if( "modelData" %in% names(modelData) ) { modelData <- modelData$modelData }

  performanceIndex <- "auc"

  if( is.null(hyperparameters) ) {

    if( algorithm == "brt" ) {
      hyperparameters.grid <- expand.grid(learning.rate=0.01 , interaction.depth=3, n.trees=1000)
    }

    if( algorithm == "adaboost" ) {
      hyperparameters.grid <- expand.grid(shrinkage=0.1 , degf=5, mstop=50)
    }

    if( algorithm == "xgboost" ) {
      hyperparameters.grid <- expand.grid(shrinkage=0.1 , gamma=1, interaction.depth=3, rounds=50)
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

    if( algorithm == "adaboost" ) {

      if( ! "shrinkage" %in% names(hyperparameters.grid) ) {
        hyperparameters.grid[,"shrinkage"] <- 0.1
      }
      if( ! "degf" %in% names(hyperparameters.grid) ) {
        hyperparameters.grid[,"degf"] <- 5
      }
      if( ! "mstop" %in% names(hyperparameters.grid) ) {
        hyperparameters.grid[,"mstop"] <- 50
      }
    }

    if( algorithm == "xgboost" ) {
      if( ! "shrinkage" %in% names(hyperparameters.grid) ) {
        hyperparameters.grid[,"shrinkage"] <- 0.1
      }
      if( ! "gamma" %in% names(hyperparameters.grid) ) {
        hyperparameters.grid[,"gamma"] <- 1
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

  varsInModelData <- colnames(modelData)[! colnames(modelData) %in% c("cvFold","rep","Set","PA","Lon","Lat")]

  if( algorithm == "maxent" & sum(grepl(" ",varsInModelData)) > 0 ) { stop("In maxent algorithm, variable names cannot have spaces.")  }

  if( ! is.null(monotonicity) ) {

    if( ! algorithm %in% c("brt","xgboost") ) { stop("Error :: The algorithm chosen does not have monotonic constrains implemented") }

    if( length(monotonicity) != length(varsInModelData) ) { stop("Error :: The length of monotonicity must be equal to the number of variables in the model") }
    monotonicity.m <- monotonicity

  }

  if( is.null(monotonicity) ) { monotonicity.m <- rep(0,length(varsInModelData)); names(monotonicity.m) <- varsInModelData }

  cvFolds <- unique(modelData$cvFold)
  paRep <- unique(modelData$rep)

  if( algorithm == "maxent" ) {
    paRep <- 1
  }

  modelsTestPerformance <- data.frame()
  modelsValidationPerformance <- data.frame()
  modelsTrainPerformance <- data.frame()
  modelsHyperparameters <- data.frame()

  allModels <- list()
  allModelsLoc <- 0

  for( i in cvFolds ) {

    innerModels <- list()

    testData <- modelData[modelData$cvFold == i & modelData$Set == "test" & modelData$rep == 1,]
    if( length(unique(testData$PA)) != 2 ) { next }

    for( j in paRep) {

      trainData <- modelData[modelData$cvFold == i & modelData$Set == "train" & modelData$rep == j,]
      validationData <- modelData[modelData$cvFold == i & modelData$Set == "validation" & modelData$rep == j,]

      trainData <- trainData[,which(colnames(trainData) %in% c("PA",varsInModelData))]
      validationData <- validationData[,which(colnames(validationData) %in% c("PA",varsInModelData))]

      if( nrow(trainData) == 0 ) { next }
      if( nrow(validationData) == 0 ) { next }

      if( length(unique(trainData$PA)) != 2 ) { next }
      if( length(unique(validationData$PA)) != 2 ) { next }

      hyperparameters.validation.m <- data.frame()
      hyperparameters.train.m <- data.frame()

      for( m in 1:nrow(hyperparameters.grid)) {

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

          tryCatch( model.m <- maxnet(p=trainData$PA,
                                      data=trainData[,which(names(trainData) != "PA")],
                                      regmult=hyperparameters.grid[m,"betamultiplier"],
                                      classes=hyperparameters.grid[m,"feature"]), error=function(e) { Error <<- TRUE })

          pred.train <- predict(model.m, newdata=trainData, type="cloglog")
          pred.validation <- predict(model.m, newdata=validationData, type="cloglog")

        }

        if( algorithm == "adaboost" ) { stop("!")}

        if(algorithm == "xgboost") {
          trainData.x <- xgb.DMatrix(data = data.matrix(trainData[,-which(names(trainData) == "PA")]), label = trainData[,which(names(trainData) == "PA")], nthread = 1)
          model.m <- xgboost(data = trainData.x, monotone_constraints=monotonicity.m[names(trainData[,which(colnames(trainData) != "PA")])],
                             max_depth = hyperparameters.grid[m,"interaction.depth"],
                             gamma=hyperparameters.grid[m,"gamma"],
                             nrounds=hyperparameters.grid[m,"rounds"],
                             verbose = 0, nthread = 1, objective="binary:logistic")

          newData.x <- trainData[model.m$feature_names]
          newData.x = xgb.DMatrix(data = data.matrix(newData.x),  nthread = 1)
          pred.train <- predict(model.m,newData.x, type="response")
          newData.x = xgb.DMatrix(data = data.matrix(validationData[,setdiff(names(validationData),c("PA"))]),  nthread = 1)
          pred.validation <- predict(model.m,newData.x, type="response")

        }

        hyperparameters.validation <- modelPerformance(validationData$PA,pred.validation,index=performanceIndex)
        hyperparameters.train <- modelPerformance(trainData$PA,pred.train,index=performanceIndex)
        hyperparameters.validation.m <- rbind(hyperparameters.validation.m, hyperparameters.validation)
        hyperparameters.train.m <- rbind(hyperparameters.train.m, hyperparameters.train)

      }

      if( class(hyperparameters.validation.m)[1] != "data.frame") { hyperparameters.validation.m <- t(hyperparameters.validation.m)}
      if( class(hyperparameters.train.m)[1] != "data.frame") { hyperparameters.train.m <- t(hyperparameters.train.m)}

      hyperparameters.validation.m <- data.frame(fold=i,rep=j,hyperparameters.validation.m)
      hyperparameters.train.m <- data.frame(fold=i,rep=j,hyperparameters.train.m)

      hyperparameters.validation.m.index <- hyperparameters.validation.m[,performanceIndex]
      m <- which( hyperparameters.validation.m.index == hyperparameters.validation.m.index[which.max(hyperparameters.validation.m.index)])[1]
      hyperparameters.m <- hyperparameters.grid[m,]
      performance.validation.m <- hyperparameters.validation.m[m,]
      performance.train.m <- hyperparameters.train.m[m,]

      modelsHyperparameters <- rbind(modelsHyperparameters,data.frame(fold=i,rep=j,hyperparameters.m))
      modelsValidationPerformance <- rbind(modelsValidationPerformance,performance.validation.m)
      modelsTrainPerformance <- rbind(modelsTrainPerformance,performance.train.m)

      # Final model

      if( algorithm == "glm" ) {

        model.m <- glm(PA ~ ., data=rbind(trainData,validationData), family=binomial)

      }

      if( algorithm == "adaboost" ) { stop("!") }

      if(algorithm == "xgboost") {

        trainData.x = rbind(trainData,validationData)
        trainData.x <- xgb.DMatrix(data = data.matrix(trainData.x[,-which(names(trainData.x) == "PA")]), label = trainData.x[,which(names(trainData.x) == "PA")], nthread = 1)
        model.m <- xgboost(data = trainData.x, monotone_constraints=monotonicity.m[names(trainData[,which(colnames(trainData) != "PA")])],
                           max_depth = hyperparameters.grid[m,"interaction.depth"],
                           gamma=hyperparameters.grid[m,"gamma"],
                           nrounds=hyperparameters.grid[m,"rounds"],
                           verbose = 0, nthread = 1, objective="binary:logistic")

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

      }

      if( algorithm == "maxent" ) {

        model.m <- maxnet(p=rbind(trainData,validationData)$PA,
                          data=rbind(trainData,validationData)[,which(names(trainData) != "PA")],
                          regmult=hyperparameters.grid[m,"betamultiplier"],
                          classes=hyperparameters.grid[m,"feature"])

      }

      allModelsLoc <- allModelsLoc + 1
      innerModels[[j]] <- model.m
      allModels[[allModelsLoc]] <- model.m

    }

    # Test final ensemble of innerModels

    if(length(innerModels) == 0) { next }

    for(e in 1:length(innerModels) ) {

      model.m <- innerModels[[e]]

      if( algorithm == "glm" ) {
        pred.test <- predict(model.m, newdata=testData, type="response")
      }

      if( algorithm == "brt" ) {
        pred.test <- predict(model.m, newdata=testData, n.trees=model.m$n.trees, type="response")
      }

      if( algorithm == "maxent" ) {
        pred.test <- predict(model.m, newdata=testData, type="cloglog")
      }

      if( algorithm == "xgboost" ) {
        newData.x <- testData[model.m$feature_names]
        newData.x = xgb.DMatrix(data = data.matrix(newData.x),  nthread = 1)
        pred.test <- predict(model.m,newData.x, type="response")
      }

      modelsTestPerformance.e <- modelPerformance(testData$PA,pred.test,index=performanceIndex)
      modelsTestPerformance <- rbind(modelsTestPerformance,data.frame(fold=i,rep=0,modelsTestPerformance.e))

    }

  }

  res <- rbind(data.frame(dataset="test",index="Average",t( apply(modelsTestPerformance[,-which(colnames(modelsTestPerformance) %in% c("fold","rep","criteria"))],2,mean) )),
               data.frame(dataset="test",index="sd",t( apply(modelsTestPerformance[,-which(colnames(modelsTestPerformance) %in% c("fold","rep","criteria"))],2,sd) )),
               data.frame(dataset="validation",index="Average",t( apply(modelsValidationPerformance[,-which(colnames(modelsValidationPerformance) %in% c("fold","rep","criteria"))],2,mean) )),
               data.frame(dataset="validation",index="sd",t( apply(modelsValidationPerformance[,-which(colnames(modelsValidationPerformance) %in% c("fold","rep","criteria"))],2,sd) )),
               data.frame(dataset="train",index="Average",t( apply(modelsTrainPerformance[,-which(colnames(modelsTrainPerformance) %in% c("fold","rep","criteria"))],2,mean) )),
               data.frame(dataset="train",index="sd",t( apply(modelsTrainPerformance[,-which(colnames(modelsTrainPerformance) %in% c("fold","rep","criteria"))],2,sd) ))
  )

  modelsHyperparameters <- aggregate(modelsHyperparameters, by=list(modelsHyperparameters$fold), FUN=mean)
  modelsHyperparameters <- modelsHyperparameters[,-which(colnames(modelsHyperparameters) %in% c("rep","fold"))]
  names(modelsHyperparameters)[1] <- "kFold"
  res <- res[,-which(colnames(res) %in% c("omission.rate","prop.correct","aicc","prevalence","Kappa"))]

  return(list(performance=res, model=allModels, hyperparameters=modelsHyperparameters))

}
