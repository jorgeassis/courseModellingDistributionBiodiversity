#' Calculate Model Performance Metrics
#'
#' Evaluates model predictions against observed presence/absence values using
#' various threshold-dependent and independent metrics. It utilizes functions
#' from the 'sdm' and 'ecospat' packages and calculates AIC and R-squared.
#' Allows selection of metrics based on a specific threshold criterion.
#'
#'
#' @export
#' @importFrom sdm evaluates
#' @importFrom stats glm binomial AIC lm na.omit var cor complete.cases kmeans dist as.matrix sd quantile aggregate family formula coefficients
#' @importFrom ecospat ecospat.boyce
#'

modelPerformance <- function(observed,predicted,index) {

  if( is.null(predicted)) {

    predicted.accuracy <- data.frame(threshold=NA, auc=NA, omission.rate=NA, sensitivity=NA, specificity=NA, prop.correct=NA, Kappa=NA, aicc=NA, deviance=NA, boyce=NA, tss=NA, prevalence=NA)

  }

  if( ! is.null(predicted)) {

    keepData <- which(!is.na(observed) & !is.na(predicted))
    observed <- observed[keepData]
    predicted <- predicted[keepData]

    if( length(unique(predicted)) == 1 & length(observed) > 1 ) { predicted[which(observed == 1)] <- predicted[which(observed == 1)] + 0.01 }

    predicted.accuracy <- sdm::evaluates( x=observed , p=predicted )
    predicted.accuracy <- data.frame(criteria=predicted.accuracy@threshold_based$criteria,
                                     threshold=predicted.accuracy@threshold_based$threshold,
                                     auc=predicted.accuracy@statistics$AUC,
                                     omission.rate=predicted.accuracy@threshold_based$ommission,
                                     sensitivity=predicted.accuracy@threshold_based$sensitivity,
                                     specificity=predicted.accuracy@threshold_based$specificity,
                                     prop.correct=predicted.accuracy@threshold_based$ommission + predicted.accuracy@threshold_based$commission - 1,
                                     Kappa=predicted.accuracy@threshold_based$Kappa,
                                     aicc=NA,
                                     deviance=NA,
                                     boyce=NA,
                                     tss=predicted.accuracy@threshold_based$sensitivity + predicted.accuracy@threshold_based$specificity - 1,
                                     prevalence=predicted.accuracy@threshold_based$prevalence)

    predicted.accuracy$aicc <- AIC(glm(observed~predicted, family = binomial))
    predicted.accuracy$deviance <- summary(lm(observed~predicted))$r.squared

    boyce.a <- predicted
    boyce.a <- boyce.a[!is.na(boyce.a)]
    boyce.p <- predicted[which(observed==1)]
    boyce.p <- boyce.p[!is.na(boyce.p)]

    if( length(boyce.p) == 1) { boyce.p <- c(boyce.p-0.01,boyce.p-0.005,boyce.p,boyce.p+0.005,boyce.p+0.01) }

    tryCatch( boyce.value <- ecospat::ecospat.boyce(boyce.a , boyce.p ,PEplot=FALSE)$Spearman.cor, error=function(e) { Error <<- TRUE })
    tryCatch( boyce.value <- ecospat::ecospat.boyce(boyce.a , boyce.p ,PEplot=FALSE)$cor, error=function(e) { Error <<- TRUE })
    if( is.null(boyce.value) ) { boyce.value <- NA }

    val <- mean(boyce.p,na.rm=T)

    while( is.na(boyce.value) ) {

      val <- val - 0.001

      tryCatch( boyce.value <- ecospat::ecospat.boyce(boyce.a , c(boyce.p, val, val) ,PEplot=FALSE)$Spearman.cor, error=function(e) { Error <<- TRUE })
      tryCatch( boyce.value <- ecospat::ecospat.boyce(boyce.a , c(boyce.p, val, val) ,PEplot=FALSE)$cor, error=function(e) { Error <<- TRUE })
      if( is.null(boyce.value) ) { boyce.value <- NA }

      if( val < 0 ) { boyce.value <- 0; break }

    }

    predicted.accuracy$boyce <- boyce.value

    if( index == "tss" | index == "auc" | index == "boyce" ) { predicted.accuracy <- predicted.accuracy[which(predicted.accuracy$criteria == "max(se+sp)"),] }
    if( index == "minimumTrain" ) { predicted.accuracy <- predicted.accuracy[which(predicted.accuracy$criteria == "P0"),] }
    if( index == "minimumTrain95" ) { predicted.accuracy <- predicted.accuracy[which(predicted.accuracy$criteria == "P5"),] }

    predicted.accuracy[predicted.accuracy$tss < 0,"tss"] <- 0
    # predicted.accuracy <- predicted.accuracy[,-1]

  }

  rownames(predicted.accuracy) <- NULL
  return(predicted.accuracy)

}
