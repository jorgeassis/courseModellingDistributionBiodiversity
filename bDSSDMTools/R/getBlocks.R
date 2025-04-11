#' Assign Spatial Block Membership using K-Means Clustering
#'
#' @export
#' @importFrom terra subset crds SpatVector rbind nrow values `values<-` `names<-` crs as.data.frame vect xyFromCell is.factor `geomtype<-` geomtype ngeom crds<- `crds<-` `crs<-` `ext<-` ext Compare maths Arithsubset extract `subset,SpatVector,numeric-method` `subset,SpatRaster,missing-method` `extract,SpatRaster,SpatVector-method` `extract,SpatRaster,matrix-method` `is.na,SpatRaster-method` `anyNA,SpatRaster-method`
#' @importFrom stats kmeans dist
#' @importFrom methods is
#'

getBlocks <- function(records, kFolds) {

  if(class(records)[1] != "SpatVector") {
    stop("Error :: records must be a SpatVector object")
  }

  presences <- records[records$PA==1,]
  presences <- data.frame(PA=rep(1,nrow(presences)),crds(presences))
  absences <- records[records$PA==0,]
  absences <- data.frame(PA=rep(0,nrow(absences)),crds(absences))

  set.seed(123) # for reproducibility of random starts
  kmeans_result <- kmeans(presences[,2:3], centers = kFolds, nstart = 25)

  centres <- kmeans_result$centers
  clusters <- kmeans_result$cluster

  presences$membership <- clusters

  absences$membership <- NA
  for(i in 1:nrow(absences)) {
    distances <- apply(centres, 1, function(center) {
      sqrt(sum((absences[i, 2:3] - center)^2))
    })
    absences$membership[i] <- which.min(distances)
  }

  combCoords <- rbind(presences, absences)
  names(combCoords) <- c("PA","Lon","Lat","membership")

  return( list(blocks=combCoords))

}

