
library(terra)
library(ecospat)
library(dismo)
library(nicheROVER)
library(ade4)

meanNicheOverlap <- function(data,nsamples,colors.sp) {
  
  species.par <- tapply(1:nrow(data), data$species, function(ii) niw.post(nsamples = nsamples, X = data[ii, 2:ncol(data)]))
  
  # Overlap calculation.  use nsamples = nprob = 10000 (1e4) for higher
  # accuracy.  the variable over.stat can be supplied directly to the
  # overlap.plot function
  
  over.stat <- nicheROVER::overlap(species.par, nreps = nsamples, nprob = 10000, alpha = c(0.95, 0.99))
  
  cat( paste0("\n"))
  cat( paste0("............................................................."))
  cat( paste0("\n"))
  cat( paste0("\n"))
  cat( paste0("The mean overlap metrics calculated across iteratations for both niche region sizes"))
  cat( paste0("\n"))
  
  over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
  print(round(over.mean[, ,  1], 2))
  
  cat( paste0("\n"))
  cat( paste0("............................................................."))
  cat( paste0("\n"))
  cat( paste0("\n"))
  cat( paste0("Credible interval"))
  cat( paste0("\n"))
  
  over.cred <- apply(over.stat * 100, c(1:2, 4), quantile, prob = c(0.025, 0.975), na.rm = TRUE)
  print(round(over.cred[, , , 1]),2)  # display alpha = .95 niche region
  
  # Overlap plot.Before you run this, make sure that you have chosen your
  # alpha ecoregionsLevel.
  
  clrs <- colors.sp  # colors for each species
  over.stat <- nicheROVER::overlap(species.par, nreps = nsamples, nprob = 1000, alpha = 0.95)
  overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE, xlab = "Overlap Probability (%) -- Niche Region Size: 95%")
  
}

plotNicheOverlap <- function(data,nsamples) {
  
  species.par <- tapply(1:nrow(data), data$species, function(ii) niw.post(nsamples = nsamples, X = data[ii, 2:ncol(data)]))
  
  clrs <- c("black", "red")  # colors for each species
  nsamples <- 10
  species.par <- tapply(1:nrow(data), data$species, function(ii) niw.post(nsamples = nsamples, X = data[ii, 2:ncol(data)]))
  species.par.data <- tapply(1:nrow(data), data$species, function(ii) X = data[ii, 2:ncol(data)])
  niche.plot(niche.par = species.par, ndens=1000, niche.data = species.par.data, pfrac = 0, col = clrs, xlab = expression("Niche overlap"))
  
}

nicheSize <- function(data,nsamples) {
  
  species.par <- tapply(1:nrow(data), data$species, function(ii) niw.post(nsamples = nsamples, X = data[ii, 2:ncol(data)]))
  
  clrs <- c("black", "red")  # colors for each species
  species.par <- tapply(1:nrow(data), data$species, function(ii) niw.post(nsamples = nsamples, X = data[ii, 2:ncol(data)]))
  species.par.data <- tapply(1:nrow(data), data$species, function(ii) X = data[ii, 2:ncol(data)])
  niche.plot(niche.par = species.par, ndens=1000, niche.data = species.par.data, pfrac = 0, col = clrs, xlab = expression("Niche overlap"))
  
  niche.size <- sapply(species.par, function(spec) {apply(spec$Sigma, 3, niche.size, alpha = .95)})
  # point estimate and standard error
  rbind(est = colMeans(niche.size), se = apply(niche.size, 2, sd))
  
}

prepareDataNicheOverlap <- function(environment,rasterSPA,rasterSPB) {

  occurrenceRecordsSpA <- as.data.frame(rasterSPA, xy=T)
  occurrenceRecordsSpA <- occurrenceRecordsSpA[occurrenceRecordsSpA[,names(rasterSPA)] == 1,c("x","y")]
  occurrenceRecordsSpB <- as.data.frame(rasterSPB, xy=T)
  occurrenceRecordsSpB <- occurrenceRecordsSpB[occurrenceRecordsSpB[,names(rasterSPB)] == 1,c("x","y")]
  
  environment.used.by.species.1 <- terra::extract(environment,occurrenceRecordsSpA,ID=F)
  environment.used.by.species.1 <- environment.used.by.species.1[complete.cases(environment.used.by.species.1),]
  environment.used.by.species.2 <- terra::extract(environment,occurrenceRecordsSpB,ID=F)
  environment.used.by.species.2 <- environment.used.by.species.2[complete.cases(environment.used.by.species.2),]
  
  data <- data.frame(  species = c(rep(names(rasterSPA),nrow(environment.used.by.species.1)),rep(names(rasterSPB),nrow(environment.used.by.species.2))) ,
                       rbind(environment.used.by.species.1,environment.used.by.species.2) )
  
  print(aggregate(data[2:ncol(data)], data[1], mean))
  return(data)
  
}

prepareDataNicheOverlap2 <- function(environment,rasterSPA,rasterSPB) {
    
    names(rasterSPA) <- "A"
    names(rasterSPB) <- "B"
  
    occurrenceRecordsSpA <- as.data.frame(rasterSPA, xy=T)
    occurrenceRecordsSpA <- occurrenceRecordsSpA[occurrenceRecordsSpA$A == 1,c("x","y")]
    occurrenceRecordsSpB <- as.data.frame(rasterSPB, xy=T)
    occurrenceRecordsSpB <- occurrenceRecordsSpB[occurrenceRecordsSpB$B == 1,c("x","y")]
    
    environment.used.by.species.1 <- terra::extract(environment,occurrenceRecordsSpA,ID=F)
    environment.used.by.species.1 <- environment.used.by.species.1[complete.cases(environment.used.by.species.1),]
    environment.used.by.species.2 <- terra::extract(environment,occurrenceRecordsSpB,ID=F)
    environment.used.by.species.2 <- environment.used.by.species.2[complete.cases(environment.used.by.species.2),]
    
    environment.used.by.species.1 <- environment.used.by.species.1[complete.cases(environment.used.by.species.1),]
    environment.used.by.species.2 <- environment.used.by.species.2[complete.cases(environment.used.by.species.2),]
    
    environment.used.by.species.1 <- environment.used.by.species.1[sample(1:nrow(environment.used.by.species.1), min(c(nrow(environment.used.by.species.1),1000)), replace = F), ]
    environment.used.by.species.2 <- environment.used.by.species.2[sample(1:nrow(environment.used.by.species.2), min(c(nrow(environment.used.by.species.2),1000)), replace = F), ]
    
    
    env.used <- rbind(environment.used.by.species.1,environment.used.by.species.2)
    # env.used <- rbind(rnd.points.sp.1.env,rnd.points.sp.2.env)
    pca.env.used <- dudi.pca(env.used,scannf=F,nf=2)
    ecospat.plot.contrib(contrib=pca.env.used$co, eigen=pca.env.used$eig)
    
    sp.vect <- c( rep("env.sp1",nrow(environment.used.by.species.1)) , 
                  rep("env.sp2",nrow(environment.used.by.species.2)) , 
                  rep("sp1",nrow(environment.used.by.species.1)) , 
                  rep("sp2",nrow(environment.used.by.species.2)) )
    
    
    # PCA scores for the whole study area
    scores.globclim <- pca.env.used$li
    # PCA scores for the species native distribution
    scores.sp.1 <- suprow(pca.env.used,environment.used.by.species.1)$li
    scores.sp.2 <- suprow(pca.env.used,environment.used.by.species.2)$li
    
    scores.clim.sp.1 <- suprow(pca.env.used,environment.used.by.species.1)$li
    scores.clim.sp.2 <- suprow(pca.env.used,environment.used.by.species.2)$li
    
    grid.clim.sp.1 <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.clim.sp.1,
                                            sp=scores.sp.1, R=100,
                                            th.sp=0)
    grid.clim.sp.2 <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.clim.sp.2,
                                            sp=scores.sp.2, R=100,
                                            th.sp=0)
    
    grid.clim <- list("grid.clim.sp.1" = grid.clim.sp.1, "grid.clim.sp.2" = grid.clim.sp.2)
    return(grid.clim)
  
}


