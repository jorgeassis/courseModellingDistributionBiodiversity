# --------------------------------------
# --------------------------------------
#
# Modelling the distribution of marine biodiversity>under global climate change
# Jorge Assis, PhD // Nord University, Norway // biodiversityDS, Centre of Marine Sciences, Portugal
#
# --------------------------------------
# --------------------------------------

# Challenge 5.1
# The objective of this challenge is to fit a brt model to occurrence (presence / absence) and environmental data.

# load libraries
library(terra)
library(bDSSDMTools)

# load a vector defining global landmasses
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")

# read occurrence records (presences and absences)
records <- read.table("../Data/Text delimited/species_presence_absence_records.csv", sep = ";",
                      header = TRUE)

# transform data.frame into a spatial object
records <- vect(records, geom = c("Lon", "Lat"))

# plot records
plot(landmass, main = "Species distribution records", col = "#D4D4D4", border = "#D4D4D4", axes = TRUE)
plot(records, y = "PA", col = c("#c29431", "#043259"), pch = 16, cex = 0.5, add = TRUE)

# import environmental data
files <- c("../Data/Raster data/Present Benthic DissolvedMolecularOxygen Mean.tif", # Oxygen
           "../Data/Raster data/Present Benthic OceanTemperature LtMax.tif", # Temperature Max
           "../Data/Raster data/Present Benthic OceanTemperature LtMin.tif", # Temperature Min
           "../Data/Raster data/Present Benthic pH Mean.tif", # pH
           "../Data/Raster data/Present Benthic TotalPrimaryProductionPhyto Mean.tif" # Productivity
)

# produce and plot raster layers
environmentalLayers <- rast(files)
plot(environmentalLayers)

# develop a study region layer
myStudyRegion <- studyRegion(rasterLayers=environmentalLayers, # SpatRaster of layers
                             presences=records, # presence records
                             distanceThreshold=200000, # distance from records in meters
                             bathymetryLayer = "../Data/Raster data/Bathymetry.tif", # bathymetry layer
                             depthPref=c(0,-500), # known vertical distribution
                             intertidalLayer = "../Data/Raster data/CoastLine.tif" ) # consider intertidal region

# inspect the study region layer
plot(myStudyRegion, main="Study region", col="#979797", axes=TRUE)

# crop the environmental layers to the extent of the study region
environmentalLayers <- crop(environmentalLayers, myStudyRegion)

# mask the environmental layers with the study region
environmentalLayers <- mask(environmentalLayers, myStudyRegion)

# plot the masked environmental layers
plot(environmentalLayers, axes = TRUE)

# extract environmental values and make an object with all information to model
modelData <- generateModelData(records, envConditions = environmentalLayers)

# fit a binomial GLM to the dataset
model <- trainModel(modelData, algorithm = "brt")

# predict the distribution with the model over a raster stack
prediction <- predictModel(model = model, newData = environmentalLayers)

# map the predicted distribution
plot(prediction$rasterLayer, main = "Predicted species distribution", col = rev(topo.colors(100)))

# --------------------------------------
# --------------------------------------

# Challenge 5.2
# The objective of this challenge is to fit a brt model to occurrence (presence / pseudo-absence) and environmental data.

# load libraries
library(terra)
library(bDSSDMTools)

# load a vector defining global landmasses
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")

# read occurrence records (presences)
presences <- read.table("../Data/Text delimited/species_presence_records.csv", sep = ";", header = TRUE)

# transform data.frame into a spatial object
presences <- vect(presences, geom = c("Lon", "Lat"))

# plot records
plot(landmass, main = "Species distribution records", col = "#D4D4D4", border = "#D4D4D4", axes = TRUE)
plot(presences, col = c( "#043259"), pch = 16, cex = 0.5, add = TRUE)

# import environmental data
files <- c("../Data/Raster data/Present Benthic DissolvedMolecularOxygen Mean.tif", # Oxygen
           "../Data/Raster data/Present Benthic OceanTemperature LtMax.tif", # Temperature Max
           "../Data/Raster data/Present Benthic OceanTemperature LtMin.tif", # Temperature Min
           "../Data/Raster data/Present Benthic pH Mean.tif", # pH
           "../Data/Raster data/Present Benthic TotalPrimaryProductionPhyto Mean.tif" # Productivity
)

# produce and plot raster layers
environmentalLayers <- rast(files)
plot(environmentalLayers)

# develop a study region layer
myStudyRegion <- studyRegion(rasterLayers=environmentalLayers, # SpatRaster of layers
                             presences=presences, # presence records
                             distanceThreshold=200000, # distance from records in meters
                             bathymetryLayer = "../Data/Raster data/Bathymetry.tif", # bathymetry layer
                             depthPref=c(-15,-300)) # known vertical distribution# consider intertidal region

# inspect the study region layer
plot(myStudyRegion, main="Study region", col="#979797", axes=TRUE)

# crop the environmental layers to the extent of the study region
environmentalLayers <- crop(environmentalLayers, myStudyRegion)

# mask the environmental layers with the study region
environmentalLayers <- mask(environmentalLayers, myStudyRegion)

# plot the masked environmental layers
plot(environmentalLayers, axes = TRUE)

# generate pseudo-absences
absences <- pseudoAbsences(rasterLayers = environmentalLayers, presences = presences, n = 1000)

# crop landmass to the region for better representation
landmassMed <- crop(landmass, myStudyRegion)

# combine presences and absences into a unique object for modelling
presences$PA <- 1
absences$PA <- 0
records <- rbind(presences, absences)
# plot records
plot(landmassMed, main = "Species records", col = "#D4D4D4", border = "#D4D4D4", axes = TRUE)
plot(records, y = "PA", col = c("#c29431", "#043259"), pch = 16, cex = 0.5, add = TRUE)

# extract environmental values and make an object with all information to model. 
# If presences are < 100, the function will generate 10 rounds of different 1:1 pseudo-absences (this is the default)
modelData <- generateModelData(records, 
                               envConditions=environmentalLayers, 
                               paMinimum = 100, 
                               paRounds = 10, 
                               paRatio = 1)

# fit a binomial GLM to the dataset
model <- trainModel(modelData, algorithm = "brt")

# predict the distribution with the model over a raster stack
prediction <- predictModel(model = model, newData = environmentalLayers)

# map the predicted distribution
plot(prediction$rasterLayer, main = "Predicted species distribution", col = rev(topo.colors(100)))

# --------------------------------------
# --------------------------------------

# Challenge 5.3
# The objective of this challenge is to fit a maxent model to occurrence (presence / background) and environmental data.

# load libraries
library(terra)
library(bDSSDMTools)

# load a vector defining global landmasses
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")

# read occurrence records (presences)
presences <- read.table("../Data/Text delimited/species_presence_records.csv", sep = ";", header = TRUE)

# transform data.frame into a spatial object
presences <- vect(presences, geom = c("Lon", "Lat"))

# plot records
plot(landmass, main = "Species distribution records", col = "#D4D4D4", border = "#D4D4D4", axes = TRUE)
plot(presences, col = c( "#043259"), pch = 16, cex = 0.5, add = TRUE)

# import environmental data
files <- c("../Data/Raster data/Present Benthic DissolvedMolecularOxygen Mean.tif", # Oxygen
           "../Data/Raster data/Present Benthic OceanTemperature LtMax.tif", # Temperature Max
           "../Data/Raster data/Present Benthic OceanTemperature LtMin.tif", # Temperature Min
           "../Data/Raster data/Present Benthic pH Mean.tif", # pH
           "../Data/Raster data/Present Benthic TotalPrimaryProductionPhyto Mean.tif" # Productivity
)

# produce and plot raster layers
environmentalLayers <- rast(files)
plot(environmentalLayers)

# In maxent algorithm, variable names cannot have spaces.
names(environmentalLayers)
names(environmentalLayers) <- gsub(" ", "_", names(environmentalLayers))

# develop a study region layer
myStudyRegion <- studyRegion(rasterLayers=environmentalLayers, # SpatRaster of layers
                             presences=presences, # presence records
                             distanceThreshold=200000, # distance from records in meters
                             bathymetryLayer = "../Data/Raster data/Bathymetry.tif", # bathymetry layer
                             depthPref=c(-15,-300)) # known vertical distribution# consider intertidal region

# inspect the study region layer
plot(myStudyRegion, main="Study region", col="#979797", axes=TRUE)

# crop the environmental layers to the extent of the study region
environmentalLayers <- crop(environmentalLayers, myStudyRegion)

# mask the environmental layers with the study region
environmentalLayers <- mask(environmentalLayers, myStudyRegion)

# plot the masked environmental layers
plot(environmentalLayers, axes = TRUE)

# generate pseudo-absences
background <- backgroundInformation(rasterLayers = environmentalLayers, n = 10000)

# crop landmass to the region for better representation
landmassMed <- crop(landmass, myStudyRegion)

# combine presences and background into a unique object for modelling
presences$PA <- 1
background$PA <- 0
records <- rbind(background,presences)
# plot records
plot(landmassMed, main = "Species records", col = "#D4D4D4", border = "#D4D4D4", axes = TRUE)
plot(records, y = "PA", col = c("#c29431", "#043259"), pch = 16, cex = 0.5, add = TRUE)

# extract environmental values and make an object with all information to model
modelData <- generateModelData(records, envConditions = environmentalLayers)

# fit a binomial GLM to the dataset
model <- trainModel(modelData, algorithm = "maxent")

# predict the distribution with the model over a raster stack
prediction <- predictModel(model = model, newData = environmentalLayers)

# map the predicted distribution
plot(prediction$rasterLayer, main = "Predicted species distribution", col = rev(topo.colors(100)))

# --------------------------------------
# --------------------------------------

# Challenge 5.4
# The main objective of this challenge is to fit a brt model to occurrence (presence / pseudo-absence) and environmental data and transfer the model to a different period.

# load libraries
library(terra)
library(bDSSDMTools)

# load a vector defining global landmasses
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")

# read occurrence records (presences)
presences <- read.table("../Data/Text delimited/species_presence_records.csv", sep = ";", header = TRUE)

# transform data.frame into a spatial object
presences <- vect(presences, geom = c("Lon", "Lat"))

# plot records
plot(landmass, main = "Species distribution records", col = "#D4D4D4", border = "#D4D4D4", axes = TRUE)
plot(presences, col = c( "#043259"), pch = 16, cex = 0.5, add = TRUE)

# import environmental data
variables <- list(var1 = c("thetao", "mean", "depthmean"), var2 = c("o2", "mean", "depthmean"),
                  var3 = c("phyc", "mean", "depthmean"))
environmentalLayersPres <- download_multiple_layers(variables = variables, experiment = c("baseline"),
                                                    decade = c(2010), longitude = c(-15, 45), latitude = c(25, 50))
# plot layers
plot(environmentalLayersPres, axes = TRUE)

# download layers for end-of-century conditions
environmentalLayersFut <- download_multiple_layers(variables = variables, experiment = c("ssp245"),
                                                   decade = c(2090), longitude = c(-15, 45), latitude = c(25, 50))
# plot layers
plot(environmentalLayersFut, axes = TRUE)

# develop a study region layer
myStudyRegion <- studyRegion(rasterLayers=environmentalLayers, # SpatRaster of layers
                             presences=presences, # presence records
                             distanceThreshold=200000, # distance from records in meters
                             bathymetryLayer = "../Data/Raster data/Bathymetry.tif", # bathymetry layer
                             depthPref=c(-15,-300)) # known vertical distribution# consider intertidal region

# inspect the study region layer
plot(myStudyRegion, main="Study region", col="#979797", axes=TRUE)

# crop the environmental layers to the extent of the study region
environmentalLayersPres <- crop(environmentalLayersPres, myStudyRegion)

# mask the environmental layers with the study region
environmentalLayersPres <- mask(environmentalLayersPres, myStudyRegion)

# plot the masked environmental layers
plot(environmentalLayersPres, axes = TRUE)

# crop the environmental layers to the extent of the study region
environmentalLayersFut <- crop(environmentalLayersFut, myStudyRegion)

# mask the environmental layers with the study region
environmentalLayersFut <- mask(environmentalLayersFut, myStudyRegion)

# plot the masked environmental layers
plot(environmentalLayersFut, axes = TRUE)

# test if names are the same
all.equal(names(environmentalLayersPres), names(environmentalLayersFut))

# assign new names to the datasets
names(environmentalLayersPres) <- c("Temperature", "Oxygen", "Productivity")
names(environmentalLayersFut) <- c("Temperature", "Oxygen", "Productivity")

# test if names are the same
all.equal(names(environmentalLayersPres), names(environmentalLayersFut))

# generate pseudo-absences
absences <- pseudoAbsences(rasterLayers = environmentalLayersPres, presences = presences, n = 1000)
# combine presences and absences into a unique object for modelling
presences$PA <- 1
absences$PA <- 0
records <- rbind(presences, absences)

# extract environmental values and make an object with all information to model
modelData <- generateModelData(records, envConditions = environmentalLayersPres)

# fit a BRT model to the dataset
model <- trainModel(modelData, algorithm = "brt")

# predict the distribution to the present-day conditions
predictionPres <- predictModel(model = model, newData = environmentalLayersPres)

# map the predicted distribution
plot(predictionPres$rasterLayer, main = "Present-day predicted distribution", col = rev(topo.colors(100)))

# predict the distribution to the present-day conditions
predictionFut <- predictModel(model = model, newData = environmentalLayersFut)

# map the predicted distribution
plot(predictionFut$rasterLayer, main = "End-of-century predicted distribution", col = rev(topo.colors(100)))

# plot difference in predicted distribution
diffPrediction <- predictionFut$rasterLayer - predictionPres$rasterLayer
plot(diffPrediction, main = "Difference in predicted distribution (SSP2-4.5)", col = rev(topo.colors(100)))

# --------------------------------------
# --------------------------------------

# Challenge 5.5
# The main objective this challenge is to fit a MaxEnt model to occurrence (presence / background) and environmental data and transfer the model to different locations.

# load libraries
library(terra)
library(bDSSDMTools)

# load a vector defining global landmasses
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")

# read occurrence records (presences)
presences <- getExternalDataGbif(taxa = "Undaria pinnatifida", getCitation = TRUE)

# transform data.frame into a spatial object
presences <- vect(presences, geom = c("Lon", "Lat"))

# load a vector defining global landmasses
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")
# plot records
plot(landmass, main = "Species Records", col = "#D4D4D4", border = "#D4D4D4", axes = TRUE)
plot(presences, col = "#043259", pch = 16, cex = 0.5, add = TRUE)

# download layers for the native range
variables <- list(var1 = c("thetao", "ltmax", "depthsurf"), var2 = c("thetao", "ltmin", "depthsurf"),
                  var3 = c("no3", "ltmin", "depthsurf"))
environmentalLayersNative <- download_multiple_layers(variables = variables, experiment = c("baseline"),
                                                      decade = c(2010), longitude = c(115, 155), latitude = c(20, 50))
# plot layers
plot(environmentalLayersNative, axes = TRUE)

# download layers for the non-native range
environmentalLayersInv <- download_multiple_layers(variables = variables, experiment = c("baseline"),
                                                   decade = c(2010), longitude = c(-15, 40), latitude = c(25, 75))
# plot layers
plot(environmentalLayersInv, axes = TRUE)

# develop a study region layer
nativeStudyRegion <- studyRegion(rasterLayers = environmentalLayersNative, presences = presences,
                                 distanceThreshold = 5e+05, intertidalLayer = "../Data/Raster data/CoastLine.tif")

# inspect the study region layer
plot(nativeStudyRegion, main = "Native study region", col = "#979797", axes = TRUE)
plot(presences, col = "#043259", pch = 16, cex = 0.5, add = TRUE)

# develop a study region layer
nonNativeStudyRegion <- studyRegion(rasterLayers = environmentalLayersInv, presences = presences,
                                    distanceThreshold = 5e+05, intertidalLayer = "../Data/Raster data/CoastLine.tif")
# inspect the study region layer
plot(nonNativeStudyRegion, main = "Non-native study region", col = "#979797", axes = TRUE)
plot(presences, col = "#043259", pch = 16, cex = 0.5, add = TRUE)

# crop the environmental layers to the extent of the study region
environmentalLayersNative <- crop(environmentalLayersNative, nativeStudyRegion)

# mask the environmental layers with the study region
environmentalLayersNative <- mask(environmentalLayersNative, nativeStudyRegion)

# plot the masked environmental layers
plot(environmentalLayersNative, axes = TRUE)

# crop the environmental layers to the extent of the study region
environmentalLayersInv <- crop(environmentalLayersInv, nonNativeStudyRegion)

# mask the environmental layers with the study region
environmentalLayersInv <- mask(environmentalLayersInv, nonNativeStudyRegion)

# plot the masked environmental layers
plot(environmentalLayersInv, axes = TRUE)

# test if names are the same
all.equal(names(environmentalLayersNative), names(environmentalLayersInv))

# generate pseudo-absences
absences <- pseudoAbsences(rasterLayers = environmentalLayersNative, presences = presences, n = 1000)

# combine presences and absences into a unique object for modelling
presences$PA <- 1
absences$PA <- 0
records <- rbind(presences, absences)

# extract environmental values and make an object with all information to model
modelData <- generateModelData(records, envConditions = environmentalLayersNative)

# fit a brt model to the dataset
model <- trainModel(modelData, algorithm = "brt")

# predict the distribution to the present-day conditions
predictionNative <- predictModel(model = model, newData = environmentalLayersNative)

# map the predicted distribution
plot(predictionNative$rasterLayer, main = "Native predicted distribution", col = rev(topo.colors(100)))

# predict the distribution to the present-day conditions
predictionInv <- predictModel(model = model, newData = environmentalLayersInv)

# map the predicted distribution
plot(predictionInv$rasterLayer, main = "Non-native predicted species distribution", col = rev(topo.colors(100)))

# map the predicted distribution
plot(predictionInv$rasterLayer, main = "Non-native predicted species distribution", col = rev(topo.colors(100)))
plot(presences, col = "#043259", pch = 16, cex = 0.5, add = TRUE)

