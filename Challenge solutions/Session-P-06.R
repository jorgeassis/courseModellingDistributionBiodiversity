# --------------------------------------
# --------------------------------------
#
# Modelling the distribution of marine biodiversity>under global climate change
# Jorge Assis, PhD // Nord University, Norway // biodiversityDS, Centre of Marine Sciences, Portugal
#
# --------------------------------------
# --------------------------------------

# Challenge 6.1
# The main objective of this challenge is to fit a BRT model and assess its predictive performance with random cross-validation.

# load libraries
library(terra)
library(bDSSDMTools)

# read occurrence records (presences) and transform to spatial object
presences <- read.table("../Data/Text delimited/species_presence_records.csv", sep = ";", header = TRUE)
presences <- vect(presences, geom = c("Lon", "Lat"))

# download layers for present conditions
variables <- list(var1 = c("thetao", "mean", "depthmean"), var2 = c("o2", "mean", "depthmean"),
                  var3 = c("phyc", "mean", "depthmean"))

environmentalLayers <- download_multiple_layers(variables = variables, experiment = c("baseline"),
                                                decade = c(2010))
# plot layers
plot(environmentalLayers, axes = TRUE)

# develop a study region layer
myStudyRegion <- studyRegion(rasterLayers=environmentalLayers, # SpatRaster of layers
                             presences=presences, # presence records
                             distanceThreshold=200000, # distance from records in meters
                             bathymetryLayer = "../Data/Raster data/Bathymetry.tif", # bathymetry layer
                             depthPref=c(-15,-300)) # known vertical distribution

# inspect the study region layer
plot(myStudyRegion, main="Study region", col="#979797", axes=TRUE)

# crop the environmental layers to the extent of the study region
environmentalLayers <- crop(environmentalLayers, myStudyRegion)

# mask the environmental layers with the study region
environmentalLayers <- mask(environmentalLayers, myStudyRegion)

# plot the masked environmental layers
plot(environmentalLayers, axes = TRUE)

# generate pseudo-absences and combine presences and absences into a unique object
absences <- pseudoAbsences(rasterLayers = environmentalLayers, presences = presences, n = 1000)
presences$PA <- 1
absences$PA <- 0
records <- rbind(absences, presences)

# generate a landmass for representation
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")
landmassMed <- crop(landmass, myStudyRegion)

# plot records
plot(landmassMed, main = "Species Records", col = "#D4D4D4", border = "#D4D4D4", axes = TRUE)
plot(records, y = "PA", col = c("#c29431", "#043259"), pch = 16, cex = 0.5, add = TRUE)

# extract environmental values and make an object with all information to model
modelData <- generateModelData(records, envConditions = environmentalLayers, method = "random",
                               kFolds = 10)

# inspect the data partitioning
modelData$plotDatasets

# fit a brt model to the dataset
model <- trainModel(modelData, algorithm = "brt")

# predict the distribution with the model over a raster stack
prediction <- predictModel(model = model, newData = environmentalLayers)

# map the predicted distribution
plot(prediction$rasterLayer, main = "Predicted species distribution", col = rev(topo.colors(100)))

# inspect model performance
model$performance

# --------------------------------------
# --------------------------------------

# Challenge 6.2
# The main objective of this challenge is to fit a BRT model and assess its predictive performance with blocks cross-validation.

# load libraries
library(terra)
library(bDSSDMTools)

# read occurrence records (presences) and transform to spatial object
presences <- read.table("../Data/Text delimited/species_presence_records.csv", sep = ";", header = TRUE)
presences <- vect(presences, geom = c("Lon", "Lat"))

# download layers for present conditions
variables <- list(var1 = c("thetao", "mean", "depthmean"),
                  var2 = c("o2", "mean", "depthmean"),
                  var3 = c("phyc", "mean", "depthmean"))

environmentalLayers <- download_multiple_layers(variables = variables, experiment = c("baseline"),
                                                decade = c(2010))

# plot layers
plot(environmentalLayers, axes = TRUE)

# develop a study region layer
myStudyRegion <- studyRegion(rasterLayers=environmentalLayers, # SpatRaster of layers
                             presences=presences, # presence records
                             distanceThreshold=200000, # distance from records in meters
                             bathymetryLayer = "../Data/Raster data/Bathymetry.tif", # bathymetry layer
                             depthPref=c(-15,-300)) # known vertical distribution

# inspect the study region layer
plot(myStudyRegion, main="Study region", col="#979797", axes=TRUE)

# crop the environmental layers to the extent of the study region
environmentalLayers <- crop(environmentalLayers, myStudyRegion)

# mask the environmental layers with the study region
environmentalLayers <- mask(environmentalLayers, myStudyRegion)

# plot the masked environmental layers
plot(environmentalLayers, axes = TRUE)

# generate pseudo-absences and combine presences and absences into a unique object
absences <- pseudoAbsences(rasterLayers = environmentalLayers, presences = presences, n = 1000)
presences$PA <- 1
absences$PA <- 0
records <- rbind(absences, presences)

# generate a landmass for representation
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")
landmassMed <- crop(landmass, myStudyRegion)

# plot records
plot(landmassMed, main = "Species Records", col = "#D4D4D4", border = "#D4D4D4", axes = TRUE)
plot(records, y = "PA", col = c("#c29431", "#043259"), pch = 16, cex = 0.5, add = TRUE)

# extract environmental values and make an object with all information to model
modelData <- generateModelData(records, envConditions = environmentalLayers,
                               method = "blocks",
                               kFolds = 10)

# inspect the data partitioning
modelData$plotDatasets

# fit a brt model to the dataset
model <- trainModel(modelData, algorithm = "brt")

# predict the distribution with the model over a raster stack
prediction <- predictModel(model = model, newData = environmentalLayers)

# map the predicted distribution
plot(prediction$rasterLayer, main = "Predicted species distribution", col = rev(topo.colors(100)))

# inspect model performance
model$performance

# --------------------------------------
# --------------------------------------

# Challenge 6.3
# The main objective of this challenge is to evaluate a model by inspecting variable relative contribution.

# load libraries
library(terra)
library(bDSSDMTools)

# read occurrence records (presences) and transform to spatial object
presences <- read.table("../Data/Text delimited/species_presence_records.csv", sep = ";", header = TRUE)
presences <- vect(presences, geom = c("Lon", "Lat"))

# download layers for present conditions
variables <- list(var1 = c("thetao", "mean", "depthmean"),
                  var2 = c("o2", "mean", "depthmean"),
                  var3 = c("phyc", "mean", "depthmean"))

environmentalLayers <- download_multiple_layers(variables = variables, experiment = c("baseline"),
                                                decade = c(2010))

# plot layers
plot(environmentalLayers, axes = TRUE)

# develop a study region layer
myStudyRegion <- studyRegion(rasterLayers=environmentalLayers, # SpatRaster of layers
                             presences=presences, # presence records
                             distanceThreshold=200000, # distance from records in meters
                             bathymetryLayer = "../Data/Raster data/Bathymetry.tif", # bathymetry layer
                             depthPref=c(-15,-300)) # known vertical distribution

# inspect the study region layer
plot(myStudyRegion, main="Study region", col="#979797", axes=TRUE)

# crop the environmental layers to the extent of the study region
environmentalLayers <- crop(environmentalLayers, myStudyRegion)

# mask the environmental layers with the study region
environmentalLayers <- mask(environmentalLayers, myStudyRegion)

# plot the masked environmental layers
plot(environmentalLayers, axes = TRUE)

# generate pseudo-absences and combine presences and absences into a unique object
absences <- pseudoAbsences(rasterLayers = environmentalLayers, presences = presences, n = 1000)
presences$PA <- 1
absences$PA <- 0
records <- rbind(absences, presences)

# generate a landmass for representation
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")
landmassMed <- crop(landmass, myStudyRegion)

# plot records
plot(landmassMed, main = "Species Records", col = "#D4D4D4", border = "#D4D4D4", axes = TRUE)
plot(records, y = "PA", col = c("#c29431", "#043259"), pch = 16, cex = 0.5, add = TRUE)

# extract environmental values and make an object with all information to model
modelData <- generateModelData(records, envConditions = environmentalLayers,
                               method = "blocks",
                               kFolds = 10)

# inspect the data partitioning
modelData$plotDatasets

# fit a brt model to the dataset
model <- trainModel(modelData, algorithm = "brt")

# list the available layers
modelVarContrib <- variableContribution(model, rasterLayers = environmentalLayers)
modelVarContrib$plot
modelVarContrib$dataFrame

# --------------------------------------
# --------------------------------------

# Challenge 6.4
# The main objective of this challenge is to evaluate a model by inspecting response curves.

# load libraries
library(terra)
library(bDSSDMTools)

# read occurrence records (presences) and transform to spatial object
presences <- read.table("../Data/Text delimited/species_presence_records.csv", sep = ";", header = TRUE)
presences <- vect(presences, geom = c("Lon", "Lat"))

# download layers for present conditions
variables <- list(var1 = c("thetao", "mean", "depthmean"), var2 = c("o2", "mean", "depthmean"),
                  var3 = c("phyc", "mean", "depthmean"))

environmentalLayers <- download_multiple_layers(variables = variables, experiment = c("baseline"),
                                                decade = c(2010))

# plot layers
plot(environmentalLayers, axes = TRUE)

# develop a study region layer
myStudyRegion <- studyRegion(rasterLayers=environmentalLayers, # SpatRaster of layers
                             presences=presences, # presence records
                             distanceThreshold=200000, # distance from records in meters
                             bathymetryLayer = "../Data/Raster data/Bathymetry.tif", # bathymetry layer
                             depthPref=c(-15,-300)) # known vertical distribution

# inspect the study region layer
plot(myStudyRegion, main="Study region", col="#979797", axes=TRUE)

# crop the environmental layers to the extent of the study region
environmentalLayers <- crop(environmentalLayers, myStudyRegion)

# mask the environmental layers with the study region
environmentalLayers <- mask(environmentalLayers, myStudyRegion)

# plot the masked environmental layers
plot(environmentalLayers, axes = TRUE)

# generate pseudo-absences and combine presences and absences into a unique object
absences <- pseudoAbsences(rasterLayers = environmentalLayers, presences = presences, n = 1000)
presences$PA <- 1
absences$PA <- 0
records <- rbind(absences, presences)

# generate a landmass for representation
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")
landmassMed <- crop(landmass, myStudyRegion)

# plot records
plot(landmassMed, main = "Species Records", col = "#D4D4D4", border = "#D4D4D4", axes = TRUE)
plot(records, y = "PA", col = c("#c29431", "#043259"), pch = 16, cex = 0.5, add = TRUE)

# extract environmental values and make an object with all information to model
modelData <- generateModelData(records, envConditions = environmentalLayers, method = "blocks",
                               kFolds = 10)

# inspect the data partitioning
modelData$plotDatasets

# fit a brt model to the dataset
model <- trainModel(modelData, algorithm = "brt")

# list the available layers
names(environmentalLayers)

partialPlotTemperature <- partialPlot(model, modelData, variablePlot = "OceanTemperature_depthMean_Baseline_2010_mean")
partialPlotTemperature$partialPlot

# inspect the potential physiological tipping points
partialPlotTemperature$tippingPoints

# --------------------------------------
# --------------------------------------
