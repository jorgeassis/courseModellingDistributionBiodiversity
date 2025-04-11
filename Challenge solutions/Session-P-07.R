# --------------------------------------
# --------------------------------------
#
# Modelling the distribution of marine biodiversity>under global climate change
# Jorge Assis, PhD // Nord University, Norway // biodiversityDS, Centre of Marine Sciences, Portugal
#
# --------------------------------------
# --------------------------------------

# Challenge 7.1
# The main objective of this challenge is to perform data-driven hyperparameter selection in a brt model.

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

# Make hyperparameter list
hyperparametersList <- list(learning.rate = c(0.1, 0.01, 0.001), 
                            interaction.depth = c(2, 3,4), 
                            n.trees = seq(100, 1000, by = 100))

# train the model with non-default hyperparameters
model <- trainModel(modelData, algorithm = "brt", hyperparameters = hyperparametersList)

# inspect the hyperparameters selected in cross-validation
model$hyperparameters

# --------------------------------------
# --------------------------------------

# Challenge 7,2
# The main objective of this challenge is to fit a brt model using monotonic constrains.

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

# monotonic constrains
monoton <- c("OceanTemperature_depthMean_Baseline_2010_mean"=1, # negative effect
             "DissolvedMolecularOxygen_depthMean_Baseline_2010_mean"=1, # positive effect
             "TotalPhytoplankton_depthMean_Baseline_2010_mean"=1) # negative effect

# train the model and inspect the partial dependence plot of temperature effect
model <- trainModel(modelData, algorithm ="brt", monotonicity=monoton)
partialPlotTemperature <- partialPlot(model,modelData,variablePlot="OceanTemperature_depthMean_Baseline_2010_mean")
partialPlotTemperature$partialPlot
