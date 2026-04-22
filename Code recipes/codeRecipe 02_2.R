## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------
##
## Pack of Recipes
## Jorge Assis, PhD // Nord University, Norway // biodiversityDS, Centre of Marine Sciences, Portugal
##
## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------

## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------

## Recipe 2: Fitting Species Distribution Models (SDMs) Using Machine Learning
## This recipe walks through the process of fitting Species Distribution Models using a machine learning algorithm with environmental predictors, presence data, and projections under climate scenarios.

# 1. Setup
# A new R project and script file are created.
# External functions are loaded via libraries.R.

# 2. Data Loading
# Species presence records are imported and converted to spatial points.
# Environmental layers (e.g., temperature, oxygen, productivity) are downloaded for present-day conditions.

#  3. Study Region Definition
# A study region is generated based on species presences, bathymetry, and depth preferences.
# Environmental layers are cropped and masked to the study region.
# Correlation among predictors is checked.

#  4. Data Preparation
# Pseudo-absences are generated and combined with presences.
# Model data is created using environmental predictors and block cross-validation (kFolds=6).

# 5. Model Training
# A Boosted Regression Tree (BRT) model is trained with hyperparameter tuning and monotonic constraints for interpretability.

# 6. Model Evaluation
# Model performance is assessed.
# Variable contributions and response curves (e.g., tipping points) are plotted.

#  7. Prediction
# The model is used to predict current species distribution across the study region.
# The prediction is mapped, and presence records are overlaid.

# 8. Climate Scenario Projections
# Future environmental layers are downloaded for SSP1-1.9 and SSP2-4.5 scenarios (end-of-century).
# The trained model is used to predict future distributions under both scenarios.
# Differences between future and current suitability are calculated and visualized.
# Suitability is reclassified into binary presence/absence using a threshold from TSS maximization.

# 9. Export
# Predicted distributions (current and future) are exported as GeoTIFF files.

## -----------------------
## -----------------------

# Set a new project
# File -> New Project...

# Create a new R file
# File -> New file...
# File -> Save As... (e.g., SDM.R)

# reset the environment
closeAllConnections()
rm(list=ls())
gc(reset=TRUE)

# Load main functions
library(terra)
library(dismo)
library(ecodist)
library(usdm)
library(bDSSDMTools)

## -----------------------
# open data

# read occurrence records (presences) and transform to spatial object
presences <- read.table('../Data/Text delimited/myRecords.csv',sep=';',header=TRUE)
presences <- vect(presences, geom=c("Lon","Lat"))

# download layers for present conditions
variables <- list(var1 = c("thetao","ltmax","depthmean"),
                  var2 = c("thetao","ltmin","depthmean"),
                  var4 = c("o2","ltmin","depthmean"),
                  var5 = c("phyc","ltmin","depthmean"),
                  var6 = c("sws","ltmin","depthmean"))

environmentalLayers <- download_multiple_layers(variables=variables, experiment="baseline", decade=2010 , longitude = c(-12,33) , latitude = c(32,46))

# change names for simplicity
names(environmentalLayers)
names(environmentalLayers) <- c("TemperatureMax","TemperatureMin","Oxygen","Productivity","CurrentVelocity")

# plot layers
plot(environmentalLayers, axes=TRUE)

## -----------------------
# define the study region

# develop a study region layer
myStudyRegion <- studyRegion(rasterLayers=environmentalLayers, # SpatRaster of layers
                             presences=presences, # presence records
                             distanceThreshold=300000, # distance from records in meters
                             bathymetryLayer = "../Data/Raster data/Bathymetry.tif", # bathymetry layer
                             depthPref=c(-15,-250)) # known vertical distribution

# inspect the study region layer
plot(myStudyRegion, col="#979797", axes=TRUE, legend=FALSE)

# crop the environmental layers to the extent of the study region
environmentalLayers <- crop(environmentalLayers, myStudyRegion)
# mask the environmental layers with the study region
environmentalLayers <- mask(environmentalLayers, myStudyRegion)
# plot the masked environmental layers
plot(environmentalLayers, axes=TRUE)

# inspect correlation between layers
correlations <- getCorrPlot(environmentalLayers)
correlations$correlationPairs
correlations$plot
vifstep(environmentalLayers, th=10)

# subset environmentalLayers if needed

environmentalLayers <- subset(environmentalLayers, c("TemperatureMax","TemperatureMin","Oxygen","Productivity","CurrentVelocity"))

# Inspect spatial autocorrelation in the layers
sacAnalysis <- sac(environmentalLayers)

sacAnalysis$plot
sacDist <- sacAnalysis$distance

## -----------------------
# generate pseudo-absence and combine presences and background into a unique object
absences <- pseudoAbsences(rasterLayers=environmentalLayers,presences, n=2500)
presences$PA <- 1
absences$PA <- 0
records <- rbind(absences,presences)

# generate a landmass for representation
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")
landmassMed <- crop(landmass, myStudyRegion)

# plot records
plot(landmassMed, main="Occurrence records", col="#D4D4D4", border="#D4D4D4", axes=TRUE)
plot(records, y = "PA", col = c("#c29431", "#043259"),pch = 16, cex = 0.5, add=TRUE)

nrow(records)

# extract environmental values and make an object with all information to model
modelData <- generateModelData(records, environmentalLayers, method="blocks", proportionTrain=0.7, proportionTest=0.15, proportionVal=0.15, bootstrapRounds = 6, paMinimum=100, paRatio=1, sacDist)

# inspect the data partitioning
modelData$plotDatasets
modelData$summary

## -----------------------
# fit a brt model using hyperparameter tuning and monotonic constrains

# define hyperparameters add monotonic constrains
hyperparametersList <- list(learning.rate=c(0.1,0.01,0.001) , 
                            interaction.depth=c(2,3,4), 
                            n.trees=seq(100,1000,by=100))

names(environmentalLayers)
monotonicConstrains <- c("TemperatureMax" = -1 , "TemperatureMin" = 1, "Oxygen" = 1, "Productivity" = 1, "CurrentVelocity" = 1)

# fit brt model
model_brt <- trainModel(modelData, algorithm ="brt", hyperparameters=hyperparametersList, monotonicity=monotonicConstrains)

# inspect hyperparameter combination selected per fold
model_brt$hyperparameters

# inspect model performance
model_brt$performance

## -----------------------
# fit a xgboost model using hyperparameter tuning and monotonic constrains

# define hyperparameters add monotonic constrains
hyperparametersList <- list(shrinkage = c(0.1, 0.01, 0.001),
                            depth = c(2, 3, 4),
                            rounds = c(10,50, 100),
                            min_split_loss = c(0.1, 0.25, 0.5, 1))

names(environmentalLayers)
monotonicConstrains <- c("TemperatureMax" = -1 , "TemperatureMin" = 1, "Oxygen" = 1, "Productivity" = 1, "CurrentVelocity" = 1)

# fit brt model
model_xgboost <- trainModel(modelData, algorithm ="xgboost", hyperparameters=hyperparametersList, monotonicity=monotonicConstrains)

# inspect hyperparameter combination selected per fold
model_xgboost$hyperparameters

# inspect model performance
model_xgboost$performance

## -----------------------
# fit a glm model

# fit brt model
model_glm <- trainModel(modelData, algorithm ="glm")

# inspect model performance
model_glm$performance

## -----------------------
# fit a maxent model using hyperparameter tuning

background <- backgroundInformation(rasterLayers = environmentalLayers, n = 1000)
presences$PA <- 1
background$PA <- 0
records <- rbind(background, presences)

modelData <- generateModelData(records, environmentalLayers, method="blocks", proportionTrain=0.7, proportionTest=0.15, proportionVal=0.15, bootstrapRounds = 6, paMinimum=1000, paRatio=1000, sacDist)

# inspect the data partitioning
modelData$plotDatasets
modelData$summary

# define hyperparameters add monotonic constrains
hyperparametersList <- list(betamultiplier = c(0.5, 1, 5, 10, 25),
                            feature = c("ht", "lp", "lq", "tc", "th", "h", "l", "q", "p", "t"))

# fit brt model
model_maxent <- trainModel(modelData, algorithm ="maxent", hyperparameters=hyperparametersList)

# inspect hyperparameter combination selected per fold
model_maxent$hyperparameters

# inspect model performance
model_maxent$performance

## -----------------------
# ensemble algorithms

prediction_brt <- predictModel(model=model_brt,newData=environmentalLayers)
prediction_brt <- prediction_brt$rasterLayer
prediction_xgboost <- predictModel(model=model_xgboost,newData=environmentalLayers)
prediction_xgboost <- prediction_xgboost$rasterLayer
prediction_maxent <- predictModel(model=model_maxent,newData=environmentalLayers)
prediction_maxent <- prediction_maxent$rasterLayer
prediction_glm <- predictModel(model=model_glm,newData=environmentalLayers)
prediction_glm <- prediction_glm$rasterLayer

probabilityRaster <- app(c(prediction_brt, prediction_xgboost, prediction_maxent, prediction_glm), fun=mean)
uncertaintyRaster <- app(c(prediction_brt, prediction_xgboost, prediction_maxent, prediction_glm), fun=sd)
  
# map the predicted distribution
plot(probabilityRaster, main="Predicted species distribution",col = rev(topo.colors(100)))
plot(uncertaintyRaster, main="Uncertainty of predicted species distribution",col = rev(topo.colors(100)))

## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------