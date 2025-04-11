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

# Load main functions
library(terra)
library(bDSSDMTools)

## -----------------------
# open data

# read occurrence records (presences) and transform to spatial object
presences <- read.table('../Data/Text delimited/species_presence_records.csv',sep=';',header=TRUE)
presences <- vect(presences, geom=c("Lon","Lat"))

# download layers for present conditions
variables <- list(var1 = c("thetao","ltmax","depthmean"),
                  var2 = c("thetao","ltmin","depthmean"),
                  var4 = c("o2","ltmin","depthmean"),
                  var5 = c("phyc","ltmin","depthmean"),
                  var6 = c("sws","ltmin","depthmean"))

environmentalLayers <- download_multiple_layers(variables=variables, experiment=c("baseline"), decade=c(2010))

# change names for simplicity
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
plot(myStudyRegion, main="Study region", col="#979797", axes=TRUE)

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

## -----------------------
# generate pseudo-absence and model data

# generate pseudo-absence and combine presences and background into a unique object
absences <- pseudoAbsences(rasterLayers=environmentalLayers,presences, n=1000)
presences$PA <- 1
absences$PA <- 0
records <- rbind(absences,presences)

# generate a landmass for representation
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")
landmassMed <- crop(landmass, myStudyRegion)
# plot records
plot(landmassMed, main="Species Records", col="#D4D4D4", border="#D4D4D4", axes=TRUE)
plot(records, y = "PA", col = c("#c29431", "#043259"),pch = 16, cex = 0.5, add=TRUE)

# extract environmental values and make an object with all information to model
modelData <- generateModelData(records, envConditions=environmentalLayers, method="blocks",kFolds=6)

# inspect the data partitioning
modelData$plotDatasets

## -----------------------
# fit a model using hyperparameter tuning and monotonic constrains

# define hyperparameters add monotonic constrains
hyperparametersList <- list(learning.rate=c(0.1,0.01,0.001) , interaction.depth=c(2,3,4), n.trees=seq(100,1000,by=100))

names(environmentalLayers)
monotonicConstrains <- c("TemperatureMax" = -1 , "TemperatureMin" = 1, "Oxygen" = 1, "Productivity" = 1, "CurrentVelocity" = 1)

# fit brt model
model <- trainModel(modelData, algorithm ="brt", hyperparameters=hyperparametersList, monotonicity=monotonicConstrains)

# inspect hyperparameter combination selected per fold
model$hyperparameters

## -----------------------
# performance of the model

# inspect model performance
model$performance

# inspect relative contribution of variables
modelVarContrib <- variableContribution(model,rasterLayers=environmentalLayers)
modelVarContrib$dataFrame
modelVarContrib$plot

# inspect response curves of variables
partialPlotTemperature <- partialPlot(model,modelData,variablePlot="TemperatureMax")
partialPlotTemperature$tippingPoints
partialPlotTemperature$partialPlot

partialPlotTemperature <- partialPlot(model,modelData,variablePlot="TemperatureMin")
partialPlotTemperature$tippingPoints
partialPlotTemperature$partialPlot

partialPlotTemperature <- partialPlot(model,modelData,variablePlot="Oxygen")
partialPlotTemperature$tippingPoints
partialPlotTemperature$partialPlot

## -----------------------
# predict distribution

# predict the distribution with the model over a raster stack
prediction <- predictModel(model=model,newData=environmentalLayers)
# map the predicted distribution
plot(prediction$rasterLayer, main="Predicted species distribution",col = rev(topo.colors(100)))

# plot the predicted distribution with the study region and records
plot(prediction$rasterLayer, main="Predicted species distribution", col = rev(topo.colors(100)))
plot(records[records$PA ==1,], col = "black",pch = 16, cex = 0.5, add=TRUE)

## -----------------------
# model transferability

# download layers for end-of-century conditions
environmentalLayersSSP119 <- download_multiple_layers(variables=variables, experiment=c("ssp119"), decade=c(2090))
environmentalLayersSSP245 <- download_multiple_layers(variables=variables, experiment=c("ssp245"), decade=c(2090))

# crop and mask the environmental layers to the extent of the study region
environmentalLayersSSP119 <- crop(environmentalLayersSSP119, myStudyRegion)
environmentalLayersSSP119 <- mask(environmentalLayersSSP119, myStudyRegion)

environmentalLayersSSP245 <- crop(environmentalLayersSSP245, myStudyRegion)
environmentalLayersSSP245 <- mask(environmentalLayersSSP245, myStudyRegion)

# redirect the names of the layers
names(environmentalLayersSSP119) <- c("TemperatureMax","TemperatureMin","Oxygen","Productivity","CurrentVelocity")
names(environmentalLayersSSP245) <- c("TemperatureMax","TemperatureMin","Oxygen","Productivity","CurrentVelocity")

# plot the masked environmental layers
plot(environmentalLayersSSP119, axes=TRUE)
plot(environmentalLayersSSP245, axes=TRUE)

# transfer the distribution model
predictionSSP119 <- predictModel(model=model,newData=environmentalLayersSSP119)
predictionSSP245 <- predictModel(model=model,newData=environmentalLayersSSP245)

# For two plots stacked vertically, we want 3 rows, 1 column.
par(mfrow = c(3, 1))
par(mar = c(4, 4, 2, 1))
plot(prediction$rasterLayer, main="Present-day distribution",col = rev(topo.colors(100)))
plot(predictionSSP119$rasterLayer, main="End-of-century distribution (SSP1-1.9)",col = rev(topo.colors(100)))
plot(predictionSSP245$rasterLayer, main="End-of-century distribution (SSP2-4.5)",col = rev(topo.colors(100)))

# reset plot
dev.off()

predictionDiffSSP119 <- predictionSSP119$rasterLayer - prediction$rasterLayer
predictionDiffSSP245 <- predictionSSP245$rasterLayer - prediction$rasterLayer

# For two plots stacked vertically, we want 3 rows, 1 column.
par(mfrow = c(3, 1))
par(mar = c(4, 4, 2, 1))
plot(prediction$rasterLayer, main="Present-day distribution",col = rev(topo.colors(100)))
plot(predictionDiffSSP119, main="End-of-century difference in suitability (SSP1-1.9)",col = rev(topo.colors(100)))
plot(predictionDiffSSP245, main="End-of-century difference in suitability (SSP2-4.5)",col = rev(topo.colors(100)))

# reset plot
dev.off()

# reclassify model
model$performance

rclmat <- data.frame(from=c(0, 0.42562500),
                     to=c(0.42562500, 1),
                     becomes=c(0,1))

predictionPresentReclass <- classify(prediction$rasterLayer, rclmat)
plot(predictionPresentReclass)

predictionSSP119Reclass <- classify(predictionSSP119$rasterLayer, rclmat)
predictionSSP245Reclass <- classify(predictionSSP245$rasterLayer, rclmat)

par(mfrow = c(3, 1))
par(mar = c(4, 4, 2, 1))
plot(predictionPresentReclass, main="Present-day distribution",col = c("#E1C177", "#043259"))
plot(predictionSSP119Reclass, main="End-of-century distribution (SSP1-1.9)",col = c("#E1C177", "#043259"))
plot(predictionSSP245Reclass, main="End-of-century distribution (SSP2-4.5)",col = c("#E1C177", "#043259"))

# reset plot
dev.off()

## -----------------------
# export final predictions

# save the raster layers to external files

writeRaster(prediction$rasterLayer, filename="myFile1.tif", overwrite=TRUE)
writeRaster(predictionSSP119$rasterLayer, filename="myFile2.tif", format="GTiff", overwrite=TRUE)
writeRaster(predictionSSP245$rasterLayer, filename="myFile3.tif", format="GTiff", overwrite=TRUE)

## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------