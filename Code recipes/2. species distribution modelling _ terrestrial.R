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
library(geodata)
library(bDSSDMTools)

## -----------------------
# open data

# read occurrence records (presences) and transform to spatial object
presences <- read.table('../Data/Text delimited/myRecords3.csv',sep=';',header=TRUE)
presences <- vect(presences, geom=c("Lon","Lat"))

# download layers for present conditions from world clim
# Data present-day (baseline)
bioPresent <- worldclim_global(var = "bio", # "tmin", "tmax", "tavg", "prec", "wind", "vapr", or "bio"
                               res = 5, # resolution: 10, 5, 2.5, or 0.5 (minutes of a degree)
                               path = "_ rasters")

names(bioPresent)
plot(bioPresent)

#select the relevant variables for the modelling (tmax, tmin, precipitation)
environmentalLayers <- bioPresent[[c(5, 6, 12)]]

# change names for simplicity
names(environmentalLayers)
names(environmentalLayers) <- c("TemperatureMax","TemperatureMin","Precipitation")

# plot layers
plot(environmentalLayers, axes=TRUE)

## -----------------------
# define the study region

# develop a study region layer
myStudyRegion <- studyRegion(rasterLayers=environmentalLayers, # SpatRaster of layers
                             presences=presences, # presence records
                             distanceThreshold=500000, # distance from records in meters
                             topographyLayer = "../Data/Raster data/Elevation.tif", # elevation layer
                             verticalPref=c(0,800)) # known elevational distribution

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

environmentalLayers <- subset(environmentalLayers, c("TemperatureMax","TemperatureMin","Precipitation"))

# Inspect spatial autocorrelation in the layers
sacAnalysis <- sac(environmentalLayers)

sacAnalysis$plot
sacDist <- sacAnalysis$distance

## -----------------------
# generate pseudo-absence and combine presences and background into a unique object
absences <- pseudoAbsences(rasterLayers=environmentalLayers,presences, n=1147)
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
monotonicConstrains <- c("TemperatureMax" = -1,"TemperatureMin" = 1,"Precipitation" = 1)

# fit brt model
model <- trainModel(modelData, algorithm ="brt", hyperparameters=hyperparametersList, monotonicity=monotonicConstrains, reclassCriteria="minimumTrain95")

# inspect hyperparameter combination selected per fold
model$hyperparameters

# inspect model performance
model$performance

# inspect relative contribution of variables
modelVarContrib <- variableContribution(model,rasterLayers=environmentalLayers)
modelVarContrib$dataFrame
modelVarContrib$plot

# inspect response curves of variables
modelPartialPlot <- partialPlot(model,modelData,variablePlot="TemperatureMax")
modelPartialPlot$tippingPoints
modelPartialPlot$partialPlot

modelPartialPlot <- partialPlot(model,modelData,variablePlot="TemperatureMin")
modelPartialPlot$tippingPoints
modelPartialPlot$partialPlot

modelPartialPlot <- partialPlot(model,modelData,variablePlot="Precipitation")
modelPartialPlot$tippingPoints
modelPartialPlot$partialPlot

## -----------------------
# predict distribution

# predict the distribution with the model over a raster stack
prediction <- predictModel(model=model,newData=environmentalLayers)

probabilityRaster <- prediction$rasterLayer
uncertaintyRaster <- prediction$rasterLayerSD

# map the predicted distribution
plot(probabilityRaster, main="Predicted species distribution",col = rev(topo.colors(100)))
plot(uncertaintyRaster, main="Uncertainty of predicted species distribution",col = rev(topo.colors(100)))

# reclass the predicted distribution into binary presence/absence using the threshold that maximizes TSS
observed <- records$PA
predicted <- extract(probabilityRaster,records, ID=FALSE)[,1]
modelPerformance(observed,predicted,index="minimumTrain95")

threshold <- 0.81559
rclmat <- data.frame(from=c(0, threshold), to=c(threshold, 1), becomes=c(0,1))
rclmat

predictionPresentReclass <- classify(prediction$rasterLayer, rclmat)
plot(predictionPresentReclass, main="Predicted species distribution",col = c("#b4b295ff", "#0c2f4dff"))

## -----------------------
# model transferability

# download layers for end-of-century conditions
environmentalLayersSSP245 <- cmip6_world(model="ACCESS-ESM1-5", ssp="245", time="2081-2100", var="bioc", download=TRUE, res=5, path="_ rasters")
environmentalLayersSSP370 <- cmip6_world(model="ACCESS-ESM1-5", ssp="370", time="2081-2100", var="bioc", download=TRUE,res=5,path="_ rasters")

# select relevant layers
environmentalLayersSSP245 <- environmentalLayersSSP245[[c(5, 6, 12)]]

# change names for simplicity
names(environmentalLayersSSP245)
names(environmentalLayersSSP245) <- c("TemperatureMax","TemperatureMin","Precipitation")

# select relevant layers
environmentalLayersSSP370 <- environmentalLayersSSP370[[c(5, 6, 12)]]

# change names for simplicity
names(environmentalLayersSSP370)
names(environmentalLayersSSP370) <- c("TemperatureMax","TemperatureMin","Precipitation")

# crop and mask the environmental layers to the extent of the study region
environmentalLayersSSP245 <- crop(environmentalLayersSSP245, myStudyRegion)
environmentalLayersSSP245 <- mask(environmentalLayersSSP245, myStudyRegion)

environmentalLayersSSP370 <- crop(environmentalLayersSSP370, myStudyRegion)
environmentalLayersSSP370 <- mask(environmentalLayersSSP370, myStudyRegion)

# plot the masked environmental layers
plot(environmentalLayersSSP245, axes=TRUE)
plot(environmentalLayersSSP370, axes=TRUE)

# Multivariate environmental similarity surfaces for extrapolation detection

mess_raster <- mess(stack(environmentalLayersSSP245), as.data.frame(environmentalLayers))
mess_raster <- rast(mess_raster)
mess_raster <- mask(mess_raster, myStudyRegion)
plot(mess_raster < 0)


# transfer the distribution model
predictionSSP245 <- predictModel(model=model,newData=environmentalLayersSSP245)
predictionSSP370 <- predictModel(model=model,newData=environmentalLayersSSP370)

# For two plots stacked vertically, we want 3 rows, 1 column.
par(mfrow = c(3, 1))
par(mar = c(4, 4, 2, 1))
plot(prediction$rasterLayer, main="Present-day distribution",col = rev(topo.colors(100)))
plot(predictionSSP245$rasterLayer, main="End-of-century distribution (SSP2-4.5)",col = rev(topo.colors(100)))
plot(predictionSSP370$rasterLayer, main="End-of-century distribution (SSP3-7.0)",col = rev(topo.colors(100)))
dev.off()

predictionDiffSSP245 <- predictionSSP245$rasterLayer - prediction$rasterLayer
predictionDiffSSP370 <- predictionSSP370$rasterLayer - prediction$rasterLayer

# For two plots stacked vertically, we want 3 rows, 1 column.
par(mfrow = c(3, 1))
par(mar = c(4, 4, 2, 1))
plot(prediction$rasterLayer, main="Present-day distribution",col = rev(topo.colors(100)))
plot(predictionDiffSSP245, main="End-of-century difference in suitability (SSP2-4.5)",col = rev(topo.colors(100)))
plot(predictionDiffSSP370, main="End-of-century difference in suitability (SSP3-7.0)",col = rev(topo.colors(100)))
dev.off()

predictionSSP245Reclass <- classify(predictionSSP245$rasterLayer, rclmat)
predictionSSP370Reclass <- classify(predictionSSP370$rasterLayer, rclmat)

par(mfrow = c(3, 1))
par(mar = c(4, 4, 2, 1))
plot(predictionPresentReclass, main="Present-day distribution",col = c("#E1C177", "#043259"))
plot(predictionSSP245Reclass, main="End-of-century distribution (SSP2-4.5)",col = c("#E1C177", "#043259"))
plot(predictionSSP370Reclass, main="End-of-century distribution (SSP3-7.0)",col = c("#E1C177", "#043259"))
dev.off()

## -----------------------
# export final predictions

# save the raster layers to external files

writeRaster(predictionPresentReclass, filename="myFile1.tif", overwrite=TRUE)
writeRaster(predictionSSP245Reclass, filename="myFile2.tif", overwrite=TRUE)
writeRaster(predictionSSP370Reclass, filename="myFile3.tif", overwrite=TRUE)

## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------