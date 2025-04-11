# --------------------------------------
# --------------------------------------
#
# Modelling the distribution of marine biodiversity>under global climate change
# Jorge Assis, PhD // Nord University, Norway // biodiversityDS, Centre of Marine Sciences, Portugal
#
# --------------------------------------
# --------------------------------------

# Challenge 2.1
# The main objective of this challenge is to map the distribution of a seagrass species based on a ESRI shapefile. For context, the plot must include the global landmass.

# load terra package
library(terra)

# load a point vector shapefile with terra package
distributionRecords <- vect("../Data/Vector data/Seagrass distribution/cymodocea_nodosa.shp")

# plot point data
plot(distributionRecords, main = "Distribution of a segrass species", col = "Black", ylab = "Latitude",
     xlab = "Longitude")

# load a polygon vector shapefile with terra package
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")

# plot spatial polygon vector with a title
plot(landmass, main = "Global landmass", col = "gray", border = "black", axes = TRUE, ylab = "Latitude",
     xlab = "Longitude")

plot(landmass, main = "Distribution of a segrass species", col = "Gray", border = "Gray", axes = TRUE,
     ylab = "Latitude", xlab = "Longitude")
plot(distributionRecords, col = "Black", pch = 20, size = 0.25, add = TRUE)

# --------------------------------------
# --------------------------------------

# Challenge 2.2
# The main objective of this challenge is to generate vector data in R and export it as an ESRI shapefile.

# Load necessary libraries
library(leaflet)
library(mapedit)

# create a base leaflet map
baseMap <- leaflet()
baseMap <- addTiles(baseMap)

# drawing vector polygon in the Viewer pane
drawnFeatures <- drawFeatures(baseMap)
class(drawnFeatures)

drawnFeatures <- vect(drawnFeatures)
class(drawnFeatures)

# load a polygon vector shapefile with terra package
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")

# plot spatial polygon vector with a title
plot(landmass, col = "Gray", border = "Gray", axes = TRUE,
     ylab = "Latitude", xlab = "Longitude")
plot(drawnFeatures, col = "Black", pch = 20, size = 0.25, add = TRUE)

writeVector(x = drawnFeatures, filename = "directory/output_name.shp", filetype = "ESRI Shapefile",
            overwrite = TRUE)

# --------------------------------------
# --------------------------------------

# Challenge 2.3
# The main objective of this challenge is to crop a global map to the Azores region and plot it.

# read landmass from an R package
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")

# generate an extent object (polygon of the region of interest, order: xmin, xmax, ymin,
# ymax)
azoresExtent <- ext(-38, -16, 32, 44)

# for inspection, plot the landmass and the extent region
plot(landmass, main = "Global landmass", col = "gray", border = "gray", axes = TRUE)
plot(azoresExtent, border = "black", add = TRUE)

# generate a new polygon based on the provided extent
azoresLandmass <- crop(landmass, azoresExtent)

# plot the new polygon
plot(azoresLandmass, main = "Azores region", col = "gray", border = "gray", axes = TRUE,
     ylab = "Latitude", xlab = "Longitude")

# --------------------------------------
# --------------------------------------

# Challenge 2.4
# The main objective of this challenge is to map the maximum global sea surface temperature.

seaSurfaceTemp <- rast("../Data/Raster data/Present Surface OceanTemperature LtMax.tif")
plot(seaSurfaceTemp, axes = TRUE, main = "Global maximum sea temperatures")

# --------------------------------------
# --------------------------------------

# Challenge 2.5
# The main objective of this challenge is to map the maximum global sea surface temperature.

seaSurfaceTemp <- rast("../Data/Raster data/Present Surface OceanTemperature LtMax.tif")

# define the extent to crop (xmin, xmax, ymin, ymax)
africaExtent <- ext(-30, 22, -37.5, 37)

# for inspection, plot the raster layer and the extent region
plot(seaSurfaceTemp, main = "Global sea surface temperatures")
plot(africaExtent, add = TRUE)

# crop with polygon (region)
africaSeaSurfaceTemp <- crop(seaSurfaceTemp, africaExtent)

# plot sst of the region of interest
plot(africaSeaSurfaceTemp, main = "African sea surface temperature")

# --------------------------------------
# --------------------------------------

# Challenge 2.6
# The main objective of this challenge is to map the global mean sea surface temperatures using the app function of terra package.

seaSurfaceTempMin <- rast("../Data/Raster data/Present Surface OceanTemperature LtMin.tif")
seaSurfaceTempMax <- rast("../Data/Raster data/Present Surface OceanTemperature LtMax.tif")

meanSeaSurfaceTemp <- app(c(seaSurfaceTempMin, seaSurfaceTempMax), fun = mean)
plot(meanSeaSurfaceTemp, main = "Mean sea surface temperature")

# --------------------------------------
# --------------------------------------

# Challenge 2.7
# The main objective of this challenge is to map the average sea surface temperatures of the Mediterranean Sea

seaSurfaceTempMin <- rast("../Data/Raster data/Present Surface OceanTemperature LtMin.tif")
seaSurfaceTempMax <- rast("../Data/Raster data/Present Surface OceanTemperature LtMax.tif")

meanSeaSurfaceTemp <- app(c(seaSurfaceTempMin, seaSurfaceTempMax), fun = mean)
plot(meanSeaSurfaceTemp, main = "Mean sea surface temperature")

# define the extent to crop (xmin, xmax, ymin, ymax)
mediterraneanExtent <- ext(-15, 45, 30, 47.5)

# for inspection, plot the raster layer and the extent region
plot(meanSeaSurfaceTemp, main = "Mean sea surface temperatures")
plot(mediterraneanExtent, add = TRUE)

# crop with polygon (region)
mediterraneanSeaSurfaceTemp <- crop(meanSeaSurfaceTemp, mediterraneanExtent)

# plot sst of the region of interest
plot(mediterraneanSeaSurfaceTemp, main = "Mediterranean mean sea surface temperature")

# --------------------------------------
# --------------------------------------

# Challenge 2.8
# The main objective of this challenge is to map the global mean sea surface temperatures between 12 and 25 degrees using the classify function.

seaSurfaceTempMax <- rast("../Data/Raster data/Present Surface OceanTemperature LtMax.tif")

# from-to-becomes matrix
rclmat <- data.frame(from = c(-Inf, 12, 25), to = c(12, 25, Inf), becomes = c(0, 1, 0))
temperatureLayerReclass <- classify(seaSurfaceTempMax, rclmat)
plot(temperatureLayerReclass)

# --------------------------------------
# --------------------------------------

# Challenge 2.9
# The main objective of this challenge is to map the global mean sea surface temperatures between 12 and 25 degrees using the classify function.

library(bDSSDMTools)

myFiles <- c("../Data/Raster data/Present Benthic OceanTemperature LtMax.tif",
             "../Data/Raster data/Present Benthic OceanTemperature LtMin.tif",
             "../Data/Raster data/Present Benthic Salinity Mean.tif",
             "../Data/Raster data/Present Benthic TotalPrimaryProductionPhyto Mean.tif")

# combining the raster layers
environmentLayers <- rast(myFiles)

# change name of layers for simplicity
names(environmentLayers) <- c("Ocean Temperature Max", "Ocean Temperature Min", "Salinity", "Productivity")

# get correlation between pairs of layers
correlations <- getCorrPlot(environmentLayers)

# plot correlations
correlations$plot

# plot correlations
correlations$correlationPairs
