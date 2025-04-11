# --------------------------------------
# --------------------------------------
#
# Modelling the distribution of marine biodiversity>under global climate change
# Jorge Assis, PhD // Nord University, Norway // biodiversityDS, Centre of Marine Sciences, Portugal
#
# --------------------------------------
# --------------------------------------

# Challenge 3.1
# The main objective of this challenge is to plot biodiversity records from a csv file, for a given region.

# load libraries
library(terra)
# read occurrence records
records <- read.table("../Data/Text delimited/species_presence_records.csv", sep = ";", header = TRUE)
# transform data.frame into a spatial object
records <- vect(records, geom = c("Lon", "Lat"))
# read a landmass polygon
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")
# crop landmass to the european Extent
mediterraneanExtent <- ext(-35, 45, 25, 55)
landmassMediterranean <- crop(landmass, mediterraneanExtent)

# plot the landmass and the records
plot(landmassMediterranean, main = "Species occurrence records", col = "Gray", border = "Gray",
     axes = TRUE)
plot(records, col = "Black", pch = 19, cex = 0.5, add = TRUE)

# --------------------------------------
# --------------------------------------

# Challenge 3.2
# The objective of this challenge is to clean records for a species of interest and produce a plot. A recipe is provided.

# Load necessary libraries
library(terra)
library(leaflet)
library(mapedit)

# get data from GBIF
recordsObis <- getExternalDataObis(taxa = "Laminaria ochroleuca", getCitation = TRUE)
# inspect the top records of the data.frame
head(recordsObis)

# create a base leaflet map
baseMap <- leaflet()
baseMap <- addTiles(baseMap)

# drawing vector polygon in the Viewer pane
myRegion <- drawFeatures(baseMap)

# inspect the generated polygon depicting the distribution of the species
plot(landmass, main = "Global landmass", col = "gray", border = "gray", axes = TRUE)
plot(myRegion, border = "black", col = "black", add = TRUE)
plot(records, col = "red", pch = 19, cex = 0.5, add = TRUE)

# transform the sf object to a terra object
myRegion <- vect(myRegion)
# intersect the records with the region
records <- intersect(records, myRegion)
# inspect the class of the object
class(records)

# plot the distribution of the species without records outside the known distribution
plot(landmass, main = "Global landmass", col = "gray", border = "gray", axes = TRUE)
plot(records, color = "black", pch = 16, cex = 0.5, add = TRUE)

# import bathymetry layer
bathymetry <- rast("../Data/Raster data/Bathymetry.tif")
# plot the bathymetry layer and the records of occurrence
plot(bathymetry, main = "Bathymetry", axes = TRUE)
plot(records, color = "black", pch = 16, cex = 0.5, add = TRUE)

# extract the values of the bathymetry based on the records
depthDistribution <- extract(bathymetry, records)
# plot the distribution of the depth records
hist(depthDistribution$Bathymetry, breaks = 500)

# extract the depth values of the records
records <- records[depthDistribution$Bathymetry > -30 & !is.na(depthDistribution$Bathymetry),
]
# extract the depth values of the records
depthDistribution <- extract(bathymetry, records)
# plot the distribution of the depth records
hist(depthDistribution$Bathymetry, breaks = 500)

# plot the final records
plot(landmass, main = "Final dataset of occurrence records", col = "gray", border = "gray", axes = TRUE)
plot(records, color = "black", pch = 16, cex = 0.5, add = TRUE)