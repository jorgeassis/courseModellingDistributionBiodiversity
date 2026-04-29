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

## Recipe 1: Download and Clean Biodiversity Records. Integrate with text delimented data.
## Collect species occurrences from GBIF, clean and filter the data, visualize distributions, and export cleaned records.

# 1. Set up the project environment
# Create a new R project and R script file to organize and document the workflow.

# 2. Load custom functions
# Source a local R script (sourceFunctions.R) containing user-defined functions for data access and processing.

# 3. Download species records
# Use functions to retrieve biodiversity records for a target species from GBIF databases.

# 4. Combine and visualize raw records
# Merge datasets and plot occurrence points to inspect overall spatial distribution.

# 5. Clean data
# Remove records with missing coordinates and eliminate duplicated entries to ensure data quality.

# 6. Convert to spatial format
# Transform the cleaned data frame into a spatial object with geographic coordinates.

# 7. Visualize records with landmass
# Load a shapefile of landmass and overlay species records to contextualize their geographic distribution.

# 8. Filter by known geographic range
# Draw a polygon representing the species’ known distribution and keep only records within this region.

# 9. Filter by known evelation range
# Use an elevation raster to remove records outside the species’ known elevation range (e.g., higher than 500m).

# 10. Plot final cleaned records
# Visualize the final dataset with base R plots and a refined ggplot2 visualization using tidyterra.

# 11. Export the final dataset
# Save the cleaned and filtered records to a CSV file for later use or modeling.

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
library(ggplot2)
library(tidyterra)
library(leaflet)
library(mapedit)
library(plyr)
library(bDSSDMTools)

## -----------------------
# 01. get records

# get records from external databases
records <- getExternalDataGbif(taxa="Quercus lusitanica",getCitation=TRUE)
names(records)

# Plot the records
plot( records[,c("Lon","Lat")], col="Black" , pch= 19, cex=0.5)

## -----------------------
# 02. remove NAs and duplicated records

# remove NA coordinates
sum( is.na(records$Lon) | is.na(records$Lat) )
records <- records[ ! is.na(records$Lon) & ! is.na(records$Lat) ,]

# remove duplicated records
sum(duplicated(records))
records <- records[!duplicated(records),]

# inspect the number of records
nrow(records)

## -----------------------
# 03. transform data.frame into a spatial object

# transform data.frame into a spatial object
records <- vect(records, geom=c("Lon","Lat"), crs="epsg:4326")

# plot the landmass and the records
landmass <- vect("../Data/Vector data/Landmass/landmass.shp")
plot(landmass, main="Species occurrence records", col='Gray', border='Gray', axes=TRUE)
plot( records, col="Black" , pch= 19, cex=0.5, add=TRUE)

# Keep only points on land
nrow(records)
records <- intersect(records, landmass)
nrow(records)

# plot the landmass and the records
plot(landmass, main="Species occurrence records", col='Gray', border='Gray', axes=TRUE)
plot( records, col="Black" , pch= 19, cex=0.5, add=TRUE)

## -----------------------
# 04. confirm that all records belong to the know geographic distribution of the species

# create a base leaflet map 
baseMap <- leaflet()
baseMap <- addTiles(baseMap)

# drawing vector polygon in the Viewer pane
myRegion <- drawFeatures(baseMap)

# inspect the generated polygon depicting the distribution of the species
plot(landmass, main = "Global landmass", col="gray", border="gray", axes = TRUE)
plot(myRegion, border="black", col="black", add=TRUE)
plot( records, col="red" , pch= 19, cex=0.5, add=TRUE)

class(myRegion)

# transform the sf object to a terra object
myRegion <- vect(myRegion)
# intersect the records with the region
records <- intersect(records, myRegion)

# inspect the distribution of the species
plot(landmass, main = "Global landmass", col="gray", border="gray", axes = TRUE)
plot( records, col="black" , pch= 19, cex=0.5, add=TRUE)

## -----------------------
# 05. confirm that all records belong to the know elevational distribution of the species

# import elevation layer
elevation <- rast("../Data/Raster data/Elevation.tif")
names(elevation) <- "Elevation"

# plot the elevation layer and the records of occurrence
plot(elevation, main = "Elevation", axes=TRUE)
plot(records, col="black", pch=16, cex=0.5,add=TRUE)

# extract the values of the elevation based on the records
elevationDistribution <- extract(elevation, records)
# plot the distribution of the depth records
hist(elevationDistribution$Elevation, breaks=500)

# extract the elevation values of the records (for a known distribution of 0 to 800 meters elevation)
records <- records[elevationDistribution$Elevation >= 0  & elevationDistribution$Elevation < 800 ,]
# extract the depth values of the records
elevationDistribution <- extract(elevation, records)
# plot the distribution of the depth records
hist(elevationDistribution$Elevation, breaks=500)

## -----------------------
# 06. plot final records

# plot records with plot function
plot(landmass, main = "Final dataset of occurrence records", col="gray", border="gray", axes = TRUE)
plot(records, col="black",pch=16, cex=0.5,add=TRUE)

# plot records with ggplot

ggplot() +
  geom_spatvector(data = landmass, fill = "gray", color = "gray") +
  geom_spatvector(data = records, color = "black", size=1.5) +
  labs(title = "Species occurrence records") +
  xlab("Longitude") + ylab("Latitude") + theme_bw()

ggplot() +
  geom_spatvector(data = landmass, fill = "#d5d4d0ff", colour = "#abababff", size = 0.2) +
  geom_spatvector(data = records, color = "#6d2801ff", size = 0.5) +
  scale_y_continuous(breaks = seq(-90,90, by=20), limits = c(30, 50)) +
  scale_x_continuous(breaks = seq(-180,180,by=20), limits = c(-15, 15)) +
  labs(title = "Species occurrence records") +
  xlab("Longitude") + ylab("Latitude") + theme_bw()

## -----------------------
# 07. export final records

writeVector(x = records, filename = "myRecords.gpkg", filetype = "GPKG", overwrite = TRUE)

finalRecords <- data.frame( Lon = geom(records)[,c("x")] , Lat = geom(records)[,c("y")] , records )

# save data frame to external file
write.table(finalRecords,file="myRecords3.csv",sep=";")

## -----------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------
