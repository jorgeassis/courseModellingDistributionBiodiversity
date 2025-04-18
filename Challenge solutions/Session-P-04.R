# --------------------------------------
# --------------------------------------
#
# Modelling the distribution of marine biodiversity>under global climate change
# Jorge Assis, PhD // Nord University, Norway // biodiversityDS, Centre of Marine Sciences, Portugal
#
# --------------------------------------
# --------------------------------------

# Challenge 4.1
# The objective of this challenge is to extract and plot a set of environmental layers from Bio-ORACLE.

# load libraries
library(terra)
library(bDSSDMTools)

# download present-day layer
variables <- list(var1 = c("thetao", "ltmax", "depthsurf"),
                  var1 = c("thetao", "ltmin", "depthsurf"))
environmentalVariables <- download_multiple_layers(variables = variables,
                                                   experiment = c("baseline"),
                                                   decade = c(2010),
                                                   longitude = c(-15, 45),
                                                   latitude = c(25, 50))
# plot layers
plot(environmentalVariables, main = "Environmental variables", axes = TRUE)

# --------------------------------------
# --------------------------------------

# Challenge 4.2
# The objective of this challenge is to determine the degree of temperature change along the northeast coast of the Pacific Ocean. Consider using two different decades (2010 and 2090) and two scenarios of climate change.

# download present-day layer
variables <- list(var1 = c("thetao", "ltmax", "depthsurf"))

# download present-day layer
temperaturePresent <- download_multiple_layers(variables = variables,
                                                   experiment = c("baseline"),
                                                   decade = c(2010),
                                                   longitude = c(-175, -100),
                                                   latitude = c(17, 65))

# download end-of-century layer
temperatureFuture <- download_multiple_layers(variables = variables,
                                              experiment = c("ssp370"),
                                              decade = c(2090),
                                              longitude = c(-175, -100),
                                              latitude = c(17, 65))

temperatureAnomaly <- temperatureFuture - temperaturePresent

# plot layers
plot(temperatureAnomaly, main = "Temperature anomaly along the NE Pacific Ocean", axes = TRUE)

# --------------------------------------
# --------------------------------------
