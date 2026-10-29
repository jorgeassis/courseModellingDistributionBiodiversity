
# Information
# https://wec.wur.nl/r/spatial/raster-data.html

library(geodata)

# Data present-day (baseline)

bioPresent <- worldclim_global(var = "bio", # "tmin", "tmax", "tavg", "prec", "wind", "vapr", or "bio"
                               res = 10, # resolution: 10, 5, 2.5, or 0.5 (minutes of a degree)
                               path = "_ test")

names(bioPresent)
plot(bioPresent)

# Mean Diurnal Range (Mean of monthly max temp - min temp)

bioPresentSubset <- bioPresent[[2]]
plot(bioPresentSubset)

# ------------

# Future data

bioFuture <- cmip6_world(model="ACCESS-ESM1-5", 
                         ssp="245", 
                         time="2041-2060",
                         var="bioc",
                         download=F,
                         res=10,
                         path="_ test")

# Inspect the SpatRaster object:

bioFutureSubset <- bioFuture[[2]]
plot(bioFutureSubset)

anomaly <- bioFutureSubset - bioPresentSubset
plot(anomaly)
