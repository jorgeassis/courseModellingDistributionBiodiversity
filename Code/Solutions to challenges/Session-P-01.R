# --------------------------------------
# --------------------------------------
#
# Modelling the distribution of marine biodiversity>under global climate change
# Jorge Assis, PhD // Nord University, Norway // biodiversityDS, Centre of Marine Sciences, Portugal
#
# --------------------------------------
# --------------------------------------

# Challenge 1.2
# The main objective of this challenge is to import data into an R object from a text file.

# load text file
distributionRecords <- read.csv("Data/Text delimited/species_presence_records.csv", sep=";")
head(distributionRecords)
