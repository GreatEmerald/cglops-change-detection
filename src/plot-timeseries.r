# Script for plotting time series for visualisation
# Based on visualisation.r, but more modular

source("utils/mosaicvi.r")
source("utils/dates.r")
source("utils/cachevi.R")
source("utils/visualise_breakpoints.r")

# Input data
WinterChange = st_read("../data/ValidationPoints-December.csv", options=c("X_POSSIBLE_NAMES=sample_x", "Y_POSSIBLE_NAMES=sample_y"))
ChangeIndices = which(WinterChange$ChangeType == "LC change")

EVIRaster = MosaicVI("EVI")
EVITS = CacheVI(WinterChange, EVIRaster, "EVI")
# Breakpoints
for (i in ChangeIndices)
    visualise_breakpoints(EVITS[i,1:length(dates)], dates = dates, ylab="EVI", main=i, sub=WinterChange$comment[i])

# 577 363 313 311*** 165 124*** 113 107 95 80*** 73 65* 52 51** 47* 43*, 28*, 27**, 22, 18*, 8***, 7
ChangeIndices=c(28, 27, 18, 8, 311, 124, 80, 65, 51, 47, 43)
ChangeIndices=c(311, 124, 80, 8)
# Investigate 483 489

# BFAST. Results differ. This is because BFAST runs breakpoints on deseasonalised and detrended data separately.
for (i in ChangeIndices)
    visualise_breakpoints(EVITS[i,1:length(dates)], dates = dates, bfast=TRUE)

# Pretty labels
for (i in ChangeIndices)
    visualise_breakpoints(EVITS[i,1:length(dates)], dates = dates, ylab="EVI", main=paste("Long:", WinterChange$sample_x[i], "Lat:", WinterChange$sample_y[i]) , sub=WinterChange$comment[i])
for (i in ChangeIndices)
    visualise_breakpoints(EVITS[i,1:length(dates)], dates = dates, ylab="EVI", sub=paste("Long:", WinterChange$sample_x[i], "Lat:", WinterChange$sample_y[i]))

