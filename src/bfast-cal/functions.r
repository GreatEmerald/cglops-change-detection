# Functions for BFAST calibration

library(sf)
library(raster)
library(lubridate)
library(strucchange)
library(bfast)
library(foreach)
library(doParallel)

source("bfast-cal/plotting.r")
source("bfast-cal/utils.r")
source("bfast-cal/01-preprocess.r")
source("bfast-cal/02-detectbreaks.r")
source("bfast-cal/03-batchprocessing.r")
source("bfast-cal/04-validation.r")
source("utils/enable_fast_bfast.r")
registerDoParallel(cores = 4)

# TODO: test different VIs. Change input to LoadVITS, select which one works best (first under defaults).
TestVIs = function(VIs, ...)
{
    #Result = NULL
    #for (vi in VIs)
    Result = foreach(vi = VIs, .combine=rbind, .multicombine = TRUE, .verbose=TRUE) %dopar%
    {
        VI = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"), vi)
        VI_full = MergeAllYears(VI, RawData)
        # Filter out 2015 and 2018, former can't change and latter can't be detected
        VI_full = VI_full[VI_full$reference_year %in% 2016:2017,]
        # Filter out burnt areas
        # IDs to filter out
        BurntIDs = GetBurntIDs(VI_full)
        BurntIDs = BurntIDs[!BurntIDs %in% ChangeIDs_UTM] # Do not exclude IDs that have changed
        VI_full = VI_full[!VI_full$sample_id %in% BurntIDs,] # Exclude burnt areas
        
        VI_full = TestMODBreakpointDetection(VI_full, plot=FALSE, ...)
        VI_full$changed = VI_full$change_at_300m == "yes"
        VI_full$vi = vi
        VI_full[,c("sample_id", "year_fraction", "changed", "bfast_guess", "vi")]
    }
    return(Result)
}

# Filter out burnt area pixels
# Only needed for the time being
PreprocessBurnt = function()
{
    # Stack 2016 and 2017
    BurntFiles = list.files("../data", pattern = glob2rx("Stratification_Season_WATER-BA-MASK_Africa_201*.tif"), full.names = TRUE)
    BurntStack = stack(BurntFiles)
    # We only care about burnt, so drop the other layers
    BurntStack = dropLayer(BurntStack, c(1,3))
    # Crop to our area
    BurntStack = crop(BurntStack, EVI_8d_16int)
    # This is terribly slow, so use a VRT
    system(paste("gdalbuildvrt -b 2 -separate ../data/burnt.vrt", paste(BurntFiles, collapse=" ")))
    BurntStack = brick("../data/burnt3.vrt")
    # Replace values so that 0 is not burnt and 1 is burnt, and
    # Combine both together
    BurntLayer = calc(BurntStack, function(x){ x[is.na(x)] = 0; x[x==255] = 0; x[1] | x[2] }, progress="text")
    # Still takes forever, so instead do two gadlbuildvrts without -separate, to set 255 and 0 to NA
    # and then to set no NA to have it set to 0, then gdalwarp to not have R crash
    BurntLayer = raster("../data/burnt3.tif")
    # Pixels are 100m, whereas we are concerned with 300m, so run a filter
    Burnt300m = focal(BurntLayer, w=matrix(1,3,3), fun=modal, filename="../data/burnt-mask.tif", progress="text", options=c("COMPRESS=DEFLATE", "NUM_THREADS=4"))
    # ...which is again too slow, so use r.neighbours from GRASS
    Burnt300m = raster("../data/burnt-grass.tif")
    return(Burnt300m)
}

GetBurntIndices = function(df, burnt=raster("../data/burnt-grass.tif"))
{
    return(which(as.logical(extract(burnt, df))))
}

GetBurntVector = function(df, burnt=raster("../data/burnt-grass.tif"))
{
    return(as.logical(extract(burnt, df)))
}

GetBurntIDs = function(df, burnt=raster("../data/burnt-grass.tif"))
{
    return(df$sample_id[as.logical(extract(burnt, df))])
}
