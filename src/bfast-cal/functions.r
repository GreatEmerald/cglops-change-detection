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
source("utils/enable_fast_bfast.r")
registerDoParallel(cores = 4)


# 4) Did BFAST predict a break in the given year? TRUE/FALSE.
# BreakTimes: all breaks detected by BFAST.
# TargetYear: time at which we want to test.
# Threshold: how much fuzziness is allowed to still consider a break detected.
# This function needs to be run for every year and every point.
# If there was no break for that year, and BFAST predicted one, should return FALSE.
IsBreakInTargetYear = function(BreakTimes, TargetYear, threshold=0.5)
{
    #if (is.na(TargetYear))
    #    TargetYear = 2016 # If there is no break, we look at whether we predicted a break in 2016
    # TODO: Needs to be updated; previously, lack of break meant that the break time would be set to NA,
    # but now it's no longer the case, it's change_at_300m = FALSE and reference_year=<year>
    return(any(BreakTimes > TargetYear - threshold & BreakTimes < TargetYear+1+threshold))
}

# 4b) Vectorised version (takes a list of MODDetectBreaks() output and a column of target years)
# Returns a column of whether BFSAT predicted the break at that time or not.
VectorisedIsBreakInTargetYear = function(BreakList, threshold=0.5, TY=TargetYears)
{
    i = 1:length(BreakList)
    return(sapply(i, function(i){return(IsBreakInTargetYear(BreakList[[i]], TY[i], threshold=threshold))}))
}

# 4c) Wrapper for running 3+4 in one.
# Result is a logical vector saying how many times we accurately predicted a break,
# and how many times not. One value per year.
# VITS is the full VI TS, all years should be there.
TestMODBreakpointDetection = function(VITS, threshold=0.5, TargetYears=TargetYears, ...)
{
    # Output into a new column
    VITS$bfast_guess = NA
    
    pbi = 0
    pb = txtProgressBar(pbi, length(unique(VITS$sample_id)), style = 3)
    # Detect the breaks in a loop over unique points (don't want to run bfast multiple times)
    for (i in unique(VITS$sample_id))
    {
        SampleMatrix = GetMatrixFromSF(VITS[VITS$sample_id == i,])
        BreakTimes = MODDetectBreaks(SampleMatrix[1,], ..., quiet=TRUE)
        
        pbi = pbi + 1
        setTxtProgressBar(pb, pbi)
        
        if (length(BreakTimes) < 2 && is.na(BreakTimes))
            next # Already set to NA
        for (year in as.data.frame(VITS)[VITS$sample_id == i,"year_fraction"])
        {
            VITS[VITS$sample_id == i & VITS$year_fraction == year, "bfast_guess"] = IsBreakInTargetYear(BreakTimes, year, threshold)
        }
    }
    close(pb)
    
    return(VITS)
}

# 5) Get statistics.
FPStats = function(predictions, truth = NULL)
{
    # If we get a data.frame, try to automagically determine what is predicted and what is true
    if (is.data.frame(predictions))
    {
        truth = predictions$change_at_300m == "yes"
        predictions = predictions$bfast_guess
    }
    
    # We predicted a break and it was a break
    TruePositiveCount = sum(predictions & truth, na.rm=TRUE)
    # We predicted a break but there wasn't one (overprediction)
    FalsePositiveCount = sum(predictions & !truth, na.rm=TRUE)
    # We predicted no break, and there were none
    TrueNegativeCount = sum(!predictions & !truth, na.rm=TRUE)
    # We predicted no break, but there was one (we missed it)
    FalseNegativeCount = sum(!predictions & truth, na.rm=TRUE)
    # Percent of true positives out of all change
    Sensitivity = TruePositiveCount / (TruePositiveCount + FalseNegativeCount) # Previously TruePositiveRate
    Specificity = TrueNegativeCount / (TrueNegativeCount + FalsePositiveCount)
    # Percent of false positive out of no change
    FalsePositiveRate = FalsePositiveCount / sum(!truth, na.rm=TRUE) # False positive rate or alpha or p-value or Type I Error
    PositiveProportion = TruePositiveCount / FalsePositiveCount
    PositiveLikelihood = TruePositiveRate / FalsePositiveRate # Likelihood Ratio for Positive Tests
    PositivePredictiveValue = TruePositiveCount / (TruePositiveCount + FalsePositiveCount)
    Accuracy = (TruePositiveCount + TrueNegativeCount) / length(truth)
    return(data.frame(TruePositiveCount, FalsePositiveCount, TrueNegativeCount, FalseNegativeCount,
                      TruePositiveRate, Specificity, FalsePositiveRate, SignalToNoise, SignalToNoiseRate,
                      PositivePredictiveValue, Accuracy))
}

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

# ParamLists is a list of lists of parameter and value pairs
TestParams = function(ParamLists, vi)
{
    Result = foreach(ParamList = ParamLists, .combine=rbind, .multicombine = TRUE, .verbose=TRUE) %dopar%
    {
        VI = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"), vi)
        VI_full = MergeAllYears(VI, RawData)
        rm(VI)
        # Filter out 2015 and 2018, former can't change and latter can't be detected
        VI_full = VI_full[VI_full$reference_year %in% 2016:2017,]
        # Filter out burnt areas
        # IDs to filter out
        BurntIDs = GetBurntIDs(VI_full)
        BurntIDs = BurntIDs[!BurntIDs %in% ChangeIDs_UTM] # Do not exclude IDs that have changed
        VI_full = VI_full[!VI_full$sample_id %in% BurntIDs,] # Exclude burnt areas
        
        callstring = paste(names(ParamList), ParamList, collapse=", ")
        ParamList$plot = FALSE
        ParamList$VITS = VI_full
        rm(VI_full)
        gc()
        Result = do.call("TestMODBreakpointDetection", ParamList)
        Result$changed = Result$change_at_300m == "yes"
        Result$vi = as.factor(vi) # Just for reference
        Result$call = callstring
        Result[,c("sample_id", "year_fraction", "changed", "bfast_guess", "call", "vi")]
    }
    return(Result)
}
