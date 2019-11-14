# Script for determining which BFAST/change detection method works the best,
# according to change validation data.

library(sf)
library(raster)
library(lubridate)
library(strucchange)
library(bfast)
library(foreach)
library(doParallel)

source("utils/enable_fast_bfast.r")

# Import the validation CSV; it contains one entry per year that gives the approximate time of the year.

# 1) Load the CSV. Ideally we only need an st_read() to get an sf object
LoadReferenceData = function(path="../data/ValidationPoints-AfricaPriority.csv")
{
    # Normally if we have a good CSV, we just do this.
    #Data = st_read(path, options=c("X_POSSIBLE_NAMES=x", "Y_POSSIBLE_NAMES=y"))
    # Now our centroids are in another file, so we need to do more.
    DataBadCoords = read.csv(path)
    CentroidData = read.csv("../data/2019_11_06_training_data_100m.csv")
    Data = merge(DataBadCoords, CentroidData, by.x="sample_id", by.y="sampleid")
    Data = st_as_sf(Data, coords = c("centroid_x", "centroid_y"), dim="XY")
    Data = Data[!is.na(Data$change_at_300m),] # Remove NAs
    
    # Which columns are numeric
    NumCols = c("rowid", "location_id", "sample_id", "bare", "burnt", "crops",
        "fallow_shifting_cultivation", "grassland", "lichen_and_moss", "shrub",
        "snow_and_ice", "tree", "urban_built_up", "water", "wetland_herbaceous",
        "not_sure", "reference_year")
    for (ColName in NumCols)
    {
        Data[[ColName]] = as.numeric(Data[[ColName]])
    }
    st_crs(Data) = 4326
    return(Data)
}

# Calculate dates from number of elements in the input
# Returns a ts object
GetTS = function(data, frequency=NULL, start=2009, years=10)
{
    stopifnot(is.vector(data)) # Univariate only
    # 8-daily frequency is 46, 16-daily is 23, we have 10 years of data
    if (is.null(frequency))
        frequency = length(data)/years
    return(ts(data, start=start, frequency = frequency))
}

# Returns a Date object
GetDates = function(...)
{
    TS = GetTS(...)
    return(as.Date(date_decimal(as.numeric(time(TS)))))
}

GetDates8d = function(...)
{
    return(GetDates(1:460))
}

GetDates16d = function(...)
{
    return(GetDates(1:230))
}

# Utility to extract only the time series table from an input full sf object
GetMatrixFromSF = function(sf)
{
    # Assumes that the pattern of columns is XNNNN.NN.NN
    TSNames = grep(glob2rx("X????.??.??"), names(sf), value = TRUE)
    as.matrix(as.data.frame(sf)[, TSNames])
}

# 2) Extract time series data from the coordinates,
# and cache it in a CSV/GPKG so that we don't need to do that again.
# This is where we select different VIs.
# TODO: Check the case when a point is in more than one UTM zone
# The output is an sf data.frame, with rows being unique points,
# and columns being timesteps, first columns being x, y, and sample_id, last being geometry.
# Data is the value of the vegetation index.
LoadVITS = function(pointlocs, vi="EVI_8d_Int16", sourcedir="/data/users/Public/greatemerald/modis-utm/input-vrt/", prefix="", force=FALSE)
{
    # Cache file. Has location_id, sample_id, x, y, geometry and the values
    VITSFile = paste0("../data/", prefix, vi, "-TS.gpkg")
    if (force || !file.exists(VITSFile))
    {
        print(paste("Cache file", VITSFile, "not found, generating..."))
        # Deduplicate the input. We shouldn't extract from the same point more than once
        # sample_id is unique, but also an option is to dedup on x/y
        pointlocs = pointlocs[!duplicated(pointlocs$sample_id),]
        
        # Input directory that contains all our VRTs is sourcedir,
        # inside we have VIs, and then UTM zones as VRTs
        InputVRTs = list.files(file.path(sourcedir, vi), full.names = TRUE)
        
        OutDF = NULL
        for (UTMfile in InputVRTs)
        {
            print(paste("Processing", UTMfile))
            VIMosaic = brick(UTMfile)
            VIMosaic = setZ(VIMosaic, GetDates(1:nlayers(VIMosaic)))
            names(VIMosaic) = GetDates(1:nlayers(VIMosaic))
            # Reproject all points to this UTM zone
            PointsUTM = st_transform(pointlocs, crs(VIMosaic))
            # Which ones are inside the UTM zone?
            PointsInZone = st_contains(st_as_sfc(st_bbox(VIMosaic)), PointsUTM)[[1]]
            if (length(PointsInZone) <= 0)
                next
            # Extract those
            print(system.time(ChangedVITS <- extract(VIMosaic, PointsUTM[PointsInZone,]))) # method="bilinear" takes too much RAM
            # Put back into sf with orginal coords
            OutDF = rbind(OutDF, cbind(pointlocs[PointsInZone,c("x", "y", "sample_id")], ChangedVITS))
        }
        
        # Finally, save to cache and not do that again
        st_write(OutDF, VITSFile, delete_dsn = TRUE)
    } else OutDF = st_read(VITSFile)
    
    NALocations = apply(GetMatrixFromSF(OutDF), 1, function(x) all(is.na(x)))
    OutDF = OutDF[!NALocations,] # Remove all that are only NAs
    OutDF = OutDF[!duplicated(OutDF$sample_id),] # Deduplicate
    return(OutDF)
}

# Merge the matrix of unique points with the reference data, to get all years info
MergeAllYears = function(df, data)
{
    return(merge(df, as.data.frame(data), "sample_id"))
}

# Example:
#MyEVI = LoadVITS(LoadReferenceData())
# plot(GetDataFromSF(MyEVI)[2,]~GetDates(1:ncol(GetDataFromSF(MyEVI))), type="l")
# bfts = GetTS(c(GetMatrixFromSF(MyEVI[2,])))
# bf = breakpoints(formula(response~trend+harmon), data=bfastpp(bfts), h=46)

RawData = LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv")
ChangeIDs = as.data.frame(RawData)[RawData$change_at_300m == "yes","sample_id"]
length(ChangeIDs) # 75: total number of (non-unique) change points
length(unique(ChangeIDs)) # 62: unique change locations
plot(RawData[RawData$change_at_300m == "yes","sample_id"]) # All Africa

EVI_8d_16int = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"))
plot(GetMatrixFromSF(EVI_8d_16int)[1,]~GetDates8d(), type="l")

EVI_8d_16int_full = MergeAllYears(EVI_8d_16int, RawData)
ChangeIDs_UTM = as.data.frame(EVI_8d_16int_full)[EVI_8d_16int_full$change_at_300m == "yes","sample_id"]
length(ChangeIDs_UTM) # 50: total number of points in UTM
length(unique(ChangeIDs_UTM)) # 39: unique change locations
length(unique(as.data.frame(EVI_8d_16int_full)[,"sample_id"])) # out of 1900
plot(EVI_8d_16int_full[EVI_8d_16int_full$change_at_300m == "yes","sample_id"]) # UTM zones only

NDMI_8d_16int = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"), "NDMI_8d_Int16")
NDMI_8d_16int_full = MergeAllYears(NDMI_8d_16int, RawData)
plot(GetMatrixFromSF(NDMI_8d_16int)[1,]~GetDates8d(), type="l")

NDMI_16d_16int = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"), "NDMI_16d_Int16")
NDMI_16d_16int_full = MergeAllYears(NDMI_16d_16int, RawData)

NIRv_8d_16int = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"), "NIRv_8d_Int16")
NIRv_8d_16int_full = MergeAllYears(NIRv_8d_16int, RawData)
plot(GetMatrixFromSF(NIRv_8d_16int)[1,]~GetDates8d(), type="l")

NIRv_16d_16int = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"), "NIRv_16d_Int16")
NIRv_16d_16int_full = MergeAllYears(NIRv_16d_16int, RawData)

NIRv_16d_byte = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"), "NIRv_16d_Byte")
NIRv_16d_byte_full = MergeAllYears(NIRv_16d_byte, RawData)

# 3) Run BFAST over the time series and get the detected breaks.
# This is the function that should be optimised;
# optional arguments should include all the tunables.
# Ideally we would run this function in an optimiser (?optim),
# so it tries all combinations of the variables and finds the best one.
# TODO: "dates" should be an argument; we need to define somehow how we want to store time.
#       BFAST works with ts(), so instead of "dates" it's better to use start and frequency.
#       e.g. ts(data, start=c(StartYear, 1), frequency=Frequency) where
#       StartYear = 2009, Frequency = 23 for 16-daily composites, 46 for 8-daily composites
# TODO: Set h = Frequency (this is a year)
# TODO: Vectorise; as input it should take the reference data.
#       Output would be a list/matrix of breaks (one column per point),
#       or maybe better a data.frame with X, Y, break time (fractional year).
# TODO: Optionally, also allow plotting, but that only makes sense for single inputs
#       (don't want 1000 plots)
# TODO: Try to add more tunables to BFAST0N
# The output is fractional years of all detected breaks, or FALSE if none detected,
# or NA if not enough observations/error in running the function.
MODDetectBreaks = function(InputTS, scrange=c(2009, 2019), scsig=0.05, breaks="LWZ", sctype="OLS-MOSUM", plot=FALSE, quiet=FALSE)
{
    # The input should be a single row of a matrix.
    InputTS = GetTS(InputTS)
    Observations = sum(!is.na(InputTS))
    if (!quiet)
        print(paste("Observations for point:", Observations))
    h = frequency(InputTS)
    if (Observations > h*2) {
        bpp = bfastpp(InputTS, order=3)
        if (!is.null(scrange)) {
            SC = sctest(efp(response ~ (harmon + trend),
                            data=bpp[bpp$time > scrange[1] & bpp$time < scrange[2],],
                            h=h/Observations, sctype=type))$p.value > scsig
            if (!is.null(SC) && !is.na(SC) && SC) return(FALSE)
        }
        
        bp = try(breakpoints(formula(response~trend+harmon), data=bpp, h=h))
        
        if ("try-error" %in% class(bp)) {
            print("An error has occurred:")
            print(bp)
            return(NA)
        }
        bpOptim = breakpoints(bp, breaks=breaks)
        
        if (plot)
        {
            plot(InputTS)
            lines(fitted(bp, breaks=breaks)~bpp$time, col="green")
            if (length(bpOptim$breakpoints) > 0) {
                abline(v=bpp[bpOptim$breakpoints, "time"], col="red")
            }
        }
        
        if (all(is.na(bpOptim$breakpoints)))
            return(FALSE)
        if (!quiet)
            print(bpp[bpOptim$breakpoints, "time"])
        return(bpp[bpOptim$breakpoints, "time"])
    } else {
        print("too cloudy")
        return(NA)
    }
}

MODDetectBreaks(GetMatrixFromSF(EVI_8d_16int)[EVI_8d_16int$sample_id == unique(ChangeIDs_UTM)[1],], breaks="LWZ", plot=TRUE)
MODDetectBreaks(GetMatrixFromSF(NDMI_8d_16int)[NDMI_8d_16int$sample_id == unique(ChangeIDs_UTM)[1],], breaks="LWZ", plot=TRUE)
MODDetectBreaks(GetMatrixFromSF(NIRv_16d_16int)[NIRv_16d_16int$sample_id == unique(ChangeIDs_UTM)[1],], breaks="LWZ", plot=TRUE)
MODDetectBreaks(GetMatrixFromSF(NIRv_16d_16int)[NIRv_16d_16int$sample_id == unique(ChangeIDs_UTM)[1],], breaks="BIC", plot=TRUE)

# 3b) t-test
TestMODttest = function(i, ChangedVITS, TargetYears=AllTargetYears, sig=0.05)
{
    Point1TS = ChangedVITS[i, ]
    TargetYear = TargetYears[i]
    if (is.na(TargetYear))
        TargetYear = 2016
    start = as.Date(paste(TargetYear, "01", "01", sep="-"))
    end = start+366
    Observations = sum(!is.na(Point1TS))
    print(paste("Observations for point:", Observations))
    if (Observations > 42) {
        Point1Mean = mean(Point1TS[dates < start], na.rm=TRUE)
        P1Remainder = Point1TS[dates > start & dates < end]
        result = t.test(P1Remainder, mu=Point1Mean)$p.value
        return(result < sig)
        # Might be useful to t-test the other params too
        #P1TS = bfastts(Point1TS, dates, "10-day")
        #P1PP = bfastpp(P1TS, order=3)
        
    } else {
        print("too cloudy")
        return(NA)
    }
}

# 3c) BFAST Monitor
TestMODMonitor = function(i, ChangedVITS, TargetYears=AllTargetYears, threshold=0.25, cloud_threshold=42, ...)
{
    MyDates = dates
    Point1TS = ChangedVITS[i, ]
    TargetYear = TargetYears[i]
    if (is.na(TargetYear))
        TargetYear = 2016
    Observations = sum(!is.na(Point1TS))
    print(paste("Observations for point:", Observations))
    #plot(Point1TS~MyDates, type="l"); abline(v=as.POSIXct(as.Date("2017-01-01")))
    if (Observations > cloud_threshold) {
        P1TS = bfastts(Point1TS, MyDates, "10-day")
        P1BM = bfastmonitor(P1TS, TargetYear-threshold, ...)
        if (is.na(P1BM$breakpoint))
            return(FALSE)
        return(P1BM$breakpoint < TargetYear+1+threshold)
    } else {
        print("too cloudy")
        return(NA)
    }
}

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

MyBreak = MODDetectBreaks(GetMatrixFromSF(NIRv_8d_16int)[NIRv_8d_16int$sample_id == unique(ChangeIDs_UTM)[2],], breaks="BIC", plot=TRUE)
MyReference = as.data.frame(RawData)[RawData$sample_id==unique(ChangeIDs_UTM)[2] & RawData$change_at_300m=="yes","year_fraction"]
abline(v=MyReference, col="blue")
IsBreakInTargetYear(MyBreak, MyReference)

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

EVI_8d_16int_defaults = TestMODBreakpointDetection(EVI_8d_16int_full[EVI_8d_16int_full$sample_id %in% ChangeIDs_UTM,], breaks="BIC", plot=FALSE)
EVI_8d_16int_LWZ = TestMODBreakpointDetection(EVI_8d_16int_full[EVI_8d_16int_full$sample_id %in% ChangeIDs_UTM,], breaks="LWZ", plot=FALSE)
NDMI_8d_16int_defaults = TestMODBreakpointDetection(NDMI_8d_16int_full[NDMI_8d_16int_full$sample_id %in% ChangeIDs_UTM,], breaks="BIC", plot=FALSE)
NDMI_8d_16int_LWZ = TestMODBreakpointDetection(NDMI_8d_16int_full[NDMI_8d_16int_full$sample_id %in% ChangeIDs_UTM,], breaks="LWZ", plot=FALSE)
NIRv_8d_16int_defaults = TestMODBreakpointDetection(NIRv_8d_16int_full[NIRv_8d_16int_full$sample_id %in% ChangeIDs_UTM,], breaks="BIC", plot=FALSE)
NIRv_8d_16int_LWZ = TestMODBreakpointDetection(NIRv_8d_16int_full[NIRv_8d_16int_full$sample_id %in% ChangeIDs_UTM,], breaks="LWZ", plot=FALSE)

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
    TruePositiveRate = TruePositiveCount / sum(truth, na.rm=TRUE) # AKA Sensitivity
    Specificity = TrueNegativeCount / (TrueNegativeCount + FalsePositiveCount)
    # Percent of false positive out of no change
    FalsePositiveRate = FalsePositiveCount / sum(!truth, na.rm=TRUE) # False positive rate or alpha or p-value or Type I Error
    SignalToNoise = TruePositiveCount / FalsePositiveCount
    SignalToNoiseRate = TruePositiveRate / FalsePositiveRate # Likelihood Ratio for Positive Tests
    PositivePredictiveValue = TruePositiveCount / (TruePositiveCount + FalsePositiveCount)
    Accuracy = (TruePositiveCount + TrueNegativeCount) / length(truth)
    return(data.frame(TruePositiveCount, FalsePositiveCount, TrueNegativeCount, FalseNegativeCount,
                      TruePositiveRate, Specificity, FalsePositiveRate, SignalToNoise, SignalToNoiseRate,
                      PositivePredictiveValue, Accuracy))
}

FPStats(EVI_8d_16int_defaults)
FPStats(EVI_8d_16int_LWZ) # Better positive predictive value
FPStats(NDMI_8d_16int_defaults)
FPStats(NDMI_8d_16int_LWZ)
FPStats(NIRv_8d_16int_defaults)
FPStats(NIRv_8d_16int_LWZ)

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

registerDoParallel(cores = 4)
# This runs with LWZ
VIResult = TestVIs(c(
    "EVI_16d_Int16", "EVI_8d_Int16", "NDMI_16d_Int16", "NDMI_8d_Int16",
    "NIRv_16d_Byte", "NIRv_16d_Int16", "NIRv_8d_Byte", "NIRv_8d_Int16")
    )
FPStats(VIResult[VIResult$vi=="EVI_16d_Int16",]$bfast_guess, VIResult[TestResult$vi=="EVI_16d_Int16",]$changed)

# Get a table of all statistics
tapply(1:nrow(VIResult), VIResult$vi, function(x)FPStats(VIResult[x,]$bfast_guess, VIResult[x,]$changed))

# Test with BIC
BICResult = TestVIs(c("EVI_16d_Int16", "NDMI_16d_Int16", "NIRv_16d_Byte", "NIRv_16d_Int16"), breaks="BIC")
tapply(1:nrow(BICResult), BICResult$vi, function(x)FPStats(BICResult[x,]$bfast_guess, BICResult[x,]$changed))

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

EVI_unburnt = EVI_8d_16int_full[-GetBurntIndices(EVI_8d_16int_full),]
ChangeIDs_UTM_unburnt = as.data.frame(EVI_unburnt)[EVI_unburnt$change_at_300m == "yes","sample_id"]
length(unique(ChangeIDs_UTM_unburnt)) # 32
IDs_burnt = GetBurntIDs(EVI_8d_16int_full)

# TODO: make an optimisation routine where we run BFAST with different parameters,
# get statistics and choose the best values based on the statistics. Maybe using optim().
# Caveat that if we do it automatically, we may run into edge cases,
# i.e. all NA and one correct break + one correct no break is best.

BICResult = TestVIs(c("EVI_16d_Int16", "NDMI_16d_Int16", "NIRv_16d_Byte", "NIRv_16d_Int16"), breaks="BIC")
tapply(1:nrow(BICResult), BICResult$vi, function(x)FPStats(BICResult[x,]$bfast_guess, BICResult[x,]$changed))

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

SCTestTest = TestParams(list(
    list(breaks="LWZ", scsig=0.25), # Sensitivity 0.44, specificity 0.90, a bit better
    list(breaks="LWZ", scsig=0.75), # Sensitivity 0.44, specificity 0.90
    list(breaks="BIC", scsig=0.01), # Sensitivity 0.56, specificity 0.39, a bit better
    list(breaks="BIC", scsig=0.001) # All non-break
    ), "NIRv_16d_Byte")
tapply(1:nrow(SCTestTest), SCTestTest$call, function(x)FPStats(SCTestTest[x,]$bfast_guess, SCTestTest[x,]$changed))

SCTestTest = TestParams(list(
    list(breaks="LWZ", scsig=0.15), # Sensitivity 0.41, specificity 0.90, worse than p<0.25
    list(breaks="BIC", scsig=0.005), # All non-break
    list(breaks="BIC", scsig=0.007), # All non-break
    list(breaks="BIC", scsig=0.003) # All non-break
), "NIRv_16d_Byte")
tapply(1:nrow(SCTestTest), SCTestTest$call, function(x)FPStats(SCTestTest[x,]$bfast_guess, SCTestTest[x,]$changed))

LWZSubset = SCTestTest[SCTestTest$call == "breaks LWZ, scsig 0.15",]
FNSamples = unique(LWZSubset[LWZSubset$changed & !LWZSubset$bfast_guess,]$sample_id)
for (i in 1:length(FNSamples))
{
    MODDetectBreaks(c(GetMatrixFromSF(NIRv_16d_byte[NIRv_16d_byte$sample_id == FNSamples[i],])), scrange=NULL, plot=TRUE)
    abline(v=LWZSubset[LWZSubset$changed & LWZSubset$sample_id == FNSamples[i],]$year_fraction, col="blue")
    title(FNSamples[i])
} # 1363551 is a mismatch of break definition
# A lot of cases like 1363169 will get detected properly in the future

# Try scrange

SCRangeTest = TestParams(list(
    list(breaks="LWZ", scrange=c(2015, 2019)), # Sensitivity 0.29, specificity 0.94, misses a lot
    list(breaks="BIC", scrange=c(2015, 2019))  # Sensitivity 0.44, specificity 0.60, it's a bad version of LWZ
), "NIRv_16d_Byte")

Result = SCRangeTest
tapply(1:nrow(Result), Result$call, function(x)FPStats(Result[x,]$bfast_guess, Result[x,]$changed))

SCRangeTest = TestParams(list(
    list(breaks="LWZ", scrange=c(2016, 2018)), # Sensitivity 0.29, specificity 0.95, misses a lot
    list(breaks="BIC", scrange=c(2016, 2018))  # Sensitivity 0.32, specificity 0.66, even worse
), "NIRv_16d_Byte")

Result = SCRangeTest
tapply(1:nrow(Result), Result$call, function(x)FPStats(Result[x,]$bfast_guess, Result[x,]$changed))

# Try different test types

SCTypeTest = TestParams(list(
    list(breaks="LWZ", sctype="Rec-CUSUM") # Sensitivity 0.38, specificity 0.94, it's even more specific
    ,list(breaks="BIC", sctype="Rec-CUSUM") # Sensitivity 0.47, specificity 0.62, better than limiting the range
), "NIRv_16d_Byte")

Result = SCTypeTest
tapply(1:nrow(Result), Result$call, function(x)FPStats(Result[x,]$bfast_guess, Result[x,]$changed))

SCTypeTest = TestParams(list(
    list(breaks="LWZ", sctype="RE") # Sensitivity 0.38, specificity 0.94, no change
    ,list(breaks="BIC", sctype="RE") # Sensitivity 0.47, specificity 0.62, no change
), "NIRv_16d_Byte")

Result = SCTypeTest
tapply(1:nrow(Result), Result$call, function(x)FPStats(Result[x,]$bfast_guess, Result[x,]$changed))

SCTypeTest = TestParams(list(
    list(breaks="LWZ", sctype="Score-CUSUM") # Sensitivity 0.38, specificity 0.94, no change
    ,list(breaks="BIC", sctype="Score-CUSUM") # Sensitivity 0.47, specificity 0.62, no change
), "NIRv_16d_Byte")

Result = SCTypeTest
tapply(1:nrow(Result), Result$call, function(x)FPStats(Result[x,]$bfast_guess, Result[x,]$changed))

SCTypePTest = TestParams(list(
    list(breaks="BIC", sctype="Rec-CUSUM", scsig=0.2) # Sensitivity 0.53, specificity 0.46, balanced but still not great, similar to OLS-MOSUM with 0.01
    ,list(breaks="BIC", sctype="Rec-CUSUM", scsig=0.1) # Sensitivity 0.47, specificity 0.55, worse
), "NIRv_16d_Byte")

Result = SCTypePTest
tapply(1:nrow(Result), Result$call, function(x)FPStats(Result[x,]$bfast_guess, Result[x,]$changed))
