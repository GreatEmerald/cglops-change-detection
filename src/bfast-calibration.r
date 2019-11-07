# Script for determining which BFAST/change detection method works the best,
# according to change validation data.

library(sf)
library(raster)
library(lubridate)
library(strucchange)
library(bfast)

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
    # Assumes 4 extra columns!
    as.matrix(as.data.frame(sf)[, paste0("X", gsub("-", ".", GetDates(1:(ncol(sf)-4))))])
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
            print(system.time(ChangedVITS <- extract(VIMosaic, PointsUTM[PointsInZone,])))
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

EVI_8d_16int_full = merge(EVI_8d_16int, as.data.frame(RawData), "sample_id")
ChangeIDs_UTM = as.data.frame(EVI_8d_16int_full)[EVI_8d_16int_full$change_at_300m == "yes","sample_id"]
length(ChangeIDs_UTM) # 50: total number of points in UTM
length(unique(ChangeIDs_UTM)) # 39: unique change locations
length(unique(as.data.frame(EVI_8d_16int_full)[,"sample_id"])) # out of 1900
plot(EVI_8d_16int_full[EVI_8d_16int_full$change_at_300m == "yes","sample_id"]) # UTM zones only

NDMI_8d_16int = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"), "NDMI_8d_Int16", force=TRUE)
plot(GetMatrixFromSF(NDMI_8d_16int)[1,]~GetDates8d(), type="l")

NIRv_8d_16int = LoadVITS(LoadReferenceData("../data/training_data_100m_20191105_V4_no_time_gaps_africa_subset.csv"), "NIRv_8d_Int16")
plot(GetMatrixFromSF(NIRv_8d_16int)[1,]~GetDates8d(), type="l")

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
MODDetectBreaks = function(InputTS, scrange=c(2009, 2019), scsig=0.05, breaks="LWZ", plot=FALSE)
{
    # The input should be a single row of a matrix.
    InputTS = GetTS(InputTS)
    Observations = sum(!is.na(InputTS))
    print(paste("Observations for point:", Observations))
    h = frequency(InputTS)
    #plot(Point1TS~MyDates, type="l"); abline(v=as.POSIXct(as.Date("2017-01-01")))
    if (Observations > h*2) {
        bpp = bfastpp(InputTS, order=3)
        if (!is.null(scrange)) {
            SC = sctest(efp(response ~ (harmon + trend),
                            data=bpp[bpp$time > scrange[1] & bpp$time < scrange[2],],
                            h=h/Observations, type="OLS-MOSUM"))$p.value > scsig
            if (SC) return(FALSE)
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
        print(bpp[bpOptim$breakpoints, "time"])
        return(bpp[bpOptim$breakpoints, "time"])
    } else {
        print("too cloudy")
        return(NA)
    }
}

MODDetectBreaks(GetMatrixFromSF(EVI_8d_16int)[EVI_8d_16int$sample_id == unique(ChangeIDs_UTM)[1],], breaks="LWZ", plot=TRUE)
MODDetectBreaks(GetMatrixFromSF(NDMI_8d_16int)[NDMI_8d_16int$sample_id == unique(ChangeIDs_UTM)[1],], breaks="LWZ", plot=TRUE)
MyBreak = MODDetectBreaks(GetMatrixFromSF(NIRv_8d_16int)[NIRv_8d_16int$sample_id == unique(ChangeIDs_UTM)[2],], breaks="BIC", plot=TRUE)
MyReference = as.data.frame(RawData)[RawData$sample_id==unique(ChangeIDs_UTM)[2] & RawData$change_at_300m=="yes","year_fraction"]
abline(v=MyReference, col="blue")

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
IsBreakInTargetYear = function(BreakTimes, TargetYear, threshold=0.25)
{
    #if (is.na(TargetYear))
    #    TargetYear = 2016 # If there is no break, we look at whether we predicted a break in 2016
    # TODO: Needs to be updated; previously, lack of break meant that the break time would be set to NA,
    # but now it's no longer the case, it's change_at_300m = FALSE and reference_year=<year>
    return(any(BreakTimes > TargetYear - threshold & BreakTimes < TargetYear+1+threshold))
}

# 4b) Vectorised version (takes a list of MODDetectBreaks() output and a column of target years)
# Returns a column of whether BFSAT predicted the break at that time or not.
VectorisedIsBreakInTargetYear = function(BreakList, threshold=0.25, TY=TargetYears)
{
    i = 1:length(BreakList)
    return(sapply(i, function(i){return(IsBreakInTargetYear(BreakList[[i]], TY[i], threshold=threshold))}))
}

IsBreakInTargetYear(MyBreak, MyReference)

# 4c) Wrapper for running 3+4 in one.
# Result is a logical vector saying how many times we accurately predicted a break,
# and how many times not. One value per year.
TestMODBreakpointDetection = function(i, threshold=0.25, TargetYears=TargetYears, ...)
{
    BreakTimes = MODDetectBreaks(i, ...)
    TargetYear = TargetYears[i]
    if (is.na(BreakTimes[1])) return(NA)
    # If MODDetectBreaks returned FALSE, we also return FALSE: no breaks were detected
    if (!BreakTimes[1]) return(FALSE)
    return(IsBreakInTargetYear(BreakTimes, TargetYear, threshold))
}

# 5) Get statistics. 
# TODO: Update with the current CSV scheme.
# TODO: Should take two inputs: BFAST predictions and truth.
FPStats = function(IsBreakInTargetYear, uncertain=TRUE)
{
    # ChangePoints is TRUE if there was a change that year, FALSE if no change.
    ChangePoints = AllPoints$ChangeType != "no LC change"
    if (!uncertain) # We may filter out if it's uncertain (NAs are ignored)
        ChangePoints[AllPoints$confidence == "-1"] = NA
    # We predicted a break (IsBreakInTargetYear == TRUE) and it was a break (ChangePoints == TRUE)
    TruePositiveCount = sum(IsBreakInTargetYear[ChangePoints], na.rm=TRUE)
    # We predicted a break but htere wasn't one
    FalsePositiveCount = sum(IsBreakInTargetYear[!ChangePoints], na.rm=TRUE)
    # We predicted no break, and there were none
    TrueNegativeCount = sum(!IsBreakInTargetYear[!ChangePoints], na.rm=TRUE)
    # We predicted no break, but there was one
    FalseNegativeCount = sum(!IsBreakInTargetYear[ChangePoints], na.rm=TRUE)
    # Percent of true positives out of all change
    TruePositiveRate = TruePositiveCount / sum(ChangePoints, na.rm=TRUE) # AKA Sensitivity
    Specificity = TrueNegativeCount / (TrueNegativeCount + FalsePositiveCount)
    # Percent of false positive out of no change
    FalsePositiveRate = FalsePositiveCount / sum(!ChangePoints, na.rm=TRUE) # False positive rate or alpha or p-value or Type I Error
    SignalToNoise = TruePositiveCount / FalsePositiveCount
    SignalToNoiseRate = TruePositiveRate / FalsePositiveRate # Likelihood Ratio for Positive Tests
    PositivePredictiveValue = TruePositiveCount / (TruePositiveCount + FalsePositiveCount)
    Accuracy = (TruePositiveCount + TrueNegativeCount) / length(ChangePoints)
    return(data.frame(TruePositiveCount, FalsePositiveCount, TrueNegativeCount, FalseNegativeCount,
                      TruePositiveRate, Specificity, FalsePositiveRate, SignalToNoise, SignalToNoiseRate,
                      PositivePredictiveValue, Accuracy))
}

# TODO: test different VIs. Change input to LoadVITS, select which one works best (first under defaults).

# TODO: make an optimisation routine where we run BFAST with different parameters,
# get statistics and choose the best values based on the statistics. Maybe using optim().
# Caveat that if we do it automatically, we may run into edge cases,
# i.e. all NA and one correct break + one correct no break is best.
