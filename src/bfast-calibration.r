# Script for determining which BFAST/change detection method works the best,
# according to change validation data.

library(sf)

# Import the validation CSV; it contains one entry per year that gives the approximate time of the year.

# 1) Load the CSV. Ideally we only need an st_read() to get an sf object
LoadReferenceData = function()
{
    
}

# 2) Extract time series data from the coordinates,
# and cache it in a CSV/GPKG so that we don't need to do that again.
# This is where we select different VIs.
# TODO: handle UTM zones. The input files are VRTs in ~6 projection systems.
#       We need to try and figure out which UTM zone the point is at, and then extract it.
#       So as an input, we can use a directory that contains all these VRTs;
#       read in each, and get from the reference data the points that fall within that image.
#       (Making sure to handle cases where we hae a point in multiple files!)
# The output is a data.frame, with rows being unique points, and columns being timesteps;
# data is the value of the vegetation index.
# TODO: This means we also need to get the dates somehow, see 3). We can get it from a ts() with
#       as.Date(date_decimal(as.numeric(time(ExampleTS))))
# TODO: Can simplify the input to just a directory path and the reference points,
#       and save the cached output with something derived from the name of the directory path.
# TODO: Check utils/cachevi.R
LoadVITS = function(vi="EVI", pointlocs=AllChanged, prefix="Changed", sourcedir="/data/mep_cg1/MOD_S10/additional_VIs_run04")
{
    VITSFile = paste0("../data/", prefix, vi, "TS.csv")
    if (!file.exists(VITSFile))
    {
        VIMosaicFile = paste0("/data/users/Public/greatemerald/modis/raw-input-mosaic-", prefix, "-", vi, ".vrt")
        if (!file.exists(VIMosaicFile))
        {
            InputFiles = list.files(sourcedir, glob2rx(paste0("MOD_S10_TOC_*_250m_C6_", vi, ".tif")), recursive = TRUE, full.names = TRUE)
            gdalbuildvrt(InputFiles, VIMosaicFile)
        }
        VIMosaic = brick(VIMosaicFile)
        ChangedVITS = extract(VIMosaic, pointlocs)
        write.csv(ChangedVITS, VITSFile)
    } else ChangedVITS = as.matrix(read.csv(VITSFile, row.names=1))
    return(ChangedVITS)
}

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
MODDetectBreaks = function(i, ChangedVITS, h = 36, start=as.Date("2009-01-01"),
                           scrange=c(2009, 2019), scsig=0.05)
{
    MyDates = dates[dates >= start]
    Point1TS = ChangedVITS[i, dates >= start]
    Observations = sum(!is.na(Point1TS))
    print(paste("Observations for point:", Observations))
    #plot(Point1TS~MyDates, type="l"); abline(v=as.POSIXct(as.Date("2017-01-01")))
    if (Observations > h*2 && Observations > 42) {
        P1TS = bfastts(Point1TS, MyDates, "10-day")
        P1PP = bfastpp(P1TS, order=3)
        if (!is.null(scrange)) {
            SC = sctest(efp(response ~ (harmon + trend),
                            data=P1PP[P1PP$time > scrange[1] & P1PP$time < scrange[2],],
                            h=h/Observations, type="OLS-MOSUM"))$p.value > scsig
            if (SC) return(FALSE)
        }
        
        P1BP = try(breakpoints(formula(response~trend+harmon), data=P1PP, h=h))
        #if (length(P1BP$breakpoints) > 0) {
        #    abline(v=date_decimal(P1PP[P1BP$breakpoints, "time"]), col="red")
        #}
        #print(P1BP)
        if ("try-error" %in% class(P1BP)) {
            print("An error has occurred:")
            print(P1BP)
            return(NA)
        }
        if (all(is.na(P1BP$breakpoints)))
            return(FALSE)
        print(P1PP[P1BP$breakpoints, "time"])
        return(P1PP[P1BP$breakpoints, "time"])
    } else {
        print("too cloudy")
        return(NA)
    }
}

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
    if (is.na(TargetYear))
        TargetYear = 2016 # If there is no break, we look at whether we predicted a break in 2016
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
