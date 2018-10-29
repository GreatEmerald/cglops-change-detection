# Tuning BFAST params to get a better true to false change detection ratio
library(raster)
library(sf)
library(lubridate)
library(strucchange)
library(bfast)

LoadTrueChange = function()
{
    breakval = st_read("../data/valgroup27_exports_20180904_LandCoverChangeDetection_1.csv", options=c("X_POSSIBLE_NAMES=sample_x", "Y_POSSIBLE_NAMES=sample_y"))
    st_crs(breakval) = 4326
    ChangedPixels = breakval[breakval$legend_name=="CGLOPS Land Cover Change Detection validation - Did the LC type at 100m x 100m change?",]
    FilteredPoints = breakval[breakval$legend_name=="CGLOPS Land Cover Change Detection validation - Mark LC change",]
    FilteredPoints = FilteredPoints[!duplicated(FilteredPoints[,c("name", "sample_x", "sample_y")]),]
    ChangedOnly = ChangedPixels[ChangedPixels$name!="No",]
    ChangedOnly$geometry %in% FilteredPoints$geometry # There are two points that are labelled as change but no subpixels
    ExtraPoints = ChangedOnly[!ChangedOnly$geometry %in% FilteredPoints$geometry,] # They are complete change in 2016, so relabel as such
    ExtraPoints$name = factor("LC change 2016", levels = levels(ChangedOnly$name))
    
    AllChanged = rbind(FilteredPoints, ExtraPoints)
    return(AllChanged)
}

# Make a DF that has 1010 points, with a column describing what year there was a change, or NA
LoadAllReference = function()
{
    breakval = st_read("../data/valgroup27_exports_20180904_LandCoverChangeDetection_1.csv", options=c("X_POSSIBLE_NAMES=sample_x", "Y_POSSIBLE_NAMES=sample_y"))
    st_crs(breakval) = 4326
    # Filter out all points to get all unique ones
    UniquePoints = breakval[breakval$legend_name == "CGLOPS Land Cover Change Detection validation - Did the LC type at 100m x 100m change?",]
    # Filter out to only change
    FilteredPoints = breakval[breakval$legend_name=="CGLOPS Land Cover Change Detection validation - Mark LC change",]
    FilteredPoints = FilteredPoints[!duplicated(FilteredPoints[,c("name", "sample_x", "sample_y")]),]
    # Filter out to the ones that actually did change
    ChangedOnly = UniquePoints[UniquePoints$name!="No",]
    # Add extra points that are full pixel change
    ExtraPoints = ChangedOnly[!ChangedOnly$geometry %in% FilteredPoints$geometry,]
    ExtraPoints$name = factor("LC change 2016", levels = levels(ChangedOnly$name))
    AllChanged = rbind(FilteredPoints, ExtraPoints)
    # Combine all points with change with points without change again
    AllPoints = UniquePoints[!UniquePoints$validation_id %in% AllChanged$validation_id,]
    AllPoints = rbind(AllPoints, AllChanged)
    # Add a "Change" column
    AllPoints$Change = NA
    AllPoints[AllPoints$name == "LC change 2016", "Change"] = 2016
    AllPoints[AllPoints$name == "LC change 2017", "Change"] = 2017
    return(AllPoints)
}

LoadVITS = function(vi="EVI", pointlocs=AllChanged, prefix="Changed")
{
    VITSFile = paste0("../data/", prefix, vi, "TS.csv")
    if (!file.exists(VITSFile))
    {
        VIMosaicFile = paste0("/data/users/Public/greatemerald/modis/raw-input-mosaic-", vi, ".vrt")
        if (!file.exists(VIMosaicFile))
        {
            InputFiles = c(paste0("/data/mep_cg1/MOD_S10/additional_VIs_new/X16Y06/MOD_S10_TOC_X16Y06_20090101-20171231_250m_C6_", vi, ".tif"),
                           paste0("/data/mep_cg1/MOD_S10/additional_VIs_new/X17Y06/MOD_S10_TOC_X17Y06_20090101-20171231_250m_C6_", vi, ".tif"),
                           paste0("/data/mep_cg1/MOD_S10/additional_VIs_new/X18Y06/MOD_S10_TOC_X18Y06_20090101-20171231_250m_C6_", vi, ".tif"),
                           paste0("/data/mep_cg1/MOD_S10/additional_VIs_new/X19Y06/MOD_S10_TOC_X19Y06_20090101-20171231_250m_C6_", vi, ".tif"))
            gdalbuildvrt(InputFiles, VIMosaicFile)
        }
        VIMosaic = brick(VIMosaicFile)
        ChangedVITS = extract(VIMosaic, pointlocs)
        write.csv(ChangedVITS, VITSFile)
    } else ChangedVITS = as.matrix(read.csv(VITSFile, row.names=1))
    return(ChangedVITS)
}

AC.df = LoadTrueChange()
class(AC.df) = "data.frame"
dates = GetDatesFromDir("/data/mep_cg1/MOD_S10/")

TargetYears = integer(nrow(AC.df))
TargetYears[as.character(AC.df[,"name"])=="LC change 2016"] = 2016
TargetYears[as.character(AC.df[,"name"])=="LC change 2017"] = 2017

ChangedEVITS = LoadVITS("EVI")
ChangedNDMITS = LoadVITS("NDMI")

TestMODBreakpointDetection = function(i, ChangedVITS, h = 36, threshold=0.25,
                                      start=as.Date("2009-01-01"), scrange=c(2009, 2019), scsig=0.05,
                                      TargetYears=TargetYears)
{
    MyDates = dates[dates >= start]
    Point1TS = ChangedVITS[i, dates >= start]
    TargetYear = TargetYears[i]
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
        if ("try-error" %in% class(P1BP))
            return(NA)
        if (all(is.na(P1BP$breakpoints)))
            return(FALSE)
        print(P1PP[P1BP$breakpoints, "time"])
        return(any(P1PP[P1BP$breakpoints, "time"] > TargetYear - threshold &
                       P1PP[P1BP$breakpoints, "time"] < TargetYear+1+threshold))
    } else {
        print("too cloudy")
        return(NA)
    }
}

M_EVI_2009 = sapply(1:nrow(AC.df), TestMODBreakpointDetection, ChangedEVITS, scrange=NULL)
mean(M_EVI_2009) # 53%
M_NDMI_2009 = sapply(1:nrow(AC.df), TestMODBreakpointDetection, ChangedNDMITS, scrange=NULL)
mean(M_NDMI_2009) # 44%
mean(M_EVI_2009 | M_NDMI_2009) # 64%

# From 2014
M_EVI_2014 = sapply(1:nrow(AC.df), TestMODBreakpointDetection, ChangedEVITS, start=as.Date("2014-01-01"), scrange=NULL)
mean(M_EVI_2014, na.rm=TRUE) # 53%
M_NDMI_2014 = sapply(1:nrow(AC.df), TestMODBreakpointDetection, ChangedNDMITS, start=as.Date("2014-01-01"), scrange=NULL)
mean(M_NDMI_2014, na.rm=TRUE) # 50%
mean(M_EVI_2014 | M_NDMI_2014, na.rm=TRUE) # 66%
# This is better and faster...

## With sctest
M_EVI_2009_SC = sapply(1:nrow(AC.df), TestMODBreakpointDetection, ChangedEVITS)
mean(M_EVI_2009_SC) # 47%
M_NDMI_2009_SC = sapply(1:nrow(AC.df), TestMODBreakpointDetection, ChangedNDMITS)
mean(M_NDMI_2009_SC) # 41%
mean(M_EVI_2009_SC | M_NDMI_2009_SC) # 61%
# Only a few discards

# From 2014
M_EVI_2014_SC = sapply(1:nrow(AC.df), TestMODBreakpointDetection, ChangedEVITS, start=as.Date("2014-01-01"))
mean(M_EVI_2014_SC, na.rm=TRUE) # 38%, a lot of discards
M_NDMI_2014_SC = sapply(1:nrow(AC.df), TestMODBreakpointDetection, ChangedNDMITS, start=as.Date("2014-01-01"))
mean(M_NDMI_2014_SC, na.rm=TRUE) # 35%
mean(M_EVI_2014_SC | M_NDMI_2014_SC, na.rm=TRUE) # 50% Nicely complimentary

# Make SC significance loose
M_EVI_2009_SC_low = sapply(1:nrow(AC.df), TestMODBreakpointDetection, ChangedEVITS, scsig=0.1)
mean(M_EVI_2009_SC_low) # 53% - no discards
M_NDMI_2009_SC_low = sapply(1:nrow(AC.df), TestMODBreakpointDetection, ChangedNDMITS, scsig=0.1)
mean(M_NDMI_2009_SC_low) # 44% - ditto
mean(M_EVI_2009_SC_low | M_NDMI_2009_SC_low) # 64%
# No discards

# From 2014
M_EVI_2014_SC_low = sapply(1:nrow(AC.df), TestMODBreakpointDetection, ChangedEVITS, start=as.Date("2014-01-01"), scsig=0.1)
mean(M_EVI_2014_SC_low, na.rm=TRUE) # 38%, same
M_NDMI_2014_SC_low = sapply(1:nrow(AC.df), TestMODBreakpointDetection, ChangedNDMITS, start=as.Date("2014-01-01"), scsig=0.1)
mean(M_NDMI_2014_SC_low, na.rm=TRUE) # 43%, a bit better
mean(M_EVI_2014_SC_low | M_NDMI_2014_SC_low, na.rm=TRUE) # 50%, same

# Targetted sctest
M_EVI_2014_SC_2 = sapply(1:nrow(AC.df), TestMODBreakpointDetection, ChangedEVITS, start=as.Date("2014-01-01"), scrange=c(2015.75, 2018.25))
mean(M_EVI_2014_SC_2, na.rm=TRUE) # 19%, a lot of discards
M_NDMI_2014_SC_2 = sapply(1:nrow(AC.df), TestMODBreakpointDetection, ChangedNDMITS, start=as.Date("2014-01-01"), scrange=c(2015.75, 2018.25))
mean(M_NDMI_2014_SC_2, na.rm=TRUE) # 16%, even more discards
mean(M_EVI_2014_SC_2 | M_NDMI_2014_SC_2, na.rm=TRUE) # 28%

# Loose targetted sctest
M_EVI_2014_SC_2l = sapply(1:nrow(AC.df), TestMODBreakpointDetection, ChangedEVITS, start=as.Date("2014-01-01"), scrange=c(2015.75, 2018.25), scsig=0.1)
mean(M_EVI_2014_SC_2l, na.rm=TRUE) # 25%, better
M_NDMI_2014_SC_2l = sapply(1:nrow(AC.df), TestMODBreakpointDetection, ChangedNDMITS, start=as.Date("2014-01-01"), scrange=c(2015.75, 2018.25), scsig=0.1)
mean(M_NDMI_2014_SC_2l, na.rm=TRUE) # 31%, much better
mean(M_EVI_2014_SC_2l | M_NDMI_2014_SC_2l, na.rm=TRUE) # 44%, pretty tolerable

###
### Try to process points of no change to see whether the false detection rate can be improved
###



# Get all data
AllPoints = LoadAllReference()

# Extract TS out of all those
AllEVITS = LoadVITS(pointlocs=AllPoints, prefix="All")
AllNDMITS = LoadVITS("NDMI", pointlocs=AllPoints, prefix="All")

# Define target year
AllTargetYears = AllPoints$Change
AllTargetYears[is.na(AllTargetYears)] = 2016

# Test the true vs false positive rates
FPStats = function(IsBreakInTargetYear)
{
    ChangePoints = !is.na(AllPoints$Change)
    TruePositiveCount = sum(IsBreakInTargetYear[ChangePoints], na.rm=TRUE)
    FalsePositiveCount = sum(IsBreakInTargetYear[!ChangePoints], na.rm=TRUE)
    TruePositiveRate = TruePositiveCount / sum(ChangePoints)
    FalsePositiveRate = FalsePositiveCount / sum(!ChangePoints)
    SignalToNoise = TruePositiveCount / FalsePositiveCount
    SignalToNoiseRate = TruePositiveRate / FalsePositiveRate
    return(data.frame(TruePositiveCount, FalsePositiveCount, TruePositiveRate, FalsePositiveRate,
                      SignalToNoise, SignalToNoiseRate))
}

# First: reference values
A_EVI_2009 = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, AllEVITS, scrange=NULL, TargetYears=AllTargetYears)
FPStats(A_EVI_2009) # 3.4%
A_NDMI_2009 = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, AllNDMITS, scrange=NULL, TargetYears=AllTargetYears)
FPStats(A_NDMI_2009) # 4.9%
FPStats(A_EVI_2009 | A_NDMI_2009) # 3.5%

# If we add the sctest
A_EVI_2009_SC = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, AllEVITS, TargetYears=AllTargetYears)
FPStats(A_EVI_2009_SC) # 3.6%
A_NDMI_2009_SC = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, AllNDMITS, TargetYears=AllTargetYears)
FPStats(A_NDMI_2009_SC) # 5.3%
FPStats(A_EVI_2009_SC | A_NDMI_2009_SC) # 3.7%

# Targetted sctest
A_EVI_2009_SCT = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, AllEVITS, TargetYears=AllTargetYears, scrange=c(2015.75, 2018.25))
FPStats(A_EVI_2009_SCT) # 3.2%
A_NDMI_2009_SCT = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, AllNDMITS, TargetYears=AllTargetYears, scrange=c(2015.75, 2018.25))
FPStats(A_NDMI_2009_SCT) # 5.0%
FPStats(A_EVI_2009_SCT | A_NDMI_2009_SCT) # 4.1%

# From 2014
A_EVI_2014 = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, AllEVITS, TargetYears=AllTargetYears, start=as.Date("2014-01-01"), scrange=NULL)
FPStats(A_EVI_2014) # 2.9%
A_NDMI_2014 = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, AllNDMITS, TargetYears=AllTargetYears, start=as.Date("2014-01-01"), scrange=NULL)
FPStats(A_NDMI_2014) # 2.9%
FPStats(A_EVI_2014 | A_NDMI_2014) # 3.0%

# Targetted sctest from 2014
A_EVI_2014_SCT = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, AllEVITS, TargetYears=AllTargetYears, start=as.Date("2014-01-01"), scrange=c(2015.75, 2018.25))
FPStats(A_EVI_2014_SCT) # 2.4%
A_NDMI_2014_SCT = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, AllNDMITS, TargetYears=AllTargetYears, start=as.Date("2014-01-01"), scrange=c(2015.75, 2018.25))
FPStats(A_NDMI_2014_SCT) # 2.4%
FPStats(A_EVI_2014_SCT | A_NDMI_2014_SCT) # 2.5%

# Loose targetted sctest from 2014
A_EVI_2014_SCTL = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, AllEVITS, TargetYears=AllTargetYears, start=as.Date("2014-01-01"), scrange=c(2015.75, 2018.25), scsig=0.1)
FPStats(A_EVI_2014_SCTL) # 2.5%
A_NDMI_2014_SCTL = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, AllNDMITS, TargetYears=AllTargetYears, start=as.Date("2014-01-01"), scrange=c(2015.75, 2018.25), scsig=0.1)
FPStats(A_NDMI_2014_SCTL) # 3.7%
FPStats(A_EVI_2014_SCTL | A_NDMI_2014_SCTL) # 3.2%

## Overall this does not seem to help. Could other methods work better?

# Try Eline's t-test
TestMODttest = function(i, ChangedVITS, TargetYears=AllTargetYears, sig=0.05)
{
    Point1TS = ChangedVITS[i, ]
    TargetYear = TargetYears[i]
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

# Monkey method
FPStats(rbinom(nrow(AllPoints), 1, 0.5)) # ~3.4%...

# t-test
T_EVI_2009 = sapply(1:nrow(AllPoints), TestMODttest, AllEVITS)
FPStats(T_EVI_2009) # 3.7%
T_NDMI_2009 = sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS)
FPStats(T_NDMI_2009) # 10.6%
FPStats(T_EVI_2009 | T_NDMI_2009) # 7.6%

# And if we change the significance: tighter
T_EVI_2016T = sapply(1:nrow(AllPoints), TestMODttest, AllEVITS, sig=0.01)
FPStats(T_EVI_2016T) # 2.3%
T_NDMI_2016T = sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS, sig=0.01)
FPStats(T_NDMI_2016T) # 11.5%
FPStats(T_EVI_2016T | T_NDMI_2016T) # 5.5%

# Looser
T_EVI_2016L = sapply(1:nrow(AllPoints), TestMODttest, AllEVITS, sig=0.1)
FPStats(T_EVI_2016L) # 5.2%
T_NDMI_2016L = sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS, sig=0.1)
FPStats(T_NDMI_2016L) # 9.2%
FPStats(T_EVI_2016L | T_NDMI_2016L) # 7.3%

# Check around it
FPStats(sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS, sig=0.07)) # 10.5% - still one of the best, with 23% true positives
FPStats(sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS, sig=0.13)) # 7.9%
FPStats(sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS, sig=0.08)) # 9.6%
FPStats(sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS, sig=0.06)) # 9.8%
FPStats(sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS, sig=0.02)) # 10.2%
FPStats(sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS, sig=0.005)) # 9.5%

# How about BFAST Monitor?

TestMODMonitor = function(i, ChangedVITS, TargetYears=AllTargetYears, threshold=0.25, cloud_threshold=42, ...)
{
    MyDates = dates
    Point1TS = ChangedVITS[i, ]
    TargetYear = TargetYears[i]
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

# Defaults
M_EVI_2009 = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS)
FPStats(M_EVI_2009) # 3.8%
M_NDMI_2009 = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS)
FPStats(M_NDMI_2009) # 3.3%
FPStats(M_EVI_2009 | M_NDMI_2009) # 3.6%

# Adjust significance
M_EVI_2009T = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS, level=0.001)
FPStats(M_EVI_2009T) # 3.1%
M_NDMI_2009T = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS, level=0.001)
FPStats(M_NDMI_2009T) # 4.3%
FPStats(M_EVI_2009T | M_NDMI_2009T) # 3.9%

# Adjust type
M_EVI_2009TC = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS, level=0.001, type="OLS-CUSUM")
FPStats(M_EVI_2009TC) # 3.6%
M_NDMI_2009TC = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS, level=0.001, type="OLS-CUSUM")
FPStats(M_NDMI_2009TC) # 4.5%
FPStats(M_EVI_2009TC | M_NDMI_2009TC) # 3.7%

# Adjust h
M_EVI_2009TH5 = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS, level=0.001, h=0.5)
FPStats(M_EVI_2009TH5) # 3.8%
M_NDMI_2009TH5 = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS, level=0.001, h=0.5)
FPStats(M_NDMI_2009TH5) # 4.3%
FPStats(M_EVI_2009TH5 | M_NDMI_2009TH5) # 4.2%

M_EVI_2009TH1 = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS, level=0.001, h=1)
FPStats(M_EVI_2009TH1) # 2.7%
M_NDMI_2009TH1 = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS, level=0.001, h=1)
FPStats(M_NDMI_2009TH1) # 7.2% - pretty good but very few detections overall!
FPStats(M_EVI_2009TH1 | M_NDMI_2009TH1) # 4.6%

# Add lag
M_EVI_2009TL = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS, level=0.001, lag=1, cloud_threshold=75)
FPStats(M_EVI_2009TL) # 5.1%
M_NDMI_2009TL = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS, level=0.001, lag=1, cloud_threshold=75)
FPStats(M_NDMI_2009TL) # 4.3%
FPStats(M_EVI_2009TL | M_NDMI_2009TL) # 4.8%

# Add slag
M_EVI_2009TS = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS, level=0.001, slag=1, cloud_threshold=75)
FPStats(M_EVI_2009TS) # 3.7%
M_NDMI_2009TS = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS, level=0.001, slag=1, cloud_threshold=75)
FPStats(M_NDMI_2009TS) # Error
FPStats(M_EVI_2009TS | M_NDMI_2009TS) # Error

# Trend only
M_EVI_2009TO = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS, level=0.001, formula=response~trend)
FPStats(M_EVI_2009TO) # 3.8%
M_NDMI_2009TO = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS, level=0.001, formula=response~trend)
FPStats(M_NDMI_2009TO) # 4.0%
FPStats(M_EVI_2009TO | M_NDMI_2009TO) # 3.5%

# Harmonics only
M_EVI_2009HO = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS, level=0.001, formula=response~harmon)
FPStats(M_EVI_2009HO) # 4.8%
M_NDMI_2009HO = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS, level=0.001, formula=response~trend)
FPStats(M_NDMI_2009HO) # 4.0%
FPStats(M_EVI_2009HO | M_NDMI_2009HO) # 3.8%

# All in all, t-test wins... But signal to noise is still only 10%
