# Tuning BFAST params to get a better true to false change detection ratio
library(raster)
library(sf)
library(lubridate)
library(strucchange)
library(bfast)
library(gdalUtils)

LoadTrueChange_2018summer = function()
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
LoadAllReference_2018summer = function()
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

LoadAllReference = function()
{
    ValPoints = st_read("../data/ValidationPoints-December.csv", options=c("X_POSSIBLE_NAMES=sample_x", "Y_POSSIBLE_NAMES=sample_y"), stringsAsFactors = FALSE)
    st_crs(ValPoints) = 4326
    ValPoints$ChangeYear = as.integer(ValPoints$ChangeYear)
    ValPoints$ChangeType = factor(ValPoints$ChangeType, levels=c("LC change", "Fire dynamics", "Water dynamics", "LC dynamics", "no LC change"))
    return(ValPoints)
}

LoadTrueChange = function()
{
    ValPoints = LoadAllReference()
    return(ValPoints[ValPoints$ChangeType != "no LC change",])
}

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

#### Below is if looking at true change only, not very useful
AC.df = LoadTrueChange()
class(AC.df) = "data.frame"

TargetYears = AC.df$ChangeYear #integer(nrow(AC.df))
#TargetYears[as.character(AC.df[,"name"])=="LC change 2016"] = 2016
#TargetYears[as.character(AC.df[,"name"])=="LC change 2017"] = 2017
# Assume that the change year of unmarked change is 2016
TargetYears[AC.df$ChangeType != "no LC change" & is.na(TargetYears)] = 2016 

#ChangedEVITS = LoadVITS("EVI")
#ChangedNDMITS = LoadVITS("NDMI")

# Function that returns whether there has been a break detected in the target year
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

IsBreakInTargetYear = function(BreakTimes, TargetYear, threshold=0.25)
{
    if (is.na(TargetYear))
        TargetYear = 2016 # If there is no break, we look at whether we predicted a break in 2016
    return(any(BreakTimes > TargetYear - threshold & BreakTimes < TargetYear+1+threshold))
}

VectorisedIsBreakInTargetYear = function(BreakList, threshold=0.25, TY=TargetYears)
{
    i = 1:length(BreakList)
    return(sapply(i, function(i){return(IsBreakInTargetYear(BreakList[[i]], TY[i], threshold=threshold))}))
}

TestMODBreakpointDetection = function(i, threshold=0.25, TargetYears=TargetYears, ...)
{
    BreakTimes = MODDetectBreaks(i, ...)
    TargetYear = TargetYears[i]
    if (is.na(BreakTimes[1])) return(NA)
    if (!BreakTimes[1]) return(FALSE)
    return(IsBreakInTargetYear(BreakTimes, TargetYear, threshold))
}

###
### Try to process points of no change to see whether the false detection rate can be improved
###



# Get all data
AllPoints = LoadAllReference()

# Extract TS out of all those
AllEVITS = LoadVITS("EVI", pointlocs=AllPoints, prefix="All")
AllNDMITS = LoadVITS("NDMI", pointlocs=AllPoints, prefix="All")
AllNIRVTS = LoadVITS("NIRV", pointlocs=AllPoints, prefix="All")

dates = GetDatesFromDir("/data/mep_cg1/MOD_S10/MOD_S10/")[1:ncol(AllEVITS)]

# Define target year
AllTargetYears = AllPoints$ChangeYear
AllTargetYears[AllPoints$ChangeType != "no LC change" & is.na(AllTargetYears)] = 2016 
#AllTargetYears[is.na(AllTargetYears)] = 2016

# Test the true vs false positive rates
FPStats = function(IsBreakInTargetYear, uncertain=TRUE)
{
    ChangePoints = AllPoints$ChangeType != "no LC change"
    if (!uncertain)
        ChangePoints[AllPoints$confidence == "-1"] = NA
    TruePositiveCount = sum(IsBreakInTargetYear[ChangePoints], na.rm=TRUE)
    FalsePositiveCount = sum(IsBreakInTargetYear[!ChangePoints], na.rm=TRUE)
    TrueNegativeCount = sum(!IsBreakInTargetYear[!ChangePoints], na.rm=TRUE)
    FalseNegativeCount = sum(!IsBreakInTargetYear[ChangePoints], na.rm=TRUE)
    TruePositiveRate = TruePositiveCount / sum(ChangePoints, na.rm=TRUE) # Sensitivity
    Specificity = TrueNegativeCount / (TrueNegativeCount + FalsePositiveCount)
    FalsePositiveRate = FalsePositiveCount / sum(!ChangePoints, na.rm=TRUE) # False positive rate or alpha or p-value or Type I Error
    SignalToNoise = TruePositiveCount / FalsePositiveCount
    SignalToNoiseRate = TruePositiveRate / FalsePositiveRate # Likelihood Ratio for Positive Tests
    PositivePredictiveValue = TruePositiveCount / (TruePositiveCount + FalsePositiveCount)
    Accuracy = (TruePositiveCount + TrueNegativeCount) / length(ChangePoints)
    return(data.frame(TruePositiveCount, FalsePositiveCount, TrueNegativeCount, FalseNegativeCount,
                      TruePositiveRate, Specificity, FalsePositiveRate, SignalToNoise, SignalToNoiseRate,
                      PositivePredictiveValue, Accuracy))
}

# Function for getting a graph for visualising how many points each algorithm identifies
GroupPlots = function(IsBreakInTargetYear, uncertain=TRUE, ...)
{
    #IsBreakInTargetYear[is.na(IsBreakInTargetYear)] = FALSE
    order = c(2, 1, 5, 3, 4)
    ylim = c(0, 250)
    ChangeTypes = levels(as.factor(AllPoints$ChangeType))
    Reference = AllPoints[!is.na(IsBreakInTargetYear), ]
    if (!uncertain)
        Reference = Reference[Reference$confidence != "-1",]
    
    MyTable = function(type, Tested)
    {
        sum(Tested[["ChangeType"]] == type, na.rm=TRUE)
    }
    
    totals = sapply(ChangeTypes, MyTable, Tested=Reference)
    totals_c = as.integer(totals)
    names(totals_c) = names(totals)
    Detected = AllPoints[IsBreakInTargetYear, ]
    if (!uncertain)
        Detected = Detected[Detected$confidence != "-1",]
    detected_c = sapply(ChangeTypes, MyTable, Tested=Detected)
    names(detected_c) = ChangeTypes
    undetected = totals_c - detected_c
    stopifnot(all(undetected >= 0))
    barmatrix = rbind(detected_c, undetected)
    barplot(barmatrix[,order], ...)
}

# First: reference values
A_EVI_2009_dates = lapply(1:nrow(AllPoints), MODDetectBreaks, AllEVITS, scrange=NULL)
A_EVI_2009B = VectorisedIsBreakInTargetYear(A_EVI_2009_dates)
A_EVI_2009 = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, AllEVITS, scrange=NULL, TargetYears=AllTargetYears)
FPStats(A_EVI_2009)
FPStats(A_EVI_2009B) # 61% OA
GroupPlots(A_EVI_2009B, main="EVI without sctest, 2009-2017")
options(strucchange.use_armadillo=TRUE)
A_EVI_2009C = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, ChangedVITS=AllEVITS, scrange=NULL, TargetYears=AllTargetYears)
all.equal(A_EVI_2009C, A_EVI_2009B) #identical
A_NDMI_2009 = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, ChangedVITS=AllNDMITS, scrange=NULL, TargetYears=AllTargetYears)
GroupPlots(A_NDMI_2009, main="NDMI without sctest, 2009-2017")
A_NDMI_2009_raw = lapply(1:nrow(AllPoints), MODDetectBreaks, ChangedVITS=AllNDMITS, scrange=NULL)
FPStats(A_NDMI_2009) # 66%
FPStats(A_EVI_2009 | A_NDMI_2009) # 67%

### Visualisation ###

A_EVI_2009_SC_raw = sapply(1:nrow(AllPoints), MODDetectBreaks, ChangedVITS=AllEVITS)

which(!A_NDMI_2009 & AllPoints$ChangeType == "LC change") # These do not detect change 
A_NDMI_2009_raw[which(!A_NDMI_2009 & AllPoints$ChangeType == "LC change")]
AllPoints[which(!A_NDMI_2009 & AllPoints$ChangeType == "LC change"),]

FalseDetections = A_EVI_2009_SC & AllPoints$ChangeType == "no LC change"
which(FalseDetections) # False detections
AllPoints[FalseDetections,]
# 1 7, ~9, 11, 13, 22, ~27, 30, 31, 34, 46, 50, 51, 57, 65, 72, 74, 76, 87
# Drought: [1] 9 50
# New cultivation: 7
# Funky: 11 [30] 13 22 34 46 [72]
# Crop variability: 31 51 [65]
# Clear change: [57] [76] [87]
TryIdx = 72
BreakDates = as.Date(lubridate::date_decimal(A_EVI_2009_SC_raw[[which(FalseDetections)[TryIdx]]]))
plot(AllEVITS[which(FalseDetections)[TryIdx],]~dates, type="l", main=AllPoints[which(FalseDetections)[TryIdx],][["comment"]], ylab="EVI")
for (i in 1:length(BreakDates)) abline(v=BreakDates[i], col="red")

###

# If we add the sctest
A_EVI_2009_SC = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, ChangedVITS=AllEVITS, TargetYears=AllTargetYears)
FPStats(A_EVI_2009_SC) # 59%
GroupPlots(A_EVI_2009_SC, main="EVI with sctest, 2009-2017")
A_NDMI_2009_SC = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, ChangedVITS=AllNDMITS, TargetYears=AllTargetYears)
FPStats(A_NDMI_2009_SC) # 66%
GroupPlots(A_NDMI_2009_SC, main="NDMI with sctest, 2009-2017")
FPStats(A_EVI_2009_SC | A_NDMI_2009_SC) # 68%
A_NIRV_2009_SC = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, ChangedVITS=AllNIRVTS, TargetYears=AllTargetYears)
FPStats(A_NIRV_2009_SC)
GroupPlots(A_NIRV_2009_SC, main="NIRv with sctest, 2009-2017")

# Targetted sctest
A_EVI_2009_SCT = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, ChangedVITS=AllEVITS, TargetYears=AllTargetYears, scrange=c(2015.75, 2018.25))
FPStats(A_EVI_2009_SCT) # 50%
GroupPlots(A_EVI_2009_SCT, main="EVI with targetted sctest, 2009-2017")
A_NDMI_2009_SCT = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, ChangedVITS=AllNDMITS, TargetYears=AllTargetYears, scrange=c(2015.75, 2018.25))
FPStats(A_NDMI_2009_SCT) # 55%
GroupPlots(A_NDMI_2009_SCT, main="NDMI with targetted sctest, 2009-2017")
FPStats(A_EVI_2009_SCT | A_NDMI_2009_SCT) # 57%

# From 2014
A_EVI_2014 = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, ChangedVITS=AllEVITS, TargetYears=AllTargetYears, start=as.Date("2014-01-01"), scrange=NULL)
FPStats(A_EVI_2014) # 62%
GroupPlots(A_EVI_2014, main="EVI without sctest, 2014-2017")
A_NDMI_2014 = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, ChangedVITS=AllNDMITS, TargetYears=AllTargetYears, start=as.Date("2014-01-01"), scrange=NULL)
FPStats(A_NDMI_2014) # 65%
GroupPlots(A_NDMI_2014, main="NDMI without sctest, 2014-2017")
FPStats(A_EVI_2014 | A_NDMI_2014) # 65%

# Targetted sctest from 2014
A_EVI_2014_SCT = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, ChangedVITS=AllEVITS, TargetYears=AllTargetYears, start=as.Date("2014-01-01"), scrange=c(2015.75, 2018.25))
FPStats(A_EVI_2014_SCT) # 54%
GroupPlots(A_EVI_2014_SCT, main="EVI with targetted sctest, 2014-2017")
A_NDMI_2014_SCT = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, ChangedVITS=AllNDMITS, TargetYears=AllTargetYears, start=as.Date("2014-01-01"), scrange=c(2015.75, 2018.25))
FPStats(A_NDMI_2014_SCT) # 57%
GroupPlots(A_NDMI_2014_SCT, main="NDMI with targetted sctest, 2014-2017")
FPStats(A_EVI_2014_SCT | A_NDMI_2014_SCT) # 59%

# Loose targetted sctest from 2014
A_EVI_2014_SCTL = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, ChangedVITS=AllEVITS, TargetYears=AllTargetYears, start=as.Date("2014-01-01"), scrange=c(2015.75, 2018.25), scsig=0.1)
FPStats(A_EVI_2014_SCTL) # 57%
GroupPlots(A_EVI_2014_SCTL, main="EVI with targetted sctest (2016-2017), 2014-2017")
A_NDMI_2014_SCTL = sapply(1:nrow(AllPoints), TestMODBreakpointDetection, ChangedVITS=AllNDMITS, TargetYears=AllTargetYears, start=as.Date("2014-01-01"), scrange=c(2015.75, 2018.25), scsig=0.1)
FPStats(A_NDMI_2014_SCTL) # 58%
GroupPlots(A_NDMI_2014_SCTL, main="NDMI with targetted sctest (2016-2017), 2014-2017")
FPStats(A_EVI_2014_SCTL | A_NDMI_2014_SCTL) # 61%

## Overall this does not seem to help. Could other methods work better?

# Try Eline's t-test
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

# Monkey method
FPStats(rbinom(nrow(AllPoints), 1, 0.5)) # ~49%

# t-test
T_EVI_2009 = sapply(1:nrow(AllPoints), TestMODttest, AllEVITS)
FPStats(T_EVI_2009) # 63%
GroupPlots(T_EVI_2009, main="EVI t-test, 0.05 significance")
T_NDMI_2009 = sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS)
FPStats(T_NDMI_2009) # 69%
GroupPlots(T_NDMI_2009, main="NDMI t-test, 0.05 significance")
FPStats(T_EVI_2009 | T_NDMI_2009) # 70%

# And if we change the significance: tighter
T_EVI_2016T = sapply(1:nrow(AllPoints), TestMODttest, AllEVITS, sig=0.01)
FPStats(T_EVI_2016T) # 60%
GroupPlots(T_EVI_2016T, main="EVI t-test, 0.01 significance")
T_NDMI_2016T = sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS, sig=0.01)
GroupPlots(T_NDMI_2016T, main="NDMI t-test, 0.01 significance")
FPStats(T_NDMI_2016T) # 67%
FPStats(T_EVI_2016T | T_NDMI_2016T) # 68%

# Looser
T_EVI_2016L = sapply(1:nrow(AllPoints), TestMODttest, AllEVITS, sig=0.1)
FPStats(T_EVI_2016L) # 65%
GroupPlots(T_EVI_2016L, main="EVI t-test, 0.1 significance")
T_NDMI_2016L = sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS, sig=0.1)
FPStats(T_NDMI_2016L) # 69%
GroupPlots(T_NDMI_2016L, main="NDMI t-test, 0.1 significance")
FPStats(T_EVI_2016L | T_NDMI_2016L) # 70%

# Check around it
FPStats(sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS, sig=0.07))  # 68.9%
FPStats(sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS, sig=0.13))  # 69.0%
FPStats(sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS, sig=0.08))  # 69.0%
FPStats(sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS, sig=0.06))  # 68.8%
FPStats(sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS, sig=0.02))  # 67.9%
FPStats(sapply(1:nrow(AllPoints), TestMODttest, AllNDMITS, sig=0.005)) # 67.2%

# How about BFAST Monitor?

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

# Defaults
M_EVI_2009 = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS)
FPStats(M_EVI_2009) # 65%
GroupPlots(M_EVI_2009)
M_NDMI_2009 = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS)
FPStats(M_NDMI_2009) # 70%
GroupPlots(M_NDMI_2009, main="BFAST Monitor 2009-2017, defaults")
FPStats(M_EVI_2009 | M_NDMI_2009) # 71%
GroupPlots(M_EVI_2009 | M_NDMI_2009)

# Adjust significance
M_EVI_2009T = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS, level=0.001)
FPStats(M_EVI_2009T) # 62%
M_NDMI_2009T = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS, level=0.001)
FPStats(M_NDMI_2009T) # 70%
GroupPlots(M_NDMI_2009T)
FPStats(M_EVI_2009T | M_NDMI_2009T) # 71%
GroupPlots(M_EVI_2009T | M_NDMI_2009T)

# Adjust type
M_EVI_2009TC = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS, level=0.001, type="OLS-CUSUM")
FPStats(M_EVI_2009TC) # 64%
M_NDMI_2009TC = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS, level=0.001, type="OLS-CUSUM")
FPStats(M_NDMI_2009TC) # 70%
GroupPlots(M_NDMI_2009TC)
FPStats(M_EVI_2009TC | M_NDMI_2009TC) # 71%
GroupPlots(M_EVI_2009TC | M_NDMI_2009TC)

# Adjust h
M_EVI_2009TH5 = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS, level=0.001, h=0.5)
FPStats(M_EVI_2009TH5) # 53%
M_NDMI_2009TH5 = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS, level=0.001, h=0.5)
FPStats(M_NDMI_2009TH5) # 63%
FPStats(M_EVI_2009TH5 | M_NDMI_2009TH5) # 64%

M_EVI_2009TH1 = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS, level=0.001, h=1)
FPStats(M_EVI_2009TH1) # 37%
GroupPlots(M_EVI_2009TH1)
M_NDMI_2009TH1 = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS, level=0.001, h=1)
FPStats(M_NDMI_2009TH1) # 44%
FPStats(M_EVI_2009TH1 | M_NDMI_2009TH1) # 47%

# Add lag
M_EVI_2009TL = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS, level=0.001, lag=1, cloud_threshold=75)
FPStats(M_EVI_2009TL) # 61%
M_NDMI_2009TL = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS, level=0.001, lag=1, cloud_threshold=75)
FPStats(M_NDMI_2009TL) # 67%
FPStats(M_EVI_2009TL | M_NDMI_2009TL) # 70%

# Add slag
M_EVI_2009TS = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS, level=0.001, slag=1, cloud_threshold=75)
FPStats(M_EVI_2009TS) # 60%
M_NDMI_2009TS = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS, level=0.001, slag=1, cloud_threshold=75)
FPStats(M_NDMI_2009TS) # 68%
FPStats(M_EVI_2009TS | M_NDMI_2009TS) # 69%

# Trend only
M_EVI_2009TO = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS, level=0.001, formula=response~trend)
FPStats(M_EVI_2009TO) # 60%
M_NDMI_2009TO = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS, level=0.001, formula=response~trend)
FPStats(M_NDMI_2009TO) # 68%
FPStats(M_EVI_2009TO | M_NDMI_2009TO) # 69%

# Harmonics only
M_EVI_2009HO = sapply(1:nrow(AllPoints), TestMODMonitor, AllEVITS, level=0.001, formula=response~harmon)
FPStats(M_EVI_2009HO) # 49%
M_NDMI_2009HO = sapply(1:nrow(AllPoints), TestMODMonitor, AllNDMITS, level=0.001, formula=response~trend)
FPStats(M_NDMI_2009HO) # 68%
FPStats(M_EVI_2009HO | M_NDMI_2009HO) # 69%

# All in all, t-test wins... But signal to noise is still only 10%

# Combinations of breakpoints and t-test
FPStats(M_NDMI_2009) # 3.2%
FPStats(M_EVI_2009) # 3.7%
FPStats(T_NDMI_2009) # 9.6%
FPStats(T_EVI_2009) # 3.6%
FPStats(M_EVI_2009 & M_NDMI_2009) # 3.6%
FPStats(M_EVI_2009 | M_NDMI_2009) # 71%
FPStats(T_EVI_2009 & T_NDMI_2009) # 3.2%
FPStats(T_EVI_2009 | T_NDMI_2009) # 7.1%
FPStats(M_EVI_2009 & T_EVI_2009) # 2.9%
FPStats(M_EVI_2009 | T_EVI_2009) # 3.7%
FPStats(M_EVI_2009 & T_NDMI_2009) # 9.5%
FPStats(M_EVI_2009 | T_NDMI_2009) # 71%
FPStats(M_NDMI_2009 & T_EVI_2009) # 6.7%
FPStats(M_NDMI_2009 | T_EVI_2009) # 3.1%
FPStats(M_NDMI_2009 & T_NDMI_2009) # 9.5%
FPStats(M_NDMI_2009 | T_NDMI_2009) # 71%
GroupPlots(M_NDMI_2009 | T_NDMI_2009, main="NDMI BFAST Monitor OR t-test")
FPStats(M_EVI_2009 & M_NDMI_2009 & T_EVI_2009) # 0
FPStats(M_EVI_2009 | M_NDMI_2009 | T_EVI_2009) # 71%
FPStats(M_EVI_2009 & M_NDMI_2009 | T_EVI_2009) # 3.7%
FPStats(M_EVI_2009 | M_NDMI_2009 & T_EVI_2009) # 3.9%
FPStats(M_EVI_2009 & M_NDMI_2009 & T_NDMI_2009) # 7.1%
FPStats(M_EVI_2009 | M_NDMI_2009 | T_NDMI_2009) # 71%
FPStats(M_EVI_2009 & M_NDMI_2009 | T_NDMI_2009) # 5.1%
FPStats(M_EVI_2009 | M_NDMI_2009 & T_NDMI_2009) # 3.9%
FPStats(T_EVI_2009 & M_NDMI_2009 & T_NDMI_2009) # 0
FPStats(T_EVI_2009 | M_NDMI_2009 | T_NDMI_2009) # 71%
FPStats(T_EVI_2009 & M_NDMI_2009 | T_NDMI_2009) # 10.0%
FPStats(T_EVI_2009 | M_NDMI_2009 & T_NDMI_2009) # 5.2%
FPStats(M_EVI_2009 & T_NDMI_2009 & T_EVI_2009) # 0
FPStats(M_EVI_2009 | T_NDMI_2009 | T_EVI_2009) # 71%
FPStats(M_EVI_2009 & T_NDMI_2009 | T_EVI_2009) # 6.6%
FPStats(M_EVI_2009 | T_NDMI_2009 & T_EVI_2009) # 3.8%

FPStats(T_NDMI_2016T) # 10.3%
