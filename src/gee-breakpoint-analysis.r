source("gee-json-to-df.r")
CompleteDF = LoadGEEJSON()

# Convert to zoos?
library(zoo)
redzoo = zoo(CompleteDF$red, CompleteDF$time)
nirzoo = zoo(CompleteDF$nir, CompleteDF$time)
ndvizoo = (redzoo - nirzoo) / (redzoo + nirzoo)

# Split the complete DF into per-point DFs
# Needs AC.df from visualisation.r
# The coordinates are off by ~3m due to float precision!..
# Point 13 is covered by two overpasses.
h = ceiling(365.25/8)
for (i in 1:nrow(AC.df))
{
    print(i)
    Point1TS = CompleteDF[(abs(CompleteDF$longitude - AC.df[i,"sample_x"]) < 0.0002) &
                              (abs(CompleteDF$latitude - AC.df[i,"sample_y"]) < 0.0002),]
    Point1TS = TemporalDailyMerge(Point1TS)
    #print(table(year(Point1TS$time)))
    #print(nrow(Point1TS))
    #next
    stopifnot(!any(duplicated(Point1TS$time)))
    print(paste("Observations for point", i, ":", nrow(Point1TS)))
    plot(evi~time, data=Point1TS, type="l", main=paste("EVI", i)); abline(v=as.POSIXct(as.Date("2017-01-01")))
    if (nrow(Point1TS) > h*2) {
        P1TS = bfastts(Point1TS$evi, Point1TS$time, "irregular")
        P1PP = bfastpp(P1TS, order=3)
        P1BP = breakpoints(formula(response~trend+harmon), data=P1PP, h=h)
        if (length(P1BP$breakpoints) > 0) {
            abline(v=date_decimal(P1PP[P1BP$breakpoints, "time"]), col="red")
        }
        print(P1BP)
    } else print("too cloudy")
    plot(ndmi~time, data=Point1TS, type="l", main=paste("NDMI", i)); abline(v=as.POSIXct(as.Date("2017-01-01")))
    if (nrow(Point1TS) > h*2) {
        P1TS = bfastts(Point1TS$ndmi, Point1TS$time, "irregular")
        P1PP = bfastpp(P1TS, order=3)
        P1BP = breakpoints(formula(response~trend+harmon), data=P1PP, h=h)
        if (length(P1BP$breakpoints) > 0) {
            abline(v=date_decimal(P1PP[P1BP$breakpoints, "time"]), col="red")
        }
        print(P1BP)
    } else print("too cloudy")
}
# Visible change in 3, 5, 8, 9 but not others!..
# Coordinates seem to be shifted by ~0.0002 degrees; likely this is the centroid of the pixel

TemporalDailyMerge = function(ts.df)
{
    Dates = unique(date(ts.df$time))
    MergeColumns = c("blue", "red", "nir", "swir", "ndvi", "ndmi", "evi")
    NewDF = NULL
    for (i in 1:length(Dates))
    {
        if (month(Dates[i]) == 2 && day(Dates[i]) == 29)
            next # Discard leap days for now
        DayRows = ts.df[date(ts.df$time) %in% Dates[i],]
        DayRows[1,MergeColumns] = colMeans(DayRows[,MergeColumns])
        DayRows = DayRows[1,] # Note: this discards the QA from the second pixel
        NewDF = rbind(NewDF, DayRows)
    }
    return(NewDF)
}

# Locations of subpixels:
cbind(AC.df[,c("subpixel_center_x", "subpixel_center_y", "name", "sample_x", "sample_y")], id=1:34)

# Per-point analysis:
#  1: Questionable break detected in 2016
#  2: Nothing obvious
#  3: Slightly early detected in NDMI, even earlier in EVI
#  4: Nothing obvious
#  5: Break detected in both
#  6: Strange seasonality, (questionable) break not detected
#  7: Questionable early break detected in both
#  8: Break detected in both early, EVI is more obvious
#  9: NDMI chaotic, EVI shows change but not detected (interrupted regrowth)
# 10: Too cloudy
# 11: Too cloudy
# 12: Questionable break detected early on both + cloudy
# 13: Break not detected, but might indeed be there (regrowth)
# 14: Break detected in both but a bit early
# 15: Break detected in both, but highly questionable
# 16: Break detected in both, EVI a bit early.
# 17: Steady regrowth. EVI break detected a bit early, break visible in NDMI but not detected.
# 18: Very cloudly, nothing obvious. Slight change in NDMI in 2017 if fewer clouds filtered.
# 19: Nothing obvious, NDMI chaotic. When displaced, a little evidence in EVI.
# 20: Nothing obvious; maybe a bit lower in 2016 but not really
# 21: Something detected in NDMI, but too cloudy to tell.
# 22: Very cloudy, hard to say anything.
# 23: Also nothing obvious. When displaced it might be changed, but it's questionable.
# 24: Detected in NDMI, though a bit too late. Obvious but not detected in EVI
# 25: Nothing obvious in both. When displaced, there is slight change evidence.
# 26: Same as 27
# 27: Nothing detected nor obvious. NDMI is chaotic, EVI shows no change.
# 28: Nothing visible nor detected. No evidence of change around either!
# 29: Detected in NDMI only, arguable by eye
# 30: Nothing is visible nor detected. Some evidence in NDMI in the neighbourhood.
# 31: Detected break correctly, though it's not obvious by eye
# 32: Detected huge break in both!
# 33: Detected huge break in both!
# 34: obvious break in 2014 in both, 2017 (undetected) in EVI only.
#     The real change is displaced a bit to the east (-10.755955, 13.098332)! And not entirely obvious.
#     And yet it's supposed to be a full pixel change...
# Overall there are a lot fewer detections, not everything is detected but it's better on average.
# Detecting regrowth is problematic

WinnerLandsat = c(T, F, T, F, T, F, T, T, F, F, F, T, F, T, T, T, T, F, F, F, T, F, F, T, F, F, F, F, T, F, T, T, T, F)
mean(WinnerLandsat) # 45% accuracy, but hugely fewer false detections (almost none)
# If we only take 2016
Winner2016 = WinnerLandsat[AC.df[,"name"] == "LC change 2016"]
mean(Winner2016) # 46% accuracy

### Quick automated check
TestBreakpointDetection = function(refpoint, vi = "ndmi", h = c("static", "dynamic"), threshold=0.25, start=as.POSIXct("2009-01-01"))
{
    h = match.arg(h)
    Point1TS = CompleteDF[(abs(CompleteDF$longitude - refpoint[["sample_x"]]) < 0.0002) &
                              (abs(CompleteDF$latitude - refpoint[["sample_y"]]) < 0.0002),]
    Point1TS = Point1TS[Point1TS$time > start,]
    Point1TS = TemporalDailyMerge(Point1TS)
    h = switch(h, static=ceiling(365.25/8), dynamic=max(table(year(Point1TS$time))))
    stopifnot(!any(duplicated(Point1TS$time)))
    print(paste("Observations for point", refpoint[["validation_id"]], ":", nrow(Point1TS)))
    #plot(evi~time, data=Point1TS, type="l", main=paste("EVI", i)); abline(v=as.POSIXct(as.Date("2017-01-01")))
    TargetYear = switch(refpoint[["name"]], "LC change 2016"=2016, "LC change 2017"=2017)
    if (nrow(Point1TS) > h*2 && nrow(Point1TS) > 42) {
        P1TS = bfastts(Point1TS[[vi]], Point1TS$time, "irregular")
        P1PP = bfastpp(P1TS, order=3)
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

# Also check only breaks in 2016
idx2016 = which(AC.df$name == "LC change 2016")

# Static h
NDMI_h46 = apply(AC.df, 1, TestBreakpointDetection)
mean(NDMI_h46, na.rm=TRUE) # 38%
mean(NDMI_h46[idx2016], na.rm=TRUE) # 46% - winner for fewest breaks, never more than 2
EVI_h46 = apply(AC.df, 1, TestBreakpointDetection, vi="evi")
mean(EVI_h46, na.rm=TRUE) # 25%
mean(EVI_h46[idx2016], na.rm=TRUE) # 33%
Both_h46 = NDMI_h46 | EVI_h46
mean(Both_h46, na.rm=TRUE) # Doesn't help
mean(Both_h46[idx2016], na.rm=TRUE) # Doesn't help
NDVI_h46 = apply(AC.df, 1, TestBreakpointDetection, vi="ndvi")
mean(NDVI_h46, na.rm=TRUE) # 28%
mean(NDVI_h46 | NDMI_h46, na.rm=TRUE) # 40%, a very slight improvement

# Dynamic h
NDMI_hd = apply(AC.df, 1, TestBreakpointDetection, h="dynamic")
mean(NDMI_hd, na.rm=TRUE) # 41%
mean(NDMI_hd[idx2016], na.rm=TRUE) # 50%
EVI_hd = apply(AC.df, 1, TestBreakpointDetection, vi="evi", h="dynamic")
mean(EVI_hd, na.rm=TRUE) # 50%
mean(EVI_hd[idx2016], na.rm=TRUE) # 63%
Both_hd = NDMI_hd | EVI_hd
mean(Both_hd, na.rm=TRUE) # 53%
mean(Both_hd[idx2016], na.rm=TRUE) # 67% - winner with stricter threshold
NDVI_hd = apply(AC.df, 1, TestBreakpointDetection, vi="ndvi", h="dynamic")
mean(NDVI_hd, na.rm=TRUE) # 31%
mean(NDVI_hd | EVI_hd, na.rm=TRUE) # 53% again
mean(NDVI_hd | EVI_hd | NDMI_hd, na.rm=TRUE) # 56%

# So lowering h detects more breaks, NDVI is useless, sometimes EVI and sometimes NDMI is better,
# combining both doesn't help much but is safe

# What it looks like with a half-year threshold
# Static h
NDMI_h46 = apply(AC.df, 1, TestBreakpointDetection, threshold=0.5)
mean(NDMI_h46, na.rm=TRUE) # 41%
mean(NDMI_h46[idx2016], na.rm=TRUE) # 46%
EVI_h46 = apply(AC.df, 1, TestBreakpointDetection, vi="evi", threshold=0.5)
mean(EVI_h46, na.rm=TRUE) # 34%
mean(EVI_h46[idx2016], na.rm=TRUE) # 38%
Both_h46 = NDMI_h46 | EVI_h46
mean(Both_h46, na.rm=TRUE) # 43%
mean(Both_h46[idx2016], na.rm=TRUE) # 50%
NDVI_h46 = apply(AC.df, 1, TestBreakpointDetection, vi="ndvi", threshold=0.5)
mean(NDVI_h46, na.rm=TRUE) # 31%
mean(NDVI_h46[idx2016], na.rm=TRUE) # 29%
mean(NDVI_h46 | NDMI_h46, na.rm=TRUE) # 44%, a very slight improvement
mean(NDVI_h46[idx2016] | NDMI_h46[idx2016], na.rm=TRUE) # 46%, no improvement
# Dynamic h
NDMI_hd = apply(AC.df, 1, TestBreakpointDetection, h="dynamic", threshold=0.5)
mean(NDMI_hd, na.rm=TRUE) # 41%
mean(NDMI_hd[idx2016], na.rm=TRUE) # 50%
EVI_hd = apply(AC.df, 1, TestBreakpointDetection, vi="evi", h="dynamic", threshold=0.5)
mean(EVI_hd, na.rm=TRUE) # 56%
mean(EVI_hd[idx2016], na.rm=TRUE) # 67%
Both_hd = NDMI_hd | EVI_hd
mean(Both_hd, na.rm=TRUE) # 59%
mean(Both_hd[idx2016], na.rm=TRUE) # 71% - winner
NDVI_hd = apply(AC.df, 1, TestBreakpointDetection, vi="ndvi", h="dynamic", threshold=0.5)
mean(NDVI_hd, na.rm=TRUE) # 38%
mean(NDVI_hd[idx2016], na.rm=TRUE) # 38%
mean(NDVI_hd | EVI_hd, na.rm=TRUE) # 59% as well
mean(NDVI_hd[idx2016] | EVI_hd[idx2016], na.rm=TRUE) # 67% as well
mean(NDVI_hd | EVI_hd | NDMI_hd, na.rm=TRUE) # 63%
mean(NDVI_hd[idx2016] | EVI_hd[idx2016] | NDMI_hd[idx2016], na.rm=TRUE) # 71%

# What if we start from 2014?
# Static h
NDMI_h46_14 = apply(AC.df, 1, TestBreakpointDetection, start=as.POSIXct("2014-01-01"))
mean(NDMI_h46_14, na.rm=TRUE) # 39%
mean(NDMI_h46_14[idx2016], na.rm=TRUE) # 45% - winner for fewest breaks, only once 3
EVI_h46_14 = apply(AC.df, 1, TestBreakpointDetection, vi="evi", start=as.POSIXct("2014-01-01"))
mean(EVI_h46_14, na.rm=TRUE) # 36%
mean(EVI_h46_14[idx2016], na.rm=TRUE) # 40%
Both_h46_14 = NDMI_h46_14 | EVI_h46_14
mean(Both_h46_14, na.rm=TRUE) # 42%
mean(Both_h46_14[idx2016], na.rm=TRUE) # 50%, much better
NDVI_h46_14 = apply(AC.df, 1, TestBreakpointDetection, vi="ndvi", start=as.POSIXct("2014-01-01"))
mean(NDVI_h46_14, na.rm=TRUE) # 29%
mean(NDVI_h46_14 | NDMI_h46_14, na.rm=TRUE) # 39%, doesn't help

# Dynamic h
NDMI_hd_14 = apply(AC.df, 1, TestBreakpointDetection, h="dynamic", start=as.POSIXct("2014-01-01"))
mean(NDMI_hd_14, na.rm=TRUE) # 44%
mean(NDMI_hd_14[idx2016], na.rm=TRUE) # 50%
EVI_hd_14 = apply(AC.df, 1, TestBreakpointDetection, vi="evi", h="dynamic", start=as.POSIXct("2014-01-01"))
mean(EVI_hd_14, na.rm=TRUE) # 43%
mean(EVI_hd_14[idx2016], na.rm=TRUE) # 50%
Both_hd_14 = NDMI_hd_14 | EVI_hd_14
mean(Both_hd_14, na.rm=TRUE) # 53%
mean(Both_hd_14[idx2016], na.rm=TRUE) # 58% - winner with stricter threshold, but not by much
NDVI_hd_14 = apply(AC.df, 1, TestBreakpointDetection, vi="ndvi", h="dynamic", start=as.POSIXct("2014-01-01"))
mean(NDVI_hd_14, na.rm=TRUE) # 34%
mean(NDVI_hd_14 | EVI_hd_14, na.rm=TRUE) # 50% still
mean(NDVI_hd_14 | EVI_hd_14 | NDMI_hd_14, na.rm=TRUE) # 56%
mean(NDVI_hd_14[idx2016] | EVI_hd_14[idx2016] | NDMI_hd_14[idx2016], na.rm=TRUE) # 58%

# Overall fewer breaks are detected
