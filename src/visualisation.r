# Scratchpad for visualising POI time series
# Run most of detect-breaks to get the TS

# Get POI
library(sf)
POI = st_read("../data/POI.geojson")
POI = st_transform(POI, crs(timeseries)@projargs)
spPOI = as(POI, "Spatial")
poi_val = extract(timeseries, as(spPOI, "SpatialPoints"), cellnumbers=TRUE)
plot(POI)
poi_val[,1] # Cell numbers; not used and somehow seem to be in reverse?!
bfts = bfastts(poi_val[1,-1], dates, type = "10-day")
plot(na.approx(bfts))
bpp = bfastpp(bfts, order=2)
#bpp$ts = bfts
bfr = breakpoints(response ~ (harmon + trend), data=bpp, h=36)
bfci = confint(bfr)
#bft = breakpoints(ts ~ (harmon + trend), data=bpp, h=GetBreakNumber(dates))
#modelterms <- terms(formula, data = data)
#X <- model.matrix(modelterms, data = data)
#n <- nrow(X)
#orig.y = eval(attr(terms(formula), "variables")[[2]], data, env)

## Parse breakpoints and breakpoints.confint into yearly layers and associated confint times
# How many layers(=years) there should be
DateRange = range(dates)
Years = year(DateRange[1]):year(DateRange[2])
OutMatrix = matrix(-1, nrow=length(Years), ncol=3, dimnames=list(Years, c("breakpoint", "confint.neg", "confint.pos")))
ConfInts = bfci$confint
BreakpointYears = as.integer(sapply(ConfInts[,"breakpoints"], BreakpointDate, bpp))
BreakpointDays = sapply(ConfInts, BreakpointToDateSinceT0, bpp)
OutMatrix[rownames(OutMatrix) %in% BreakpointYears,] = BreakpointDays
c(t(OutMatrix))

BreakpointToDateSinceT0 = function(breakpoint_index, bpp)
{
    as.integer(as.Date(date_decimal(BreakpointDate(breakpoint_index, bpp))) - t0)
}

# The date of the breakpoint in decimal years
BreakpointDate = function(breakpoint_index, bpp)
{
    bpp$time[breakpoint_index]
}



# Can visualise a bpp
plot(read.zoo(bpp))

## For 20350392
# pixel[35] = mean(pixel[34], pixel[36])
# pixel[79] = mean(pixel[78], pixel[80])
# bfts = bfastts(pixel, dates, type = "16-day")
visualise_breakpoints = function(point_id, order=3, timestep="regular", ...)
{
    pixel = poi_val[point_id,-1]
    if (timestep == "regular")
    {
        z = zoo(pixel,dates) # make a zoo
        bfts = as.ts(z)
        print(bfts)
        dates = as.Date(as.numeric(time(bfts)))
    } else {
        bfts = bfastts(pixel, dates, type = timestep)
    }
    bpp = .bfastpp.full(bfts, order=order)
    #plot(read.zoo(bpp))
    bf = breakpoints(response ~ (harmon + trend), data=bpp, h=GetBreakNumber(dates))
    print(bf$breakpoints)
    print(as.integer(bpp$time[[max(bf$breakpoints)]] - as.Date("2014-03-12"))) ## Gives NA for irregular because dates do not correspond to bpp$trend in this case
    print(length(dates))
    print(length(bfts))
    print(nrow(bpp))
    plot(response/10000~time, data=bpp, col=ifelse(getSceneinfo(names(timeseries))$sensor == "OLI", "blue", "green"), ...)#main=POI$comment[point_id], ...)
    for (i in 1:length(bf$breakpoints))
        abline(v=bpp$time[[bf$breakpoints[[i]]]], col="red")
    return(bf)
}
bf = visualise_breakpoints(3, type="l", timestep="regular", ylab="NDVI")
for (i in 1:length(dates))
        abline(v=dates[[i]], col="blue")

breaks = bfast(na.approx(bfts), season="harmonic", max.iter=2, h=GetBreakNumber(getZ(timeseries)))
plot(breaks)


## Test making a regular 8-day time series
pixel = poi_val[1,-1]
	z <- zoo(pixel,dates)
	yr <- as.numeric(format(time(z), "%Y"))
	jul <- as.numeric(format(time(z), "%j"))
	delta <- min(unlist(tapply(jul, yr, diff))) # 8
	zz <- aggregate(z, yr + (jul - 1) / delta / (365.25/8))
	plot(yr+(jul-1)/365.25)
	(tso <- as.ts(zz))
plot(ts(pixel, frequency = 365.25/8, start = c(yr[1], jul[1]/8)))
## jul is not divisible by 8... But the delta is 8

## Solution: use a ts that is in days since 1970 with a frequency of 1/8
z = zoo(pixel,dates) # make a zoo
tso = as.ts(z) # make a ts: this fills in the gaps!
ts_dates = as.Date(as.numeric(time(rts))) # Convert days back into dates to get a new "dates" object
# Plot time series with correct date labels
plot(ts_dates, tso, type="l")

# Sizes of the objects:
object_size(tso) # 2.22 kB
object_size(bfts) # 14.6 kB
as.numeric(object_size(bfts) / object_size(tso)) # 6.6 times smaller!

## Overwrite bfastpp with a function that does not use frequency()
.bfastpp.full <- function(data, order = 3,
                    lag = NULL, slag = NULL, na.action = na.omit,
                    stl = c("none", "trend", "seasonal", "both"))
{
  ## double check what happens with 29-02 if that happens...
  ## we should keep it simple an remove the datum if that happens
  
  if(!is.ts(data)) data <- as.ts(data)
  
  ## STL pre-processing to try to adjust for trend or season
  stl <- match.arg(stl)
  if(stl != "none") {
    stl_adjust <- function(x) {
      x_stl <- stats::stl(x, s.window = "periodic")$time.series
      switch(stl,
             "trend" = x - x_stl[, "trend"],
             "seasonal" = x - x_stl[, "seasonal"],
             "both" = x - x_stl[, "trend"] - x_stl[, "seasonal"])
    }
    if(NCOL(data) > 1L) {
      for(i in 1:NCOL(data)) data[,i] <- stl_adjust(data[,i])
    } else {
      data <- stl_adjust(data)
    }
  }
  
  ## check for covariates
  if(NCOL(data) > 1L) {
    x <- coredata(data)[, -1L]
    y <- data[, 1L]
  } else {
    x <- NULL
    y <- data
  }
  
  freq <- frequency(y)
  dec_dates = decimal_date(as.Date(as.numeric(time(y))))
  ## data with trend and season factor
  rval <- data.frame(
    time = DateFromTS(y),#as.numeric(time(y)),
    response = y,
    trend = 1:NROW(y),
    season = if(freq >= 1) factor(cycle(y)) else dec_dates %% 1
  )
  
  ## set up harmonic trend matrix as well
  
  #order <- min(freq, order)
  if (freq >= 1) {
    harmon <- outer(2 * pi * as.vector(time(y)), 1:order)
  } else {
    harmon <- outer(2 * pi * dec_dates, 1:order)
  }
  harmon <- cbind(apply(harmon, 2, cos), apply(harmon, 2, sin))
  colnames(harmon) <- if(order == 1) {
    c("cos", "sin")
  } else {
    c(paste("cos", 1:order, sep = ""), paste("sin", 1:order, sep = ""))
  }
  if((2 * order) == freq) harmon <- harmon[, -(2 * order)]
  rval$harmon <- harmon
  
  ## add lags
  nalag <- function(x, k) c(rep(NA, k), head(x, -k))
  if(!is.null(lag)) {
    rval$lag <- sapply(lag, function(k) nalag(as.vector(y), k))
    colnames(rval$lag) <- lag
  }
  if(!is.null(slag)) {
    rval$slag <- sapply(slag * freq, function(k) nalag(as.vector(y), k))
    colnames(rval$slag) <- slag
  }
  
  ## add regressors
  rval$xreg <- x
  
  ## omit missing values
  rval <- na.action(rval)
  
  ## return everything
  return(rval)
}

library(lubridate)
dates_dec = decimal_date(ts_dates)
harmon <- outer(2 * pi * dates_dec, 1:3)

## Get a Date vector from a ts object.
DateFromTS = function(timeseries)
{
    # Get the timestamps; we don't yet know what they mean.
    timestamps = time(timeseries)
    
    # Case 1: the dates are fractional years, like in the "irregular" case.
    if (frequency(timeseries) > 1)
    {
        return(as.Date(date_decimal(as.numeric(timestamps))))
    } else { # Case 2: the dates are the number of days since 1970
        return(as.Date(as.numeric(timestamps)))
    }
}

#######################
# Check the results of the breakpoint validation

breakval = st_read("../data/valgroup27_exports_20180904_LandCoverChangeDetection_1.csv", options=c("X_POSSIBLE_NAMES=sample_x", "Y_POSSIBLE_NAMES=sample_y"))
st_crs(breakval) = 4326
ChangedPixels = breakval[breakval$legend_name=="CGLOPS Land Cover Change Detection validation - Did the LC type at 100m x 100m change?",]
plot(ChangedPixels["name"])
levels(ChangedPixels$name)
plot(ChangedPixels$name)
nrow(ChangedPixels[ChangedPixels$name!="No",])
plot(ChangedPixels[ChangedPixels$name!="No",]$name)
plot(ChangedPixels[ChangedPixels$name!="No","name"])

# All subpixels with marked change
FilteredPoints = breakval[breakval$legend_name=="CGLOPS Land Cover Change Detection validation - Mark LC change",]
plot(FilteredPoints$name)
unique(FilteredPoints$name)
FilteredPoints
# Get one value per MODIS pixel
FilteredPoints = FilteredPoints[!duplicated(FilteredPoints[,c("name", "sample_x", "sample_y")]),]
duplicated(FilteredPoints[,c("sample_x", "sample_y")]) # No points where changed happened during two years

ChangedOnly = ChangedPixels[ChangedPixels$name!="No",]
ChangedOnly$geometry %in% FilteredPoints$geometry # There are two points that are labelled as change but no subpixels
ExtraPoints = ChangedOnly[!ChangedOnly$geometry %in% FilteredPoints$geometry,] # They are complete change in 2016, so relabel as such
ExtraPoints$name = factor("LC change 2016", levels = levels(ChangedOnly$name))

AllChanged = rbind(FilteredPoints, ExtraPoints)
any(duplicated(AllChanged[,c("sample_x", "sample_y")])) # No duplicates
rm(breakval)

# Save back a trimmed file with relevant parts
st_write(ChangedPixels[c("name", "sample_x", "sample_y")], "../data/breakpoint_validation.csv")
st_write(AllChanged[c("legend_name", "name", "confidence", "sample_x", "sample_y")], "../data/breakpoint_changed_pixels.csv")

# Check which model performs the best
#                                  | First is bad, others are almost accurate
WinsEVI2   = c(T, F, F, F, T, T, F, F, T, F, F, T, T, F, F, F, T, F, T, F, F, F, T, F, F, T, T)
WinsEVI3   = c(T, F, F, F, T, T, T, F, T, T, F, T, T, F, F, T, T, F, T, F, F, F, F, T, T, F, T)
WinsNDMI2  = c(T, F, F, F, F, F, F, F, T, F, F, F, F, T, T, F, T, F, F, T, F, T, T, F, F, T, T)
WinsNDMI3  = c(T, F, F, F, F, F, T, T, T, F, F, F, F, T, T, F, T, F, F, T, F, T, T, T, T, T, T)
WinsHansen = c(F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F, F)
#                    |     | Detected as 2014 instead 
#                       |           | Detected correctly but 130/160 m away
#                                                     | Detected as 2013

Winners = cbind(WinsEVI2, WinsEVI3, WinsNDMI2, WinsNDMI3)
colSums(Winners) # Draw between EVI3 and NDMI3
colMeans(Winners) # 52% accuracy
mean(apply(Winners, 1, any)) # If we take all, we get 78% accuracy
mean(apply(cbind(WinsEVI3, WinsNDMI3), 1, any)) # And we lose nothing if we omit all 2nd order models

# Load the relevant raster time series
source("utils/dates.r")
library(gdalUtils)
library(raster)
library(bfast)
library(strucchange)
library(zoo)
EVIMosaicFile = "/data/users/Public/greatemerald/modis/raw-input-mosaic-EVI.vrt"
if (!file.exists(EVIMosaicFile))
{
    InputFiles = c("/data/mep_cg1/MOD_S10/additional_VIs_new/X16Y06/MOD_S10_TOC_X16Y06_20090101-20171231_250m_C6_EVI.tif",
                   "/data/mep_cg1/MOD_S10/additional_VIs_new/X17Y06/MOD_S10_TOC_X17Y06_20090101-20171231_250m_C6_EVI.tif",
                   "/data/mep_cg1/MOD_S10/additional_VIs_new/X18Y06/MOD_S10_TOC_X18Y06_20090101-20171231_250m_C6_EVI.tif",
                   "/data/mep_cg1/MOD_S10/additional_VIs_new/X19Y06/MOD_S10_TOC_X19Y06_20090101-20171231_250m_C6_EVI.tif")
    gdalbuildvrt(InputFiles, EVIMosaicFile)
}
EVIMosaic = brick(EVIMosaicFile)
EVITSFile = "../data/ChangedEVITS.csv"
if (!file.exists(EVITSFile))
{
    ChangedEVITS = extract(EVIMosaic, AllChanged)
    write.csv(ChangedEVITS, EVITSFile)
} else ChangedEVITS = as.matrix(read.csv(EVITSFile, row.names=1))

NDMIMosaicFile = "/data/users/Public/greatemerald/modis/raw-input-mosaic-ndmi.vrt"
if (!file.exists(NDMIMosaicFile))
{
    InputFiles = c("/data/mep_cg1/MOD_S10/additional_VIs_new/X16Y06/MOD_S10_TOC_X16Y06_20090101-20171231_250m_C6_NDMI.tif",
                   "/data/mep_cg1/MOD_S10/additional_VIs_new/X17Y06/MOD_S10_TOC_X17Y06_20090101-20171231_250m_C6_NDMI.tif",
                   "/data/mep_cg1/MOD_S10/additional_VIs_new/X18Y06/MOD_S10_TOC_X18Y06_20090101-20171231_250m_C6_NDMI.tif",
                   "/data/mep_cg1/MOD_S10/additional_VIs_new/X19Y06/MOD_S10_TOC_X19Y06_20090101-20171231_250m_C6_NDMI.tif")
    gdalbuildvrt(InputFiles, NDMIMosaicFile)
}
NDMIMosaic = brick(NDMIMosaicFile)
NDMITSFile = "../data/ChangedNDMITS.csv"
if (!file.exists(NDMITSFile))
{
    ChangedNDMITS = extract(NDMIMosaic, AllChanged)
    write.csv(ChangedNDMITS, NDMITSFile)
} else ChangedNDMITS = as.matrix(read.csv(NDMITSFile, row.names=1))
# The extraction is correct, no off-by-one issues here, checked in QGIS with original data

library(lubridate)
AC.df = AllChanged
class(AC.df) = "data.frame"
dates = GetDatesFromDir("/data/mep_cg1/MOD_S10/MOD_S10")
plot(ChangedEVITS[1,]~dates, type="l", main=paste(1, AC.df[1,"confidence"])); abline(v=as.Date("2016-01-01")); abline(v=as.Date("2017-01-01"))

# Plots
for (i in 1:34) {
    if (!all(is.na(ChangedEVITS[i,])))
        plot(ChangedEVITS[i,]~dates, type="l", main=paste(i, AllChanged$name[i], "EVI", AC.df[i,"confidence"])); abline(v=as.Date("2016-01-01")); abline(v=as.Date("2017-01-01"))
}
for (i in 1:34) {
    if (!all(is.na(ChangedNDMITS[i,])))
        plot(ChangedNDMITS[i,]~dates, type="l", main=paste(i, AllChanged$name[i], "NDMI", AC.df[i,"confidence"])); abline(v=as.Date("2016-01-01")); abline(v=as.Date("2017-01-01"))
}

scpvals = function(ts, from=0)
{
    if (all(is.na(ts)))
        return(NA)
    bfts = bfastts(ts, dates, type="10-day")
    bpp = bfastpp(bfts, order=3)
    bpp = bpp[bpp$time > from,]
    return(sctest(efp(response ~ (harmon + trend), data=bpp, h=GetBreakNumber(dates), type="OLS-MOSUM"))$p.value)
}

EVI_pvals = apply(ChangedEVITS, 1, scpvals)
NDMI_pvals = apply(ChangedNDMITS, 1, scpvals)
EVI_p_sig = EVI_pvals  < 0.05 # Considered breaks
NDMI_p_sig = NDMI_pvals < 0.05 # Ideally all of these would be true; except that in most cases we can't even see the break by eye
sum(EVI_p_sig) # One more
sum(NDMI_p_sig)
which(EVI_p_sig != NDMI_p_sig) # Disagreements
# Those that are FALSE will never be winners

EVI_pvals_14 = apply(ChangedEVITS, 1, scpvals, from=2014)
NDMI_pvals_14 = apply(ChangedNDMITS, 1, scpvals, from=2014)
EVI_p_sig_14 = EVI_pvals_14  < 0.05 # Considered breaks
NDMI_p_sig_14 = NDMI_pvals_14 < 0.05 # Ideally all of these would be true; except that in most cases we can't even see the break by eye
sum(EVI_p_sig_14) # Two more
sum(NDMI_p_sig_14)
which(EVI_p_sig_14 != NDMI_p_sig_14) # Disagreements

# Filter out 2017
EVI_pvals_2016 = EVI_pvals[AllChanged$name=="LC change 2016"]
length(EVI_pvals_2016) == length(WinsEVI3)
length(AllChanged$name)
length(EVI_pvals)

WinsEVI3[!EVI_pvals_2016 < 0.05]

## Draw EVI point 9 with harmonics
MyTS = bfastts(ChangedEVITS[9,], dates, "10-day")
#MyTS = na.approx(MyTS) # This step is optional, just for visualisation
MyBF = bfast(MyTS, GetBreakNumber(dates), season="harmonic", max.iter=3)
MyIter = length(MyBF$output)
plot(MyBF)
plot(MyBF, "seasonal")
plot(MyTS)
plot(MyBF$output)
MyFit = MyBF$output[[MyIter]]$St + MyBF$output[[MyIter]]$Tt
MyTrend = MyBF$output[[MyIter]]$Tt
plot(MyFit)
MyTimes = as.numeric(time(MyBF$output[[MyIter]]$Tt))
MyTimes = MyTimes[!is.na(MyBF$output[[MyIter]]$Tt)]
MyBreaks = MyTimes[MyBF$output[[MyIter]]$bp.Vt$breakpoints]
plot(MyTS, ylab="EVI", main="Pixel with a confirmed break in 2016"); lines(MyFit, col="blue"); lines(MyTrend, col="green"); abline(v=MyBreaks, col="red")

MyPlotBfast = function(MyBF, PlotData=TRUE, ...)
{
    MyIter = length(MyBF$output)
    MyFit = MyBF$output[[MyIter]]$St + MyBF$output[[MyIter]]$Tt
    MyTrend = MyBF$output[[MyIter]]$Tt
    MyTimes = as.numeric(time(MyBF$output[[MyIter]]$Tt))
    MyTimes = MyTimes[!is.na(MyBF$output[[MyIter]]$Tt)]
    if (any(!is.na(MyBF$output[[MyIter]]$bp.Vt)))
    {
        MyBreaks = MyTimes[MyBF$output[[MyIter]]$bp.Vt$breakpoints]
    } else {
        MyBreaks = NULL
    }
    if (PlotData)
    {
        plot(MyBF$Yt, ...)
        lines(MyFit, col="blue")
    } else {
        plot(MyFit, col="blue", ...)
    }
    lines(MyTrend, col="green")
    abline(v=MyBreaks, col="red")
}

for (i in which(EVI_p_sig != NDMI_p_sig)) {
    if (!all(is.na(ChangedEVITS[i,])))
    {
        MyBF = try(bfast(bfastts(na.approx(ChangedEVITS[i,]), dates, "10-day"), GetBreakNumber(dates), season="harmonic", max.iter=1))
        if (class(MyBF) == "bfast")
            MyPlotBfast(MyBF, ylab="EVI", main=paste(i, AllChanged$name[i], "EVI", AC.df[i,"confidence"])); abline(v=2016); abline(v=2017)
    }
}

for (i in which(EVI_p_sig != NDMI_p_sig)) {
    if (!all(is.na(ChangedNDMITS[i,])))
    {
        MyBF = try(bfast(bfastts(na.approx(ChangedNDMITS[i,]), dates, "10-day"), GetBreakNumber(dates), season="harmonic", max.iter=1))
        if (class(MyBF) == "bfast")
            MyPlotBfast(MyBF, ylab="NDMI", main=paste(i, AllChanged$name[i], "NDMI", AC.df[i,"confidence"])); abline(v=2016); abline(v=2017)
    }
}

# Interesting cases:
#   32: clear break, but a lot more spurious breaks are detected
#   29: outliers
#   26: a lot of false change, and some true change at the end
#   24: Clear change at the end
#   22/23: No change lately
#   28/21: no breaks
#   16/19: all breaks during peaks, this is very suspect. Not real change
#   15: nice big change at the right time; but none in Landsat!
#   14: too many breaks + actual real break
#   13: A lot of variation, good example
#    9: too many breaks, but a clear example
#    6/8: Also a lot of variation
#    2: strong example without too many breaks

MyBF = bfast(bfastts(na.approx(ChangedEVITS[1,]), dates, "10-day"), GetBreakNumber(dates), season="harmonic", max.iter=1)
MyPlotBfast(MyBF, TRUE, ylab="EVI", main="BFAST time series segmentation")

####################
# LPS19 poster

AllChangedPoints = LoadAllReference()
DeforestationPoint = AllChangedPoints[AllChangedPoints$location_id == 2187641,]
DeforestationTS = extract(EVIMosaic, DeforestationPoint)
DamPoints = AllChangedPoints[AllChangedPoints$location_id %in% c(2188595, 2187772, 2188293),]
DamTS = extract(EVIMosaic, DamPoints)

poi_val = DeforestationTS[1,]

plot(poi_val, type="l")

poi_val = DamTS[1,]
dates = dates[1:length(DeforestationTS)]
bts = bfastts(poi_val, dates, type="10-day")
bpp = bfastpp(bts, 3)
bf = breakpoints(response ~ (harmon + trend), data=bpp, h=GetBreakNumber(dates))
print(bf$breakpoints)

pdf("../dam-ts.pdf", width=4, height=4)
plot(na.approx(bts)/10000, ylab="EVI", ylim=c(0, 1))
dev.off()
#bplm = lm(response ~ (harmon + trend), data=bpp)
#lines(fitted(bplm)/10000~bpp$time, col="blue")
for (i in bf$breakpoints)
{
    print(bpp$time[[i]])
    abline(v=bpp$time[[i]], col="red")
}

poi_val = DeforestationTS[1,]
dates = dates[1:length(DeforestationTS)]
bts = bfastts(poi_val, dates, type="10-day")
bpp = bfastpp(bts, 3)
bf = breakpoints(response ~ (harmon + trend), data=bpp, h=GetBreakNumber(dates))
print(bf$breakpoints)

pdf("../deforestation-ts.pdf", width=4, height=4)
plot(na.approx(bts)/10000, ylab="EVI", ylim=c(0, 1))
dev.off()
#bplm = lm(response ~ (harmon + trend), data=bpp)
#lines(fitted(bplm)/10000~bpp$time, col="blue")
for (i in bf$breakpoints)
{
    print(bpp$time[[i]])
    abline(v=bpp$time[[i]], col="red")
}


