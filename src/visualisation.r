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
