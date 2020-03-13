# Break detection functions

#' 3) Run BFAST over the time series and get the detected breaks.
#'
#' This is the function that should be optimised;
#' optional arguments should include all the tunables.
#' Ideally we would run this function in an optimiser (?optim),
#' so it tries all combinations of the variables and finds the best one.
#'
#' @param scsig significance value at which the structure change test is considered successful
#' @param scrange range over which the sctest should be run
#' @param sctype type of sctest (?efp)
#' @param maginterval interval (% time series) over which to compute breakpoint magnitudes
#' @param magcomponent components (possibly multiple for which to calculate the magnitudes (trend/harmonsin1/harmoncos3/...)
#' @param magstat statistic for magnitude thresholding (diff/RMSD/MAD/MD)
#' @param magthreshold threshold above which the breaks are kept. Breaks lower than the threshold are discarded
#' @param coefcomponent component (one!) for coefficient difference calculation
#' @param coefthresholds min and max, if coefficient difference is within the range, the breakpoint will be discarded.
#' @param plot Whether to call plot.bfast0n() on the output
#' @param quiet Suppress print output
#' @param order Harmonic order
#' @param formula Formula passed to sctest() and bfastpp()
#' @param TargetYears Year fractions of expected breaks for plotting
#
# The output is fractional years of all detected breaks, or FALSE if none detected,
# or NA if not enough observations/error in running the function.
MODDetectBreaks = function(InputTS, scrange=c(2009, 2019), scsig=0.05, breaks="LWZ",
                           sctype="OLS-MOSUM", maginterval=0.1, magcomponent="trend",
                           magstat="RMSD", magthreshold=-Inf, coefcomponent="trend",
                           coefthresholds=c(0, 0), plot=FALSE, quiet=FALSE, order=3,
                           formula=response ~ harmon + trend, TargetYears=NULL, ...)
{
    # The input should be a single row of a matrix.
    InputTS = GetTS(InputTS) # Convert into a ts object
    Observations = sum(!is.na(InputTS))
    if (!quiet)
        print(paste("Observations for point:", Observations))
    h = frequency(InputTS) # Set h to number of observations per year, i.e. frequency of time series
    if (Observations > h*2) {
        bpp = bfastpp(InputTS, order=order) # Preprocess the ts into a data.frame
        # Run the sctest first to determine whether to run bfast0n
        if (!is.null(scrange)) {
            SC = sctest(efp(formula,
                            data=bpp[bpp$time > scrange[1] & bpp$time < scrange[2],],
                            h=h/Observations, sctype=type))$p.value > scsig
            if (!is.null(SC) && !is.na(SC) && SC) return(FALSE) # No break detected, return FALSE
        }
        
        # Break was detected, so run bfast0n
        bp = try(breakpoints(formula, data=bpp, h=h))
        
        if ("try-error" %in% class(bp)) {
            print("An error has occurred:")
            print(bp)
            return(NA)
        }
        
        bpOptim = breakpoints(bp, breaks=breaks) # Get breakpoint time
        
        if (length(bpOptim$breakpoints) > 0 && !all(is.na(bpOptim$breakpoints))) {
            # Get magnitudes and coefficients of each break
            bpMag = magnitude(bp, breaks=breaks, interval=maginterval, component=magcomponent)$Mag
            coefs = coef(bp, breaks=breaks)
            bpCoef = coefs[2:nrow(coefs),coefcomponent] - coefs[1:nrow(coefs)-1,coefcomponent]
            names(bpCoef) = NULL
        } else {
            bpMag = NULL
            bpCoef = NULL
        }
        
        if (plot)
        {
            plot.bfast0n(bp, bpp, breaks=breaks, bpMag=bpMag, ...)
            abline(v=TargetYears, col="red")
        }
        
        if (all(is.na(bpOptim$breakpoints))) # bfast0n didn't find any breaks, return FALSE
            return(FALSE)
        if (!quiet)
        {
            print("Breakpoints before filtering:")
            print(bpp[bpOptim$breakpoints, "time"])
            print("Magnitudes of each breakpoint:")
            print(bpMag)
            print(paste("Differences in coefficients of", coefcomponent))
            print(bpCoef)
        }
        
        # Convert the output into dates of break
        Result = bpp[bpOptim$breakpoints, "time"]
        
        if (!is.null(bpMag) && !is.null(bpCoef))
        {
            # Keep breakpoints that are above the threshold of magnitude
            MagFilter = abs(bpMag[,magstat]) > magthreshold
            # Keep breakpoints where the coefficient difference is big enough
            CoefFilter = bpCoef < min(coefthresholds) | bpCoef > max(coefthresholds)
            Result = Result[MagFilter & CoefFilter]
        }
        if (length(Result) < 1)
            return(FALSE) # If we filtered out all results, return FALSE again
        if (plot) # For ones that got kept, plot black on top
            abline(v=Result, col="black")
        return(Result)
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
