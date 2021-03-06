library(strucchange)
library(bfast)
source("../src/utils/enable_fast_bfast.r")
source("../src/bfast-cal/plotting.r")
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
#' @param seasonfreq Multiplier for how many bins to use in a season. 1 means one bin per observation, 0.5 means two observations per bin.
#'
#' @return Fractional years of all detected breaks, or FALSE if none detected,
#' or NA if not enough observations/error in running the function.
MODDetectBreaks = function(InputTS, scrange=c(2009, 2019), scsig=0.05, breaks="LWZ",
                           sctype="OLS-MOSUM", maginterval=0.1, magcomponent="trend",
                           magstat="RMSD", magthreshold=-Inf, coefcomponent="trend",
                           coefthresholds=c(0, 0), plot=FALSE, quiet=FALSE, order=3,
                           formula=response ~ harmon + trend, TargetYears=NULL,
                           seasonfreq=0.5, breaknumthreshold=Inf, altformula=NULL, ...)
{
    # The input should be a single row of a matrix.
    if (!is.ts(InputTS))
        InputTS = GetTS(InputTS) # Convert into a ts object
    Observations = sum(!is.na(InputTS))
    if (!quiet)
        print(paste("Observations for point:", Observations))
    h = ceiling(frequency(InputTS)) # Set h to number of observations per year, i.e. frequency of time series
    if (Observations > h*2) {
        bpp = bfastpp(InputTS, order=order, sbins=seasonfreq) # Preprocess the ts into a data.frame
        # Fix season when it's not an integer
        myseason = as.numeric(as.character(bpp$season)) # Deparse season again
        if (!all(is.na(myseason))) # If we failed to deparse, then we're using a fixed bfastpp already, no need to do anything
            bpp$season = cut(myseason, frequency(InputTS)*seasonfreq, ordered_result = TRUE) # Rebin all

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
            # Are we overfitting? If so, try running with altformula
            if (length(bpOptim$breakpoints) > breaknumthreshold)
            {
                return(MODDetectBreaks(InputTS=InputTS, scrange=scrange, scsig=scsig,
                    breaks=breaks, sctype=sctype, maginterval=maginterval,
                    magcomponent=magcomponent, magstat=magstat, magthreshold=magthreshold,
                    coefcomponent=coefcomponent, coefthresholds = coefthresholds, plot=plot,
                    quiet=quiet, order=order, formula = altformula, TargetYears=TargetYears,
                    seasonfreq=seasonfreq, breaknumthreshold = Inf, ...))
            }
            
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
        if (!quiet)
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

#' 3c) BFAST Monitor break detection function
#' 
#' @param InputTS matrix or ts, input time series
#' @param monitor_years vector of starts of monitoring periods
#' @param monitor_length Length of monitoring period, in years
#' @param cloud_threshold Minimal number of observations to run the algorithm
#' @param TargetYears Time of breaks as reported by reference, for plotting.
#' @param quiet Suppress output
#' @param plot Plot diagnostic plots, one per monitor_years
#' @param ... Additional arguments for bfastmonitor
TestMODMonitor = function(InputTS, monitor_years=2016:2018, monitor_length=1.25,
                          cloud_threshold=42, TargetYears=NULL, quiet=FALSE, plot=FALSE, ...)
{
    # The input should be a single row of a matrix.
    if (!is.ts(InputTS))
        InputTS = GetTS(InputTS) # Convert into a ts object
    
    Observations = sum(!is.na(window(InputTS, end=min(monitor_years)+monitor_length)))
    if (!quiet)
        print(paste("Observations for point:", Observations))
    if (Observations < cloud_threshold)
    {
        if (!quiet)
            print("too cloudy")
        return(NA)
    }
    
    Result = NULL
    for (StartYear in monitor_years)
    {
        # Cut time series to length
        ShortTS = window(InputTS, end=StartYear+monitor_length)
        BM = try(bfastmonitor(ShortTS, StartYear, ...))
        
        if ("try-error" %in% class(BM)) {
            print("An error has occurred:")
            print(BM)
            next # Assume no break
        }
        
        if (plot)
        {
            plot(BM, xlab="Time", ylab="NIRv * 255")
            abline(v=TargetYears, col="blue")
        }
        # It can happen that the break is detected early next year, then we just discard those
        if (!is.na(BM$breakpoint) && BM$breakpoint < StartYear+1)
            Result = c(Result, BM$breakpoint)
    }
    if (length(Result) < 1)
        return(FALSE)
    return(Result)
}

#' Run the original BFAST on the VI time series.
#' @param InputTS A matrix row of vegetation indices or a time series.
#' @param stlplus Whether to use stlplus for decomposition.
TestMODBFAST = function(InputTS, stlplus=TRUE, plot=FALSE, ...)
{
    if (!is.ts(InputTS))
        InputTS = GetTS(InputTS) # Convert into a ts object
    Observations = sum(!is.na(InputTS))
    if (Observations < 10)
        return(NA)
    h = ceiling(frequency(InputTS))
    
    # Try to use stlplus first
    if (stlplus)
    {
        bp = try(bfast(InputTS, decomp="stlplus", ...))
        if ("try-error" %in% class(bp))
        {
            print(paste("Couldn't run stlplus, will retry without. Error was:", bp))
        } else {
            if (plot)
                plot(bp)
            return(GetBFASTOutput(bp))
        }
    }
    
    # Else interpolate and use stl
    InputTS = na.approx(InputTS)
    bp = try(bfast(InputTS, decomp="stl", ...))
    if ("try-error" %in% class(bp))
    {
        print(paste("Couldn't run bfast, error was:", bp))
        return(NA)
    }
    if (plot)
        plot(bp)
    return(GetBFASTOutput(bp))
}

# Util to extract output from an original BFAST model
GetBFASTOutput = function(bp, ...)
{
    Result = NULL
    LastIter = bp$output[[length(bp$output)]]
    if (!bp$nobp$Vt)
        Result = c(Result, time(LastIter$Vt)[!is.na(LastIter$Vt)][LastIter$bp.Vt$breakpoints])
    if (!bp$nobp$Wt)
        Result = c(Result, time(LastIter$Wt)[!is.na(LastIter$Wt)][LastIter$bp.Wt$breakpoints])
    if (is.null(Result)) return(FALSE)
    return(Result)
}

# Util to merge breaks less than 1 (or threshold) year(s) away
# Recursive, control max recursion depth with mergemax
MergeCloseBreaks = function(breakyears, threshold=1, mergemax=Inf)
{
    if (mergemax < 1 || length(breakyears) < 2)
        return(breakyears)
    
    Dists = dist(breakyears)
    
    # Everything is already fine
    if (min(Dists)+1e-8 > threshold)
        return(breakyears)
    
    MergeYears = which(as.matrix(Dists)==min(Dists), arr.ind=TRUE)
    # Replace first value with mean of two
    breakyears[MergeYears[1,1]] = mean(breakyears[MergeYears[1,]])
    # Remove second value
    breakyears = breakyears[-MergeYears[1,2]]
    # Continue removing...
    return(MergeCloseBreaks(breakyears, threshold, mergemax-1))
}

# Function for running sctest with extra additions
MySctest = function(InputTS,
                    formula=response ~ harmon + trend, h=NULL, order=3,
                    scrange=c(2009, 2020), scsig=0.05, sctype="OLS-MOSUM")
{
    # The input should be a single row of a matrix.
    # If not, make it one.
    InputTS = c(GetMatrix(InputTS))
    if (!is.ts(InputTS))
        InputTS = GetTS(InputTS) # Convert into a ts object
    if (is.null(h))
        h = ceiling(frequency(InputTS)) # Yearly
    
    bpp = bfastpp(InputTS, order=order)
    Observations = sum(!is.na(InputTS))
    if (Observations < 42)
    {
        message(paste("Too few observations:", Observations))
        return(NA)
    }
    
    SCp = try(sctest(efp(formula,
                     data=bpp[bpp$time > scrange[1] & bpp$time < scrange[2],],
                     h=h/Observations, sctype=type))$p.value)
    if (class(SCp) == "try-error")
        return(NA)
    return(SCp)
}

# Stand-alone sctest with boolean output
IsThereABreak = function(..., scsig=0.05)
{
    SCp = MySctest(scsig=0.05, ...)
    message(paste("p-value is", SCp))
    SC = SCp > scsig
    if (is.null(SC) || is.na(SC)) return(NA)
    if (SC) return(FALSE) # No break detected, return FALSE
    return(TRUE)
}
