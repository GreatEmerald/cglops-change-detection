library(lubridate)
library(bfast)
library(strucchange)

# Module only containing the bfast0n function

EnableFastBfast = function()
{
    if (exists("set_fast_options"))
    {
        print("Using fast BFAST")
        set_fast_options()
        # Disable bfastts modifications due to issue #2
        options(bfast.use_bfastts_modifications=FALSE)
    } else print("Using reference BFAST, install appelmar/bfast for a speed increase")
}

BFAST0NBreaks = function(pixel, DateStart=2009, DateFrequency=23, DateOffset=8, Order=3, t0 = as.Date("2014-01-01"), NoBreakValue = -9999)
{
    # Utility functions: here so that the scope is correct for SparkR
    GetBreakNumberWhole = function(bfts)
    {
        return(round(frequency(bfts)))
    }
    
    GetBreakNumber = function(dates)
    {
        1/((as.numeric(difftime(max(dates), min(dates), units="weeks")))/52.18)
    }
    
    BreakpointToDateSinceT0 = function(breakpoint_index, bpp, t0, offset=0)
    {
        result = as.integer(as.Date(date_decimal(BreakpointDate(breakpoint_index, bpp))) - t0) + offset
        if (is.numeric(result) && !is.na(result) && !is.nan(result) &&
            (result < as.integer(min(dates) - t0) || result > as.integer(max(dates) - t0)))
        {
            cat("Warning: breakpoint date out of valid range!\n") # -1900 to 1376
            cat(c("Note: breakpoint index: ", breakpoint_index ,"\n"))
            cat(c("Note: calculated days since t0: ", result ,"\n"))
            cat(c("Note: breakpoint date: ", BreakpointDate(breakpoint_index, bpp) ,"\n"))
        }
        return(result)
    }
    
    # DOY of breakpoint; NOTE: problems if the year doesn't match!
    BreakpointToDOY = function(breakpoint_index, bpp)
    {
        result = yday(date_decimal(BreakpointDate(breakpoint_index, bpp)))
        return(result)
    }
    
    # The date of the breakpoint in decimal years
    BreakpointDate = function(breakpoint_index, bpp)
    {
        if (!is.numeric(breakpoint_index) || is.nan(breakpoint_index) || is.na(breakpoint_index) ||
            length(breakpoint_index) > 1 || breakpoint_index <= 0 || breakpoint_index > nrow(bpp))
        {
            cat("Warning: asked to calculate date for an invalid breakpoint index!\n")
            cat(c("Note: The index was:", breakpoint_index, "\n"))
            return(NA)
        }
        return(bpp$time[breakpoint_index])
    }
    
    # Utility: In case we can't calculate anything, return NA values for all years.
    ReturnNAs = function()
    {
        rep(NA, length(Years)*3)
    }
    # Same but if there is no break
    ReturnNoBreak = function()
    {
        rep(NoBreakValue, length(Years)*3)
    }
    
    # Check whether we have enough non-NA pixels for running breakpoints.full, without doing preprocessing.
    # The right hand side formula calculates the columns in the bfastpp object.
    #if (floor(sum(!is.na(pixel)) * GetBreakNumber(dates)) <= 4+(Order-1)*2 )
    #    return(rep(NA, length(Years)*3)) # Too many NAs
    
    # Return NA is we have all NA pixels
    if (all(is.na(pixel)))
        return(ReturnNAs())
    
    # Do not process pixels that have too low VI values (most likely bare soil/desert)
    # WARNING: need to check whether the threshold makes sense for non-EVI
    #if (mean(pixel, na.rm=TRUE) < 500)
    #    return(ReturnNoBreak())
    
    #bfts = bfastts(pixel, dates, type=TSType)
    bfts = ts(pixel, start=DateStart, frequency = DateFrequency)
    
    # Handling dates
    dates = as.Date(date_decimal(as.numeric(time(bfts))))
    DateRange = range(dates)
    Years = lubridate::year(DateRange[1]):lubridate::year(DateRange[2])
    
    # Use integers
    if (GetBreakNumberWhole(bfts) <= 4+(Order-1)*2 || GetBreakNumberWhole(bfts) >= floor(sum(!is.na(pixel))/2))
        return(ReturnNAs()) # Too many NAs
    
    bpp = bfastpp(bfts, order=Order)
    
    testforabreak = sctest(efp(response ~ (harmon + trend), data=bpp, h=GetBreakNumber(dates), type="OLS-MOSUM"))
    if (is.null(testforabreak) || is.null(testforabreak$p.value) || !is.finite(testforabreak$p.value)) {
        cat("Warning: sctest did not return a valid value!\n")
        print(str(testforabreak))
        print(testforabreak)
    } else if (testforabreak$p.value > 0.05) # If test says there should be no breaks
        return(ReturnNoBreak())
    
    bf = tryCatch(breakpoints(response ~ (harmon + trend), data=bpp, h=GetBreakNumberWhole(bfts)),
                  error = function(e){print(e); traceback(e); cat(c("Note: pixel values were: ", pixel, "\n")); return(NULL)})
    
    if (is.null(bf))
    {
        cat("Warning: failed to run breakpoints, returning NA!\n")
        return(ReturnNAs())
    }
    
    # Direct returns without calling functions
    if (all(is.na(bf$breakpoints)))
        return(ReturnNoBreak())
    
    # Make a matrix for the output
    OutMatrix = matrix(NoBreakValue, nrow=length(Years), ncol=3, dimnames=list(Years, c("confint.neg", "breakpoint", "confint.pos")))
    ConfInts = confint(bf)$confint # Get confidence interval
    BreakpointYears = as.integer(sapply(ConfInts[,"breakpoints"], BreakpointDate, bpp)) # Get years at which breakpoints happened
    if (any(duplicated(BreakpointYears))) # Sanity check: should never be true
        cat(c("ERROR: Duplicate breakpoint years! Years:", BreakpointYears, "Dates:", sapply(ConfInts[,"breakpoints"], BreakpointDate, bpp), "Breakpoints:", ConfInts[,"breakpoints"], "\n"))
    BreakpointDays = sapply(ConfInts, BreakpointToDateSinceT0, bpp, t0, DateOffset) # Convert indices to days since t0
    OutMatrix[as.character(BreakpointYears),] = BreakpointDays # Put it into our matrix in the right years
    return(c(t(OutMatrix))) # Flatten matrix
}
