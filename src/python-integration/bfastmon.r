library(lubridate)
library(bfast)

BFMBreaks = function(pixel, MonitorStart=year(today(tzone="Europe/Brussels"))-1,
                     DateStart=2009, DateFrequency=23, NoBreakValue = -9999,
                     formula=response ~ trend, level=0.001, history="all", ...)
{
    Results = NA
    if (all(is.na(pixel)))
        return(Results)
    
    # This is optimised for the dates that we have
    PixelTS = ts(pixel, start=DateStart, frequency=DateFrequency)
    result = tryCatch(bfastmonitor(PixelTS, year, formula=formula, level=level,
                                   history=history, ...),
        error = function(e){print(e); traceback(e); cat(c("Note: pixel values were: ", toString(pixel), "\n")); return(NA)})
    if (all(is.na(result)) || is.null(result[["breakpoint"]])) {
        Results = NA
    } else {
        BreakDay = as.integer(difftime(date_decimal(result$breakpoint),
                            date_decimal(MonitorStart), units="days")) # Fractional year, so convert to day of year
        Results = ifelse(!is.na(BreakDay), BreakDay, NoBreakValue)
    }
    return(Results)
}
