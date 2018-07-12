#!/usr/bin/env Rscript
# Postprocess a breakpoint output file

library(optparse)
library(raster)
library(lubridate)

# Input: 3 layers per year: confidence interval start, expected value of breakpoint, confidence interval end
# Output:
#   - Two layers per year: 1) Day of year of expected value; 2) Length of the confidence interval in days;
#     plus number of breakpoints and the last breakpoint. -9999 replaced with NA.
#   - In the future, potentially a binary mask for each time step indicating whether there was a (potential) break or not.

parser = OptionParser()
parser = add_option(parser, c("-o", "--output"), type="character",
                    help="Output filename. (Default: input filename prefixed with 'PP_')", metavar="file")
parser = add_option(parser, c("-i", "--input"), type="character",
                    help="Input file containing breakpoint detection output.", metavar="file")
args = parse_args(parser)

Input = brick(args[["input"]])
if (!is.null(args[["output"]])) {
    OutputFile = args[["output"]]
} else {
    OutputFile = file.path(dirname(args[["input"]]), paste0("PP_", basename(args[["input"]])))
    print(paste("Will write output to", OutputFile))
}

# Raise memory usage to avoid I/O; NB: these are in cells
rasterOptions(maxmemory=2*1024^3/16, chunksize=0.5*1024^3/16)

VisualPostprocessor = function(pixels)
{
    DayOfYear = function(absolute_date, t0 = as.Date("2014-03-16"))
    {
        AbsoluteDate = as.Date(absolute_date, origin=t0)
        return(yday(AbsoluteDate))
    }
    
    Result = matrix(pixels, nrow=3) # Columns are years
    Result[Result == -9999] = NA # Filter out -9999
    ExpectedValues = Result[2,]
    BreakpointCount = sum(!is.na(ExpectedValues))
    MaxValue = ifelse(all(is.na(ExpectedValues)), NA, max(ExpectedValues, na.rm=TRUE)) # Max breakpoint
    IntervalSum = function(YearValues)
    {
        ExpectedValue = DayOfYear(YearValues[2]) # Convert to days of the year
        IntervalLength = YearValues[3] - YearValues[1]
        return(c(ExpectedValue, IntervalLength))
    }
    Result = apply(Result, 2, IntervalSum)
    return(c(Result, MaxValue, BreakpointCount))
    #return(c(Result))
}

calc(Input, VisualPostprocessor, filename=OutputFile, options="COMPRESS=DEFLATE", progress="text", datatype="INT2S")
