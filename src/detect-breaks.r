# R functions for break detection
# to be moved to a subdirectory when a frontend script is made

library(SparkR)
library(rgdal)
library(raster)

# Load from an external package archive for Spark
.libPaths(c("/data/users/Public/greatemerald/r-packages", .libPaths()))
library(bfast)
library(strucchange)
library(lubridate)

# Load time series: expects two files, `timeseries.vrt` and `layernames.txt` generated from bash
# i.e. `ls *.tif | sort -k 1.11 > layernames.txt; gdalbuildvrt -separate timeseries.vrt $(< layernames.txt)`
# This makes sure that the layer names do not desync from the vrt file (inputs follow the same order)
# Returns a RasterBrick (which is in fact a GDAL-based RasterStack)
LoadTimeSeries = function(input_dir)
{
    timeseries = brick(file.path(input_dir, "timeseries.vrt"))
    names = scan(file.path(input_dir, "layernames.txt"), what="character", sep="\n")
    names(timeseries) = names
    
    # Hack: determine type based on the first characters of the first name
    library(bfastSpatial)
    if (substr(names[1], 1, 3) == "MOD")
    {
        timeseries = setZ(timeseries, getMODISinfo(names)$date)
    } else if (substr(names[1], 1, 1) == "L")
    {
        timeseries = setZ(timeseries, getSceneinfo(names)$date)
    } else
    {
        warning("Could not determine sensor type; recheck the input filenames! They must start with MOD or L.");
    }
    
    return(timeseries)
}

GetDatesFromDir = function(dir)
{
    dirnames = list.dirs(dir, FALSE, FALSE)
    dates = parse_date_time(grep(glob2rx("????????"), dirnames, value=TRUE), "ymd")
    return(as.Date(dates))
}

# Util function: calculate the size of breaks so that the minimum time between them amounts to a year
GetBreakNumber = function(dates)
{
    1/(as.numeric(difftime(max(dates), min(dates), units="weeks"))/52.25)
}

# Util function: if available, enable the use of a faster BFAST
EnableFastBfast = function()
{
    if (exists("set_fast_options"))
    {
        print("Using fast BFAST")
        set_fast_options()
    } else print("Using reference BFAST, install appelmar/bfast for a speed increase")
}

GetChunkSize = function(input_raster, mem_usage=0.9*1024^3, overhead_mult=9)
{
    # Block size: given the target memory usage, calculate how much we need
    bpp = as.integer(substr(dataType(input_raster), 4, 4))
    TotalSize = length(input_raster) * bpp * 2 # Need to keep both input and output in memory
    TotalSize = TotalSize * overhead_mult # Empirical values show that cropping uses 8 times more space(?!), 9 to be safe
    
    NumChunks = ceiling(TotalSize / mem_usage)
    BlockSize = nrow(input_raster) / NumChunks
    return(c(NumChunks = NumChunks, BlockSize=BlockSize))
}

# Utility function to generate filenames for each chunk
GetChunkFilename = function(filename, identifier, length)
{
    file.path(dirname(filename), paste0("Chunk_", 1:length, "_", identifier, "_", basename(filename)))
}

# foreach-based mc.calc
ForeachCalc = function(input_raster, fx, filename, mem_usage=0.9*1024^3, threads=12, ...)
{
    library(foreach)
    library(doParallel)
    library(tools)
    ChunkInfo = GetChunkSize(input_raster, mem_usage)
    NumChunks = ChunkInfo["NumChunks"]
    BlockSize = ChunkInfo["BlockSize"]
    
    # Lists of chunk filenames: input and output
    ChunkFilenames = GetChunkFilename(filename, "Input", NumChunks)
    ResultFilenames = GetChunkFilename(filename, "Output", NumChunks)

    registerDoParallel(cores=threads)
    foreach(i=NumChunks:1, .inorder=FALSE, .verbose=TRUE) %dopar%
    {
        psnice(value = min(threads - 1, 19))
        ChunkStart = 1+BlockSize*(i-1)
        ChunkEnd = BlockSize*i
        ChunkExtent = extent(input_raster, r1=ChunkStart, r2=ChunkEnd)
        
        if (file.exists(ChunkFilenames[i]))
        {
            print(paste0("Chunk ", i, "/", NumChunks, ": File ", ChunkFilenames[i], " already exists, reusing"))
            Chunk = brick(ChunkFilenames[i])
        } else {
            print(paste0("Chunk ", i, "/", NumChunks, ": cropping to ", ChunkFilenames[i]))
            Chunk = crop(input_raster, ChunkExtent, filename=ChunkFilenames[i], progress="text")
            print(paste0("Chunk ", i, "/", NumChunks, ": cropping complete."))
        }
        Chunk = setZ(Chunk, getZ(input_raster))
        names(Chunk) = names(input_raster)
        
        print(paste0("Chunk ", i, "/", NumChunks, ": processing to ", ResultFilenames[i]))
        ResultChunk = calc(x=Chunk, fun=fx, filename=ResultFilenames[i], ...)
        print(paste0("Chunk ", i, "/", NumChunks, ": processing complete."))
        print(paste0("unlink(", ChunkFilenames[i], ")"))
        unlink(ChunkFilenames[i])
        
        ResultChunk
    }
    
    b_metrics = gdalUtils::mosaic_rasters(gdalfile=ResultFilenames, dst_dataset=filename,
        output_Raster=TRUE, verbose=TRUE, ot="Int16")
    print(paste0("unlink(", ResultFilenames, ")"))
}

# SparkR-based mc.calc
SparkCalc = function(input_raster, fx, filename, mem_usage=0.9*1024^3, datatype=NULL, options=NULL)
{
    ChunkInfo = GetChunkSize(input_raster, mem_usage)
    NumChunks = ChunkInfo["NumChunks"]
    BlockSize = ChunkInfo["BlockSize"]
    
    # Lists of chunk filenames: input and output
    ChunkFilenames = GetChunkFilename(filename, "Input", NumChunks)
    ResultFilenames = GetChunkFilename(filename, "Output", NumChunks)
    LogFilenames = GetChunkFilename(filename, "Log", NumChunks)
    
    # The actual function that SparkR runs: crop and process a block
    scalc = function(Index)
    {
        # Set up the log
        LogFile = LogFilenames[Index]
        if (!file.exists(LogFile))
            file.create(LogFile)
        LogFile = file(LogFile, "a")
        sink(LogFile, TRUE)
    
        library(raster)
        # Crop the block
        ChunkStart = 1+BlockSize*(Index-1)
        ChunkEnd = BlockSize*Index
        ChunkExtent = extent(input_raster, r1=ChunkStart, r2=ChunkEnd)
        
        if (file.exists(ChunkFilenames[Index]))
        {
            print(paste0("Chunk ", Index, "/", NumChunks, ": File ", ChunkFilenames[Index], " already exists, reusing"))
            Chunk = brick(ChunkFilenames[Index])
        } else {
            print(paste0("Chunk ", Index, "/", NumChunks, ": cropping to ", ChunkFilenames[Index]))
            Chunk = crop(input_raster, ChunkExtent, filename=ChunkFilenames[Index])
            print(paste0("Chunk ", Index, "/", NumChunks, ": cropping complete."))
        }
        Chunk = setZ(Chunk, getZ(input_raster))
        names(Chunk) = names(input_raster)
        
        # Process the block
        print(paste0("Chunk ", Index, "/", NumChunks, ": processing to ", ResultFilenames[Index]))
        .libPaths(c("/data/users/Public/greatemerald/r-packages", .libPaths()))
        library(strucchange)
        library(bfast)
        library(lubridate)
        if (!is.null(datatype))
        {
            if (!is.null(options))
                ResultChunk = calc(x=Chunk, fun=fx, filename=ResultFilenames[Index], datatype=datatype, options=options)
            else
                ResultChunk = calc(x=Chunk, fun=fx, filename=ResultFilenames[Index], datatype=datatype)
        }
        else
        {
            if (!is.null(options))
                ResultChunk = calc(x=Chunk, fun=fx, filename=ResultFilenames[Index], options=options)
            else
                ResultChunk = calc(x=Chunk, fun=fx, filename=ResultFilenames[Index])
        }
        print(paste0("Chunk ", Index, "/", NumChunks, ": processing complete."))
        print(paste0("unlink(", ChunkFilenames[Index], ")"))
        unlink(ChunkFilenames[Index])
        
        # Stop logging
        sink()
    }
    sparkR.session()
    spark.lapply(1:length(ChunkFilenames), scalc)
    sparkR.session.stop()
    
    b_metrics = gdalUtils::mosaic_rasters(gdalfile=ResultFilenames, dst_dataset=filename,
        output_Raster=TRUE, verbose=TRUE, ot="Int16")
}

# Utility functions
BreakpointToDateSinceT0 = function(breakpoint_index, bpp, t0)
{
    as.integer(as.Date(date_decimal(BreakpointDate(breakpoint_index, bpp))) - t0)
}

# The date of the breakpoint in decimal years
BreakpointDate = function(breakpoint_index, bpp)
{
    bpp$time[breakpoint_index]
}

# Get last break in a pixel time series
GetLastBreakInTile = function(pixel)
{
    # Check whether we have enough non-NA pixels for running breakpoints.full, without doing preprocessing.
    # The right hand side formula calculates the columns in the bfastpp object.
    if (floor(sum(!is.na(pixel)) * GetBreakNumber(dates)) <= 4+(Order-1)*2 )
        return(NA) # Too many NAs
    
    bfts = bfastts(pixel, dates, type=TSType)
    bpp = bfastpp(bfts, order=Order)
    
    if (sctest(efp(response ~ (harmon + trend), data=bpp, h=GetBreakNumber(dates), type="OLS-MOSUM"))$p.value > 0.05) # If test says there should be no breaks
        return(rep(-9999, length(Years)*3))
    
    bf = breakpoints(response ~ (harmon + trend), data=bpp, h=GetBreakNumber(dates))
    
    # Direct returns without calling functions
    if (all(is.na(bf$breakpoints)))
        return(rep(-9999, length(Years)*3))
    
    # Make a matrix for the output
    OutMatrix = matrix(-1, nrow=length(Years), ncol=3, dimnames=list(Years, c("confint.neg", "breakpoint", "confint.pos")))
    ConfInts = confint(bf)$confint # Get confidence interval
    BreakpointYears = as.integer(sapply(ConfInts[,"breakpoints"], BreakpointDate, bpp)) # Get years at which breakpoints happened
    if (any(duplicated(BreakpointYears))) # Sanity check: should never be true
        cat(c("Duplicate breakpoint years!", ConfInts[,"breakpoints"]))
    BreakpointDays = sapply(ConfInts, BreakpointToDateSinceT0, bpp, t0) # Convert indices to days sinec t0
    OutMatrix[rownames(OutMatrix) %in% BreakpointYears,] = BreakpointDays # Put it into our matrix in the right years
    return(c(t(OutMatrix))) # Flatten matrix
    
    #return(as.integer(as.Date(date_decimal(bpp$time[max(bf$breakpoints)])) - t0))
    
    #print("BF:")
    #print(bf)
    
    #if (inherits(bf, "breakpoints"))
    #    return(GetLastBreak(bf))
    
    # For AnalysePixel() compat
    #if (inherits(bf, "bfast"))
    #{
        # No breaks, e.g. 16252501
    #    if (bf$nobp$Vt && bf$nobp$Wt)
    #        return(IdToProbaVEpoch(1))
    #    return(GetLastBreak(bf$output[[length(bf$output)]]$bp.Vt))
    #}
    
    #return(bf)
}

#timeseries = LoadTimeSeries("../../landsat78/mosaics-since2013/ndvi")
#dates = getZ(timeseries) # This is needed in GetLastBreakInTile, otherwise data is lost; no way to get around using the environment unless we want to re-read names on each pixel process
dates = GetDatesFromDir("/data/mep_cg1/MOD_S10/")
timeseries = brick("/data/mep_cg1/MOD_S10/additional_VIs/X16Y06/MOD_S10_TOC_X16Y06_20090101-20171231_250m_C6_EVI.tif")
timeseries = setZ(timeseries, dates)

DateRange = range(dates)
Years = year(DateRange[1]):year(DateRange[2])
t0 = as.Date("2014-03-16")

TSType = "10-day" # Type of time series
Order = 3 # Which harmonic order to use

EnableFastBfast()
#ForeachCalc(timeseries, GetLastBreakInTile, "../../landsat78/breaks/ndvi/breaks-ndvi-since2013.tif", datatype="INT2S",
#    progress="text", options="COMPRESS=DEFLATE")

SparkCalc(timeseries, GetLastBreakInTile, "/data/users/Public/greatemerald/modis/breaks/evi/breaks-X16Y06-order3.tif", datatype="INT2S", options="COMPRESS=DEFLATE")
