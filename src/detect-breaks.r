# R functions for break detection
# to be moved to a subdirectory when a frontend script is made

library(raster)
library(bfastSpatial)
library(strucchange)
library(foreach)
library(doParallel)
library(tools)
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
GetChunkFilename = function(filename, prefix, length)
{
    file.path(dirname(filename), paste0(prefix, "Chunk_", 1:length, "_", basename(filename)))
}

# foreach-based mc.calc
ForeachCalc = function(input_raster, fx, filename, mem_usage=0.9*1024^3, threads=12, ...)
{
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
        setZ(Chunk, getZ(input_raster))
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
SparkCalc = function(input_raster, fx, filename, mem_usage=0.9*1024^3, threads=12, ...)
{
    ChunkInfo = GetChunkSize(input_raster, mem_usage)
    NumChunks = ChunkInfo["NumChunks"]
    BlockSize = ChunkInfo["BlockSize"]
    
    # Lists of chunk filenames: input and output
    ChunkFilenames = GetChunkFilename(filename, "Input", NumChunks)
    ResultFilenames = GetChunkFilename(filename, "Output", NumChunks)
    Filenames = data.frame(ChunkFilenames, ResultFilenames)
    
    scalc = function(Filename, BlockSize, NumChunks)
    {
        ChunkStart = 1+BlockSize*(i-1)
        ChunkEnd = BlockSize*i
        ChunkExtent = extent(input_raster, r1=ChunkStart, r2=ChunkEnd)
        
        if (file.exists(Filename$ChunkFilenames)) # TODO: think of how to iterate over a data.frame
        {
            print(paste0("Chunk ", i, "/", NumChunks, ": File ", ChunkFilenames[i], " already exists, reusing"))
            Chunk = brick(ChunkFilenames[i])
        } else {
            print(paste0("Chunk ", i, "/", NumChunks, ": cropping to ", ChunkFilenames[i]))
            Chunk = crop(input_raster, ChunkExtent, filename=ChunkFilenames[i], progress="text")
            print(paste0("Chunk ", i, "/", NumChunks, ": cropping complete."))
        }
        setZ(Chunk, getZ(input_raster))
        names(Chunk) = names(input_raster)
        
        print(paste0("Chunk ", i, "/", NumChunks, ": processing to ", ResultFilenames[i]))
        ResultChunk = calc(x=Chunk, fun=fx, filename=ResultFilenames[i], ...)
        print(paste0("Chunk ", i, "/", NumChunks, ": processing complete."))
        print(paste0("unlink(", ChunkFilenames[i], ")"))
        unlink(ChunkFilenames[i])
        
        ResultChunk
    }
    spark.lapply(Filenames, scalc, BlockSize, NumChunks)
    
    b_metrics = gdalUtils::mosaic_rasters(gdalfile=ResultFilenames, dst_dataset=filename,
        output_Raster=TRUE, verbose=TRUE, ot="Int16")
}

# Get last break in a pixel time series
GetLastBreakInTile = function(pixel)
{
    # Check whether we have enough non-NA pixels for running breakpoints.full, without doing preprocessing.
    # The right hand side formula calculates the columns in the bfastpp object.
    if (floor(sum(!is.na(pixel)) * GetBreakNumber(dates)) <= 4+(Order-1)*2 )
        return(NA) # Too many NAs
    
    bfts = bfastts(pixel, dates, type="irregular")
    bpp = bfastpp(bfts, order=Order)
    t0 = as.Date("2014-03-16")
    
    if (sctest(efp(response ~ (harmon + trend), data=bpp, h=GetBreakNumber(dates), type="OLS-MOSUM"))$p.value > 0.05) # If test says there should be no breaks
        return(as.integer(dates[1] - t0))
    
    bf = breakpoints(response ~ (harmon + trend), data=bpp, h=GetBreakNumber(dates))
    
    # Direct returns without calling functions
    if (all(is.na(bf$breakpoints)))
        return(as.integer(dates[1] - t0))
    
    return(as.integer(as.Date(date_decimal(bpp$time[max(bf$breakpoints)])) - t0))
    
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
    
    return(bf)
}

timeseries = LoadTimeSeries("../../landsat78/mosaics-since2013/ndvi")
dates = getZ(timeseries) # This is needed in GetLastBreakInTile, otherwise data is lost; no way to get around using the environment unless we want to re-read names on each pixel process
Order = 3 # Which harmonic order to use
EnableFastBfast()
ForeachCalc(timeseries, GetLastBreakInTile, "../../landsat78/breaks/ndvi/breaks-ndvi-since2013.tif", datatype="INT2S",
    progress="text", options="COMPRESS=DEFLATE")
