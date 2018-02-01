# CGLOPS
library(bfastSpatial)
library(raster)
library(foreach)
library(pryr)
library(strucchange)
library(doParallel)
library(tools)

modis_dir = "/mnt/storage/cglops/modis"

PlotRawData = function(filename=file.path(modis_dir, "mod13q1/MOD13Q1.A2018001.h16v07.006.2018017224110.hdf"))
{
    modis_example = brick(filename)
    sds = get_subdatasets(filename)
    ndvi = raster(sds[1])
    plot(ndvi)
}

# This is too slow, use gdalbuildvrt from bash
MakeTimeBrick = function()
{
    timeStackMODIS(file.path(modis_dir, "tif/"), pattern=glob2rx("*.tif"),
        filename="userdata/cglops/modis/terra-timeseries-sahara.tif", datatype="INT2S", progress="text",
        options="COMPRESS=DEFLATE")
}

LoadTimeSeries = function()
{
    timeseries = brick(file.path(modis_dir, "tif/timeseries.vrt"))
    names = scan(file.path(modis_dir, "tif/layernames.txt"), what="character", sep=" ")
    names(timeseries) = names
    timeseries = setZ(timeseries, getMODISinfo(names)$date)
    return(timeseries)
}

# Gives H as 1/years in ts
GetBreakNumber = function(dates)
{
    1/(as.numeric(difftime(max(dates), min(dates), units="weeks"))/52.25)
}

TestBfast = function()
{
    # Test
    bfm <- bfmPixel(timeseries, start=c(2009, 1), interactive=TRUE)
    # Display a pixel
    pixel = as.vector(timeseries[17204058])
    
    bfts = bfastts(pixel, getZ(timeseries), type = "16-day")
    breaks = bfast(bfts, season="harmonic", max.iter=0xffffffff, h=GetBreakNumber(getZ(timeseries))) # Argument is of length zero if max.iter is not given, weird
    # type=RE/ME/Score-{MO,CU}SUM returns a nicer St but does not change Tt
    plot(breaks)
}

# Function that does bfast analysis on a pixel, returns a bfast
AnalysePixel = function(pixel)
{
    NAs = is.na(pixel)
    if (all(NAs))
    {
        return(NA)
    } else if (any(NAs))
    {
        #warning(paste("A pixel with most but not all NAs!"))
        return(-16380)
    } else
    {
        #print("Entering bfastts")
        #print(mem_used())
        bfts = bfastts(pixel, dates, type = "16-day")
        #print("Entering bfast")
        #print(mem_used())
        # Over 3 iterations seems to be overkill, most take 2
        # bfast is slow, statistics for a single pixel using hte fast implementation:
        # Rec-CUSUM: 2: 2.8, 3: 4.5
        # OLS-CUSUM: 2: 3.0, 3: 4.3 (fewer breaks)
        # Rec-MOSUM: 2: 3.0, 3: 4.5
        # OLS-MOSUM: 2: 2.7, 3: 4.4
        #        RE: 2: 3.5, 3: 6.0
        #        ME: 2: 3.6, 3: 5.6
        # Scr-CUSUM: crash
        # Scr-MOSUM: 2: 3.0, 3: 5.7
        #print(system.time(
        breaks <- bfast(bfts, season="harmonic", max.iter=3,  h=GetBreakNumber(dates))
            #))
        # Argument is of length zero if max.iter is not given, weird
        #print("bfast complete")
        print(mem_used())
        return(breaks)
    }
}

PlotRandomPixel = function()
{
    pixelid = as.integer(runif(1, 0, 4800*4800))
    print(paste("Examining pixel", pixelid))
    pixel = as.vector(timeseries[pixelid])
    breaks = AnalysePixel(pixel)
    if (class(breaks) == "bfast")
        plot(breaks)
    return(breaks)
}

# Breakpoint estimation only, subset of bfast
GetBreakpoints = function(Yt, h=GetBreakNumber(getZ(timeseries)), sctest=TRUE)
{
    ti <- time(Yt)
    St <- stl(Yt, "periodic")$time.series[, "seasonal"]
    Vt <- Yt-St # Deseasoning
    if (sctest)
    {
        p.Vt <- sctest(efp(Vt ~ ti, h=h, type="OLS-MOSUM")) # Is the deseasoned time series constant?
        if (p.Vt$p.value > 0.05)
            print(paste("Warning: p value insignificant:", p.Vt$p.value))
    }
    return(breakpoints(Vt ~ ti, h=h))
}

# Convert an index into days since the first PROBA-V acquisition (2014-03-12)
IdToProbaVEpoch = function(idx)
{
    #print("Image index:")
    #print(idx)
    Absolute = dates[idx]
    #print("Break date:")
    #print(as.integer(Absolute - as.Date("2014-03-12")))
    return(as.integer(Absolute - as.Date("2014-03-12")))
}

# Get last break number
GetLastBreak = function(bp)
{
    #print("Max breakpoint:")
    #print(max(bp$breakpoints))
    if (is.na(bp$breakpoints))
        return(IdToProbaVEpoch(1)) # No breakpoints: use full time series

    return(IdToProbaVEpoch(max(bp$breakpoints)))
}

#bf = PlotRandomPixel()
#GetLastBreak(bf)

GetLastBreakInTile = function(pixel)
{
    if (sum(is.na(pixel)) * GetBreakNumber(dates) >= 8 ) # 8 is the number of columns in a bfastpp formula
        return(NA) # Too many NAs
    
    bfts = bfastts(pixel, dates, type = "16-day")
    # TODO: See if sctest can be ported for this to make it not waste time on things with no breaks
    bf = breakpoints(response ~ (harmon + trend), data=bfastpp(bfts, order=3), h=GetBreakNumber(dates))
    
    # Direct returns without calling functions
    if (all(is.na(bf$breakpoints)))
        return(as.integer(dates[1] - as.Date("2014-03-12")))
    
    return(as.integer(dates[max(bf$breakpoints)] - as.Date("2014-03-12")))
    
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

# Cannot store 17 GiB in memory: need to do chunking

# foreach-based mc.calc
ForeachCalc = function(input_raster, fx, filename, mem_usage=0.9*1024^3, threads=12, ...)
{
    # Block size: given the target memory usage, calculate how much we need
    bpp = as.integer(substr(dataType(input_raster), 4, 4))
    TotalSize = length(input_raster) * bpp * 2 # Need to keep both input and output in memory
    TotalSize = TotalSize * 9 # Empirical values show that cropping uses 8 times more space(?!), 9 to be safe
    print(paste("Chunk total size:", TotalSize / 1024 / 1024, "MiB, currently used:"))
    print(mem_used())
    NumChunks = ceiling(TotalSize / mem_usage)
    BlockSize = nrow(input_raster) / NumChunks
    
    # Lists of chunk filenames: input and output
    ChunkFilenames = file.path(dirname(filename), paste0(1:NumChunks, "_input_", basename(filename)))
    ResultFilenames = file.path(dirname(filename), paste0(1:NumChunks, "_output_", basename(filename)))

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
            print(paste0("Chunk ", i, "/", NumChunks, ": cropping complete, memory used:"))
            print(mem_used())
        }
        #print(object.size(Chunk))
        #rm(Chunk)
        #gc()
        #Chunk = brick(ChunkFilenames[i])
        setZ(Chunk, getZ(input_raster))
        names(Chunk) = names(input_raster)
        #print(object.size(Chunk))
        
        print(paste0("Chunk ", i, "/", NumChunks, ": processing to ", ResultFilenames[i]))
        ResultChunk = calc(x=Chunk, fun=fx, filename=ResultFilenames[i], ...)
        print(paste0("Chunk ", i, "/", NumChunks, ": processing complete, memory used:"))
        print(mem_used())
        
        print(paste0("unlink(", ChunkFilenames[i], ")"))
        
        ResultChunk
    }
    
    #b_metrics <- gdalUtils::mosaic_rasters(gdalfile=ResultFilenames, dst_dataset=filename,
    #    output_Raster=TRUE, verbose=TRUE, ot="Uint16")
    print(paste0("unlink(", ResultFilenames, ")"))
}

# If we are on the appelmar fork, use it
if (exists("set_fast_options"))
{
    print("Using fast BFAST")
    set_fast_options()
} else print("Using reference BFAST, install appelmar/bfast for a speed increase")


timeseries = LoadTimeSeries()
dates = getZ(timeseries) # This is needed in AnalysePixel
ForeachCalc(timeseries, GetLastBreakInTile, file.path(modis_dir, "timeseries/breaks.tif"), datatype="INT2S",
    progress="text", options="COMPRESS=DEFLATE")

BfastBenchmark = function()
{
    # Benchmark the four options
    library(microbenchmark)
    # 15146057: one break at 379
    # 11167368: breaks at 55, 130, 190 from second iteration
    #           Partial returns "33  56  79 107 130 171 196 286 313 378"
    #           Bp returns "44 103 190"
    # 10971385: breaks at 80, 148, 239, 271, 294, 358
    p = as.vector(timeseries[23039998])
    bfts = bfastts(p, dates, type = "16-day")
    
    # full bfast
    GetBfast = function() bfast(bfts, season="harmonic", max.iter=2,  h=GetBreakNumber(dates))$output[[2]]$bp.Vt$breakpoints
    
    # only breakpoints
    GetBp = function()
    {
        bpp = bfastpp(bfts, order=3)
        return(breakpoints(response ~ (harmon + trend), data=bpp, h=GetBreakNumber(dates))$breakpoints)
    }
    
    GetBfastPartial = function() GetBreakpoints(bfts, sctest=TRUE)$breakpoints
    GetBfastPartialAlways = function() GetBreakpoints(bfts, sctest=FALSE)$breakpoints
    
    my_check <- function(values) {
        print(values)
      all(sapply(values[-1], function(x) identical(values[[1]], x)))
    }
    mb = microbenchmark(GetBfastPartial(), GetBfastPartialAlways(), GetBfast(), GetBp() )#, check=my_check)
    print(mb)
    
    # On reference bfast, GetBp() is the slowest, then GetBfast(), then GetBfastPartial{,Always}()
    # On fast bfast, GetBp() is the fastest (twice/thrice!), then GetBfastPartial{,Always}(), then GetBfast()
    # So use Partial on slow and GetBp on fast version
    # The workstation is 3 times faster than the VM!
}

#BfastBenchmark()
