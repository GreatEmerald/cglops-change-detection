#!/usr/bin/env RScript
# R functions for break detection
# to be moved to a subdirectory when a frontend script is made

print(version) # Make sure the R version matches the packages'

# Load from an external package archive for Spark
.libPaths(c("/data/users/Public/greatemerald/r-packages", .libPaths()))
#print(.libPaths())

library(rgdal)
library(raster)
library(lubridate)
library(optparse)

# Parse command line options
parser = OptionParser()
parser = add_option(parser, c("-t", "--tile"), type="character", default="X16Y06",
                    help="Proba-V tile to process. (Default: %default)", metavar="tile")
parser = add_option(parser, c("-v", "--vegetation-index"), type="character", default="NDMI",
                    help="Vegetation index to process. Case sensitive to input files. (Default: %default)", metavar="VI")
parser = add_option(parser, c("-i", "--input-brick"), type="character",
                    help="Input brick file of MODIS vegetation index time series.", metavar="file")
parser = add_option(parser, c("-c", "--crop-only"), type="logical", action="store_true",
                    help="Run only the step of cropping the input into chunks.")
parser = add_option(parser, c("-l", "--log"), type="logical", action="store_true",
                    help="Output log files next to the input and output chunks. If not specified, everything is output to stdout.")
parser = add_option(parser, c("-m", "--method"), type="character", default="SparkR",
                    help="Method to use for multithreading: SparkR, foreach, none. (Default: %default)", metavar="method")
parser = add_option(parser, c("-y", "--year"), type="integer", default=2014,
                    help="The year since which we are detecting breaks. (Default: %default)", metavar="year")
args = parse_args(parser)

if (args[["method"]] == "SparkR")
    library(SparkR)

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
    file.path(dirname(filename), paste0(identifier, "_Chunk_", 1:length , "_", basename(filename)))
}

# Utility function for faster/more memory efficient cropping
FastCrop = function(input_raster, crop_extent, filename, reference=FALSE, ...)
{
    if (reference)
        return(crop(input_raster, crop_extent, filename=filename, ...))
    else
        return(gdalwarp(input_raster@file@name, filename, te=c(crop_extent@xmin, crop_extent@ymin, crop_extent@xmax, crop_extent@ymax),
                        output_Raster=TRUE))
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
SparkCalc = function(input_raster, fx, filename, mem_usage=0.15*1024^3, datatype=NULL, options=NULL)
{
    ChunkInfo = GetChunkSize(input_raster, mem_usage)
    NumChunks = ChunkInfo["NumChunks"]
    BlockSize = ChunkInfo["BlockSize"]
    
    # Lists of chunk filenames: input and output
    ChunkFilenames = GetChunkFilename(filename, "Input", NumChunks)
    ResultFilenames = GetChunkFilename(filename, "Output", NumChunks)
    LogFilenames = GetChunkFilename(filename, "Log", NumChunks)
    TempResultFilenames = GetChunkFilename(file.path("tmp", basename(filename)), "Output", NumChunks)
    #TempResultGrds = as.character(unlist(as.data.frame(strsplit(TempResultFilenames, "[.]"))[1,]))
    if (!dir.exists(dirname(filename)))
        dir.create(dirname(filename), recursive=TRUE)
    
    # The actual function that SparkR runs: crop and process a block
    scalc = function(Index)
    {
        # Set up the log
        if (!is.null(args[["log"]]))
        {
            LogFile = LogFilenames[Index]
            if (!file.exists(LogFile))
                file.create(LogFile)
            LogFile = file(LogFile, "a")
            sink(LogFile, TRUE)
            sink(LogFile, TRUE, "message") # This sinks stderr: potentially dangerous!
        }
    
        # Disable file checks, because we now check it outside the loop
        # if (file.exists(ResultFilenames[Index]))
        # {
        #     cat(c("Output file ", ResultFilenames[Index], " exists, not recalculating. Delete the file to recalculate."))
        #     sink(type="message")
        #     sink()
        #     return()
        # }
        
        # Set up raster options    
        options(warn=1)
        .libPaths(c("/data/users/Public/greatemerald/r-packages", .libPaths()))
        suppressPackageStartupMessages(library(gdalUtils, quietly=TRUE))
        suppressPackageStartupMessages(library(raster, quietly=TRUE))
        if (!dir.exists("tmp"))
            dir.create("tmp")
        rasterOptions(tmpdir="tmp", chunksize=0.9*1024^3)
        
        
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
            Chunk = FastCrop(input_raster, ChunkExtent, filename=ChunkFilenames[Index])
            print(paste0("Chunk ", Index, "/", NumChunks, ": cropping complete."))
        }
        if (!is.null(args[["crop-only"]]))
            return()
        Chunk = setZ(Chunk, getZ(input_raster))
        names(Chunk) = names(input_raster)
        
        # Process the block
        print(paste0("Chunk ", Index, "/", NumChunks, ": processing to ", ResultFilenames[Index]))
        suppressPackageStartupMessages(library(strucchange, quietly=TRUE))
        suppressPackageStartupMessages(library(bfast, quietly=TRUE))
        suppressPackageStartupMessages(library(lubridate, quietly=TRUE))
        # Unclean tmp dir: make sure to remove in order to not need to overwrite
        if (file.exists(TempResultFilenames[Index]))
            file.remove(TempResultFilenames[Index])
        tryCatch(
        if (!is.null(datatype))
        {
            if (!is.null(options))
                ResultChunk = calc(x=Chunk, fun=fx, filename=TempResultFilenames[Index], datatype=datatype, options=options)
            else
                ResultChunk = calc(x=Chunk, fun=fx, filename=TempResultFilenames[Index], datatype=datatype)
        }
        else
        {
            if (!is.null(options))
                ResultChunk = calc(x=Chunk, fun=fx, filename=TempResultFilenames[Index], options=options)
            else
                ResultChunk = calc(x=Chunk, fun=fx, filename=TempResultFilenames[Index])
        },
        error=function(e){
            print(paste("Error writing output file", TempResultFilenames[Index], e))
            print("List of files in ./tmp:")
            print(list.files("tmp"))
            print("Trying to create a file ./tmp/test:")
            file.create("tmp/test")
            print("List of files in ./tmp:")
            print(list.files("tmp"))
            #print("Trying to write an example raster:")
            #r1 <- raster(nrows=108, ncols=21, xmn=0, xmx=10)
            #print(r1)
            #writeRaster(r1, filename="tmp/test.tif") # Fails
            print("List of files in ./tmp:")
            print(list.files("tmp"))
            print("ResultChunk is:")
            print(ResultChunk)
            print("Trying to write it in GRD:")
            GrdFilename = sub("[.]tif", ".grd", TempResultFilenames[Index])
            writeRaster(ResultChunk, filename=GrdFilename)
            print("List of files in ./tmp:")
            print(list.files("tmp"))
            print("Trying to convert it to a GeoTIFF:")
            ResultChunk = gdal_translate(GrdFilename, TempResultFilenames[Index], ot="Int16", verbose=TRUE, output_Raster=TRUE)
            print("List of files in ./tmp:")
            print(list.files("tmp"))
        })
        print(paste0("Chunk ", Index, "/", NumChunks, ": processing complete."))
        cat(c("Moving file ", TempResultFilenames[Index], " to ", ResultFilenames[Index], "\n"))
        file.copy(TempResultFilenames[Index], ResultFilenames[Index])
        #cat(c("Writing file ", TempResultGrds[Index], " to ", ResultFilenames[Index], "\n"))
        #writeRaster(ResultChunk, filename=ResultFilenames[Index], datatype=datatype)
        cat(c("Cleaning up ", TempResultFilenames[Index], "\n"))
        file.remove(TempResultFilenames[Index])
        print(paste0("unlink(", ChunkFilenames[Index], ")"))
        unlink(ChunkFilenames[Index])
        
        # Stop logging
        if (!is.null(args[["log"]]))
        {
            sink(type="message")
            sink()
        }
    }
    
    # Check which files don't exist and process only those chunks, so that the cluster doesn't need to spawn a VM for doing nothing
    FileIndices = which(!sapply(ResultFilenames, file.exists))
    
    if (args[["method"]] == "SparkR")
    {
        # Workaround for a timeout bug in SparkR 1.4-2.0, fixed in 2.1
        connectBackend.orig <- getFromNamespace('connectBackend', pos='package:SparkR')
        connectBackend.patched <- function(hostname, port, timeout = 3600*48) {
            connectBackend.orig(hostname, port, timeout)
        }
        assignInNamespace("connectBackend", value=connectBackend.patched, pos='package:SparkR')
        
        sparkR.session()
        spark.lapply(FileIndices, scalc)
        sparkR.session.stop()
    } else {
        lapply(FileIndices, scalc)
    }
    
    if (is.null(args[["crop-only"]]))
    {
        print(paste("Starting mosaicking to", filename))
        b_metrics = gdalUtils::mosaic_rasters(gdalfile=ResultFilenames, dst_dataset=filename,
            output_Raster=TRUE, verbose=TRUE, ot="Float32", co="COMPRESS=DEFLATE")
        print(paste("Success, removing chunks."))
        unlink(ResultFilenames)
    }
}

# Get yearly breaks using a t-test
# The schema is different: single layer per year that gives the p-value.
# This gets converted into a binary yearly mask in postprocessing.
TTestBreaks = function(pixel)
{
    endyear = year(today(tzone="Europe/Brussels"))
    Results = numeric(endyear-startyear)
    Results[] = NA
    if (all(is.na(pixel)))
        return(Results)
    
    #t0 = as.Date("2014-01-01")
    
    for (year in startyear:(endyear-1))
    {
        start = as.Date(paste(year, "01", "01", sep="-"))
        end = as.Date(paste(year+1, "01", "01", sep="-"))
        HistoryMean = mean(pixel[dates < start], na.rm=TRUE)
        MonitoringValues = pixel[dates > start & dates < end]
        if (!is.nan(HistoryMean) && sum(!is.na(MonitoringValues)) > 1 && diff(range(MonitoringValues, na.rm=TRUE)) != 0)
        {
            result = tryCatch(t.test(MonitoringValues, mu=HistoryMean, na.rm=TRUE)$p.value,
                  error = function(e){print(e); traceback(e); cat(c("Note: pixel values were: ", pixel, "\n")); return(NA)})
            Results[year+1-startyear] = result
        }
    }
    return(Results)
}

#timeseries = LoadTimeSeries("../../landsat78/mosaics-since2013/ndvi")
#dates = getZ(timeseries) # This is needed in GetLastBreakInTile, otherwise data is lost; no way to get around using the environment unless we want to re-read names on each pixel process
#tile = "X16Y06"
#Vindex = "NDMI"
tile = args[["tile"]]
Vindex = args[["vegetation-index"]]
startyear = args[["year"]]
dates = GetDatesFromDir("/data/mep_cg1/MOD_S10/")
if (!is.null(args[["input-brick"]])) {
    timeseries = brick(args[["input-brick"]])
} else {
    timeseries = brick(paste0("/data/mep_cg1/MOD_S10/additional_VIs_new/", tile, "/MOD_S10_TOC_", tile, "_20090101-20171231_250m_C6_", Vindex, ".tif"))
}
timeseries = setZ(timeseries, dates)

DateRange = range(dates)
Years = lubridate::year(DateRange[1]):lubridate::year(DateRange[2])

TSType = "10-day" # Type of time series

SparkCalc(timeseries, TTestBreaks, file.path("/data/users/Public/greatemerald/modis/breaks", Vindex, tile, "breaks-ttest.tif"), datatype="FLT4S", mem_usage=1024^3)#, options="COMPRESS=DEFLATE")
