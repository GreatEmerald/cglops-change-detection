#!/usr/bin/env Rscript
library(optparse)
library(gdalUtils)
# Script for mosaicking tiles that have been cropped into a single mosaic file, while not running out of RAM


parser = OptionParser()
parser = add_option(parser, c("-t", "--tile"), type="character", default="X16Y06",
                    help="Proba-V tile to process. (Default: %default)", metavar="tile")
parser = add_option(parser, c("-v", "--vegetation-index"), type="character", default="NDMI",
                    help="Vegetation index to process. Case sensitive to input files. (Default: %default)", metavar="VI")
args = parse_args(parser)

OutputDir = file.path("/data/users/Public/greatemerald/modis/breaks", args[["vegetation-index"]], args[["tile"]])
OutputFile = file.path(OutputDir, "breaks-order3.tif")

ChunkNames = list.files(OutputDir, pattern=glob2rx("Output_Chunk*.tif"), full.names=TRUE)

# Mosaic in steps of 100 chunks
if (length(ChunkNames) > 100*100)
    stop("Mosaicking over 10k files is not implemented.")

MosaicNames = NULL
for (Centenary in 1:ceiling(length(ChunkNames) / 100))
{
    StartIndex = (Centenary-1)*100+1
    StopIndex = min(Centenary*100, length(ChunkNames))
    VrtFilename = file.path(OutputDir, paste0("Mosaic_Chunk_", Centenary, ".vrt"))
    gdalbuildvrt(ChunkNames[StartIndex:StopIndex], VrtFilename)
    TiffFilename = file.path(OutputDir, paste0("Mosaic_Chunk_", Centenary, ".tif"))
    gdalwarp(VrtFilename, TiffFilename, wm=6000)
    MosaicNames = c(MosaicNames, TiffFilename)
}

# Mosaic the mosaics
gdalwarp(MosaicNames, OutputFile, wm=6000)
