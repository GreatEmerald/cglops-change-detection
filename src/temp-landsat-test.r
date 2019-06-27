# Test Landsat + BFAST
library(raster)
library(bfast)
library(strucchange)

rasterfile = "/data/mep_cg1/MOD_S10/LSAT_SENEGAL/BY_BAND/monthly_composite_2009-01-01_to_2017-12-31_EVI_Senegal_extended_BAND.vrt"
gi = gdalinfo(rasterfile)
EVIBrick = brick(rasterfile)
dates = seq.Date(as.Date("2009-01-01"), as.Date("2017-12-01"), along.with=1:nlayers(EVIBrick))
plot(ts(1:nlayers(EVIBrick), start=c(2009, 1), frequency=12))
Vals = getValuesBlock(EVIBrick, 10000, 1, 10000, 1)
Vals = c(Vals)
Vals[Vals==0] = NA
SampleTS = ts(Vals, start=c(2009, 1), frequency=12)
plot(SampleTS)
SampleBFM = bfastmonitor(SampleTS, start = 2017)
plot(SampleBFM)
SampleBPP = bfastpp(SampleTS, order = 3)
SampleBP = breakpoints(response ~ (harmon + trend), data=SampleBPP, h = 12)
plot(SampleTS)
abline(v=SampleBPP$time[SampleBP$breakpoints], col="red")

# Can we read in a block without loading everything in memory?
ncol(EVIBrick) = 1
EVIBrick

# gdalUtils
library(gdalUtils)
TopLeftExtent = c(xmin=extent(EVIBrick)@xmin,ymin=extent(EVIBrick)@ymin,xmax=extent(EVIBrick)@xmin+30*128,ymax=extent(EVIBrick)@ymin+30*128)
gdalbuildvrt(gdalfile = rasterfile, 
             output.vrt = "../data/tmp.vrt", 
             te = TopLeftExtent)
MyBlock = brick("../data/tmp.vrt")
plot(MyBlock[[3]])
inMemory(MyBlock)
system.time(MyBlock <- readAll(MyBlock))
inMemory(MyBlock)
pryr::object_size(MyBlock) # 7.1 MiB; 4*128*128*108/1024/1024, i.e. it's always FLT4S + overhead of 0.35 MiB
MyBlock[MyBlock == 0] = NA
plot(MyBlock[[3]])

# Try two chunks
TopTwoExtent = c(xmin=extent(EVIBrick)@xmin,ymin=extent(EVIBrick)@ymin,xmax=extent(EVIBrick)@xmin+2*30*128,ymax=extent(EVIBrick)@ymin+30*128)
gdalbuildvrt(gdalfile = rasterfile, output.vrt = "../data/tmp2.vrt", te = TopTwoExtent) # immediate
MyBlock2 = brick("../data/tmp2.vrt") # immediate
system.time(plot(MyBlock2[[3]])) # 0.57s
system.time(MyBlock2 <- readAll(MyBlock2)) # 1.69s - equivalent to one block
inMemory(MyBlock2)
pryr::object_size(MyBlock2) # 14.2 MiB, literally double, so overhead multiplies as well

# One row, or just about, to fill 1G
floor(1024/7.1) # 144 chunks
RowXMax = extent(EVIBrick)@xmin+144*30*128
RowXMax < extent(EVIBrick)@xmax # Not out of bounds yet, but actually fairly close, could try whole rows
TopRowExtent = c(xmin=extent(EVIBrick)@xmin,ymin=extent(EVIBrick)@ymin,xmax=RowXMax,ymax=extent(EVIBrick)@ymin+30*128)
system.time(gdalbuildvrt(gdalfile = rasterfile, output.vrt = "../data/tmp3.vrt", te = TopRowExtent)) # 0.6s
MyBlockRow = brick("../data/tmp3.vrt") # 2359296 cells
ncell(EVIBrick)/ncell(MyBlockRow) # 220 tiles
system.time(plot(MyBlockRow[[3]])) # 9s
system.time(MyBlockRow <- readAll(MyBlockRow)) # 31s - super slow
pryr::object_size(MyBlockRow) # 1.02 GiB as expected
rm(MyBlockRow)
gc()

# Blocks by pixel are 256x256
floor(1024/(7.1*4)) # 36 chunks
RowXMax = extent(EVIBrick)@xmin+36*30*256
RowXMax < extent(EVIBrick)@xmax # Almost exactly half!
TopRowExtent = c(xmin=extent(EVIBrick)@xmin,ymin=extent(EVIBrick)@ymin,xmax=RowXMax,ymax=extent(EVIBrick)@ymin+30*256)
system.time(gdalbuildvrt(gdalfile = rasterfile, output.vrt = "../data/tmp4.vrt", te = TopRowExtent)) # 0.6s
MyBlockRow = brick("../data/tmp4.vrt") # 2359296 cells
ncell(EVIBrick)/ncell(MyBlockRow) # 220 tiles
system.time(plot(MyBlockRow[[3]])) # 6.7s, very slow
system.time(MyBlockRow <- readAll(MyBlockRow)) # 32s - just as slow as before
pryr::object_size(MyBlockRow) # 1.02 GiB as expected
rm(MyBlockRow)
gc()

# What if we try a line
TopRowExtent = c(xmin=extent(EVIBrick)@xmin,ymin=extent(EVIBrick)@ymin,xmax=extent(EVIBrick)@xmax,ymax=extent(EVIBrick)@ymin+30*128)
system.time(gdalbuildvrt(gdalfile = rasterfile, output.vrt = "../data/tmp5.vrt", te = TopRowExtent)) # 0.6s
MyBlockRow = brick("../data/tmp5.vrt") # 3110784 cells
ncell(EVIBrick)/ncell(MyBlockRow) # 166 tiles
system.time(plot(MyBlockRow[[3]])) # 6s
system.time(MyBlockRow <- readAll(MyBlockRow)) # DIES
pryr::object_size(MyBlockRow) # 1.02 GiB as expected
rm(MyBlockRow)
gc()

# the above was by band, but it's probably more efficient if it's tiled by pixel. Build a VRT
PixelFiles = list.files("/data/mep_cg1/MOD_S10/LSAT_SENEGAL/BY_PIXEL/", pattern=glob2rx("*.tif"), full.names = TRUE)
gdalbuildvrt(gdalfile = PixelFiles, output.vrt = "../data/tmp-landsat-pixel.vrt")
PixelFile = "../data/tmp-landsat-pixel.vrt"

# Try by pixel
floor(1024/(7.1*4)) # 36 chunks
RowXMax = extent(EVIBrick)@xmin+36*30*256
RowXMax < extent(EVIBrick)@xmax # Almost exactly half!
TopRowExtent = c(xmin=extent(EVIBrick)@xmin,ymin=extent(EVIBrick)@ymin,xmax=RowXMax,ymax=extent(EVIBrick)@ymin+30*256)
system.time(gdalbuildvrt(gdalfile = PixelFile, output.vrt = "../data/tmp6.vrt", te = TopRowExtent)) # 0.6s
MyBlockRow = brick("../data/tmp6.vrt") # 2359296 cells
ncell(EVIBrick)/ncell(MyBlockRow) # 220 tiles
system.time(plot(MyBlockRow[[3]])) # 50s, very slow
system.time(MyBlockRow <- readAll(MyBlockRow)) # 5144s - ridiculously slow
pryr::object_size(MyBlockRow) # 1.02 GiB as expected
rm(MyBlockRow)
gc()

# Looks like VRTs report incorrect block sizes, as they don't correspnd with raw TIFFs. Let's work on a single TIFF instead.
SingleBand = brick("/data/mep_cg1/MOD_S10/LSAT_SENEGAL/BY_BAND/monthly_composite_2009-01-01_to_2017-12-31_NDVI_Senegal_extended-0000018432-0000018432_BAND.tif")
ncell(SingleBand) * 4 * 108 / 1024 / 1024 / 1024 # 5 GiB

# Blocks are 4608x1, i.e. row-wise
TopRowExtent = c(xmin=extent(SingleBand)@xmin,ymin=extent(SingleBand)@ymin,xmax=extent(SingleBand)@xmax,ymax=extent(SingleBand)@ymin+30*256*2)
system.time(gdalbuildvrt(gdalfile = SingleBand@file@name, output.vrt = "../data/tmp7.vrt", te = TopRowExtent)) # 0.6s
MyBlockRow = brick("../data/tmp7.vrt") # 2359296 cells
ncell(SingleBand)/ncell(MyBlockRow) # 6 tiles
system.time(plot(MyBlockRow[[3]])) # 7s
system.time(MyBlockRow <- readAll(MyBlockRow)) # 30s; slightly faster but not by much
pryr::object_size(MyBlockRow) # 1.02 GiB as expected
rm(MyBlockRow)
gc()

# And by pixel is 256x256
SinglePixel = brick("/data/mep_cg1/MOD_S10/LSAT_SENEGAL/BY_PIXEL/monthly_composite_2009-01-01_to_2017-12-31_NDVI_Senegal_extended-0000018432-0000018432.tif")
TopRowExtent = c(xmin=extent(SinglePixel)@xmin,ymin=extent(SinglePixel)@ymin,xmax=extent(SinglePixel)@xmax,ymax=extent(SinglePixel)@ymin+30*256*2)
system.time(gdalbuildvrt(gdalfile = SinglePixel@file@name, output.vrt = "../data/tmp8.vrt", te = TopRowExtent)) # 0.6s
MyBlockRow = brick("../data/tmp8.vrt") # 1179648 cells
ncell(SingleBand)/ncell(MyBlockRow) # 6 tiles
system.time(plot(MyBlockRow[[3]])) # 18s
system.time(MyBlockRow <- readAll(MyBlockRow)) # 518s; this is very slow
pryr::object_size(MyBlockRow) # 1.02 GiB as expected
rm(MyBlockRow)
gc()

# stars
library(stars)

# By band
system.time(MyBlockRow <- read_stars(SingleBand@file@name, RasterIO = list(nYSize=256))) # 62s, so twice as slow
pryr::object_size(MyBlockRow) # 4 MiB for a single row, so 1 GiB for 256 rows; 200s if we scale linearly
system.time(plot(MyBlockRow, rgb=c(3,3,3))) # 8s
rm(MyBlockRow); gc()

# by pixel
system.time(MyBlockRow <- read_stars(SinglePixel@file@name, RasterIO = list(nYSize=256))) # 64s, consistent
pryr::object_size(MyBlockRow) # 4 MiB for a single row, so 1 GiB for 256 rows; 200s if we scale linearly
system.time(plot(MyBlockRow, rgb=c(3,3,3))) # 8s, also consistent

# Raster by band wins:
MyBlockRow = brick("../data/tmp7.vrt")
MyBlockRow[[3]] # But is it valid?
plot(MyBlockRow[[3]])
TestPoint = st_point(c(-12.490501,12.290171))
TP = st_sfc(TestPoint, crs=4326)
TP = st_set_geometry(data.frame(type="point"), TP)
MyVal = extract(MyBlockRow, TP) # Values check out, looks the same on QGIS too; be careful to make sure it's NDVI
MyVal[MyVal==0] = NA
MyVal = c(MyVal) # unmatrix
MyTS = ts(MyVal, start = c(2009, 1), frequency = 12)
plot(MyTS)
bpp = bfastpp(MyTS)
bp = breakpoints(response ~ (harmon + trend), data=bpp, h=12)

# Try a random pixel
RandomPixel = round(runif(1, 0, ncell(MyBlockRow)))
MyVal = extract(MyBlockRow, RandomPixel)
MyVal[MyVal==0] = NA
MyVal = c(MyVal) # unmatrix
MyTS = ts(MyVal, start = c(2009, 1), frequency = 12)
bpp = bfastpp(MyTS)
bp = breakpoints(response ~ (harmon + trend), data=bpp, h=12)
bp
plot(MyTS)
for (i in 1:length(bp$breakpoints))
    abline(v=bpp[bp$breakpoints[i],"time"], col="red")

# Try imputing missing values manually
library(Amelia)
MyTSLong = data.frame(Y=as.matrix(MyTS), date=time(MyTS))
a = amelia(MyTSLong, ts="date", polytime=3)
plot(a$imputations$imp1$Y, col="red", type="l")
plot(a$imputations$imp1$Y, col="red", type="l")
lines(a$imputations$imp2$Y, col="blue")
lines(a$imputations$imp3$Y, col="green")
lines(a$imputations$imp4$Y, col="orange")
lines(a$imputations$imp5$Y, col="purple")
lines(MyTSLong$Y) # Not great
