# Compare breakpoint detection products and make a summary of whether a break is agreed upon or not

# We have multiple years and there may be a shift in the timing of the breaks.
# We consider the break to match if it is within a year

# Some rules:
# 1) Pixels with yearly breaks (>=6) should be removed (NA), they probably are too unstable
# 2) Only consider breakpoints since 2015
# 3) Hansen is a special case, since it reports a single break, and the precision is just a year
# In the end this should be a stratification map, by which we will get validation data;
# most interesting is from places where there is overlap between several products, because it's likely that those are correct.

# Simulated data after preprocessing (PP_*):
# OneBreakA = c(NA, NA, 193, 111, rep(NA, 14), -1343, 1)
# OneBreakB = c(NA, NA, 274, 92, rep(NA, 14), -1262, 1)
# TwoBreaks = c(NA, NA, 193, 61, 325, 163, rep(NA, 12), -836, 2)
# PseudoHansen = 10
# length(OneBreak)
# 
# # Same, but let's say it's 2015
# OneBreakA = c(rep(NA, 12), 193, 111, rep(NA, 4), 500, 1)
# OneBreakB = c(rep(NA, 12), 274, 92, rep(NA, 4), 500, 1)
# TwoBreaks = c(rep(NA, 12), 193, 61, 325, 163, rep(NA, 4), 1000, 2)
# PseudoHansen = 15

# First, implement a simple case: just check whether htere is a break in a particular year.
# Cases:
# F F | No break
# T F | Break in 2015
# F T | Break in 2016
# T T | Break in both
# ---
# 2 2 4 | How much overlap there is, also in total; the max is always Years * Layers
# 2 2 4 | Overlap of Fs
# We store: 2 2 2 2 4 4

# In the case of NA:
# F N
# N F
# N N
# T T
# ---
# 1 1 2
# 1 1 2

# We create layers yearly; for each year:
# - Stack all input layers by building a VRT that contains the needed layers (multiple -bs)
# - If the layer is NA, leave it
# - If the layer is -9999, set it to F, else set to T
# - Count total TRUEs, count total FALSEs
# - Output intermediary file: count of Ts and Fs for a particular year (INT1U)
# Then, load in intermediary files, stack layers, and calculate total Ts and Fs, save final file

library(gdalUtils)
library(raster)

GetAgreement = function(pixels)
{
    if (all(is.na(pixels)))
        return(c(NA, NA))
    ThereIsNoBreak = sum(pixels == -9999, na.rm=TRUE)
    ThereIsAbreak = sum(pixels != -9999, na.rm=TRUE)
    return(c(ThereIsAbreak, ThereIsNoBreak))
}

YearToBand = function(year)
{
    return((year-2009)*3+2)
}

# Directory containing breakpoint files
BaseDir = "/data/users/Public/greatemerald/modis/breaks/"
OutDir = "/data/users/Public/greatemerald/modis/stratification/"
# Which tiles to process; for now X15-23, Y05-06
TileCombinations = expand.grid(paste0("X", 15:23), paste0("Y0", 5:6), stringsAsFactors=FALSE)
TileList = paste0(TileCombinations$Var1, TileCombinations$Var2)
for (Tile in TileList)
{
    # Get list of layers to compare
    FilesToCompare = c(paste0(BaseDir, "EVI/", Tile, "/breaks-order", 2:3, ".tif"),
                       paste0(BaseDir, "NDMI/", Tile, "/breaks-order", 2:3, ".tif"))
    FileExistence = file.exists(FilesToCompare)
    if (!all(FileExistence))
    {
        print(paste("Skipping tile", Tile, "because input file does not exist:"))
        print(FilesToCompare[!FileExistence])
    } else {
        for (Year in 2015:2016)
        {
            # Band 23 is break in 2016, Band 20 is break in 2015
            VRTFile = paste0(OutDir, Tile, "-", Year, ".vrt")
            OutputFile = paste0(OutDir, Tile, "-", Year, ".tif")
            
            if (file.exists(VRTFile) || file.exists(OutputFile))
                print(paste("Output file exists, skipping:", paste0(OutDir, Tile, "-", Year)))
            else
            {
                gdalbuildvrt(FilesToCompare, VRTFile, separate=TRUE)
                
                # Some bug in gdalbuildvrt results in -b just dropping all bands after the number, which is not what we want.
                # Fix the resulting VRT manually:
                system(paste0('sed -i "s|<SourceBand>1</SourceBand>|<SourceBand>', YearToBand(Year) ,'</SourceBand>|" ', VRTFile))
                # Another bug is pixels displaying off by one, due to rounding errors. Something like this fixes it:
                # sed s/\.000.*?"/"/g && sed s/79\.999.*?"/80"/g && sed s/39\.999.*?"/40"/g
                # Read in the raster
                VRTRaster = brick(VRTFile)
                
                print(paste("Postprocessing file:", OutputFile))
                Result = calc(VRTRaster, GetAgreement, filename=OutputFile, datatype="INT1U", options="COMPRESS=DEFLATE")
            }
        }
    }
}

