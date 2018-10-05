# Convert Hansen's data into the same format as our own generated breakpoints
library(raster)
library(gdalUtils)
library(foreach)
library(doParallel)

# Download the data
StorageDir = "/home/greatemerald/shared/hansen"

# Tiles to download
TileCombinations = expand.grid(paste0(1:3, "0N"), c(paste0("0", 1:3, "0W"), paste0("0", 0:5, "0E")), stringsAsFactors=FALSE)
TileList = paste(TileCombinations$Var1, TileCombinations$Var2, sep="_")
for (Tile in TileList)
{
    OutputFile = file.path(StorageDir, paste0(Tile, ".tif"))
    if (!file.exists(OutputFile))
        download.file(paste0("https://storage.googleapis.com/earthenginepartners-hansen/GFC-2017-v1.5/Hansen_GFC-2017-v1.5_lossyear_", Tile, ".tif"),
                  OutputFile)
}
DownloadedFiles = list.files(StorageDir, pattern=glob2rx("?0N_0?0?.tif"), full.names=TRUE)
HansenSourceFile = file.path(StorageDir, "HansenMosaic.vrt")
if (!file.exists(HansenSourceFile))
    gdalbuildvrt(DownloadedFiles, HansenSourceFile)
HansenSource = raster(HansenSourceFile)
#spplot(HansenSource)

PVTileCombinations = expand.grid(paste0("X", 15:23), paste0("Y0", 5:6), stringsAsFactors=FALSE)
PVTileList = paste0(PVTileCombinations$Var1, PVTileCombinations$Var2)

registerDoParallel(cores = 4)
foreach (Tile=iter(PVTileList), .inorder=FALSE) %dopar%
{
    SampleRaster = raster(paste0("/data/mep_cg1/MOD_S10/20090101/MOD_S10_TOC_", Tile, "_20090101_250M_C6.tiff"))
    OutputFile = file.path(StorageDir, "ProbaTiles", paste0(Tile, ".tif"))
    if (!file.exists(OutputFile))
    {
        CroppedFile = file.path(StorageDir, "Cropped", paste0(Tile, ".tif"))
        CroppedTile = crop(HansenSource, SampleRaster, snap="out",
                 progress="text", filename=CroppedFile, datatype="INT1U", options="COMPRESS=DEFLATE")
        resample(CroppedTile, SampleRaster, method="ngb",
                 progress="text", filename=OutputFile, datatype="INT1U", options="COMPRESS=DEFLATE")
    }
}
