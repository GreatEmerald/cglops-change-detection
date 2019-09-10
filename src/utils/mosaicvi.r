library(raster)
library(gdalUtils)

source("utils/dates.r")

# Mosaic a given VI into a VRT. Returns a raster of the mosaic.
MosaicVI = function(VI, InputDir="/data/mep_cg1/MOD_S10/additional_VIs_2009-2018/", VIMosaicFile=paste0("../data/modis-mosaic-", VI, ".vrt"), Reset=FALSE)
{
    if (!Reset && file.exists(VIMosaicFile))
        return(brick(VIMosaicFile))
    
    FilesToMosaic = list.files(InputDir, glob2rx(paste0("*", VI, ".tif")), recursive = TRUE, full.names = TRUE)
    gdalbuildvrt(FilesToMosaic, VIMosaicFile)
    return(brick(VIMosaicFile))
}
