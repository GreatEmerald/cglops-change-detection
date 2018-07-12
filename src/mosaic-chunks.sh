#!/bin/bash
# Script for mosaicking tiles that have been cropped into a single mosaic file, while not running out of RAM
# $1 is VI and $2 is tile
OutputDir="/data/users/Public/greatemerald/modis/breaks/$1/$2"
RAMUse=4000
GDALOptions="-wm ${RAMUse} -ot Int16 -co COMPRESS=DEFLATE -co ZLEVEL=9 -co SPARSE_OK=TRUE"
ChunkStep=2000
ChunkNames=($OutputDir/Output_Chunk_*.tif)
if [[ ${#ChunkNames[@]} -gt ${ChunkStep} ]]; then
    for c in `seq 1 $((${#ChunkNames[@]} / $ChunkStep))`; do
        ChunkStart=$(($c * $ChunkStep))
        ChunksToMosaic=${ChunkNames[@]:$ChunkStart:$ChunkStep}
        VrtName=${OutputDir}/Mosaic_Chunk_$c.vrt
        TiffName=${OutputDir}/Mosaic_Chunk_$c.tif
        gdalbuildvrt ${VrtName} ${ChunksToMosaic} || exit 1
        GDAL_CACHE_MAX=${RAMUse} gdalwarp ${GDALOptions} ${VrtName} ${TiffName} || exit 1
    done
    # Mosaic the mosaics
    GDAL_CACHE_MAX=${RAMUse} gdalwarp ${GDALOptions} ${OutputDir}/Mosaic_Chunk_*.tif ${OutputDir}/breaks-order3-mosaic.tif
else
    gdalbuildvrt ${OutputDir}/breaks-order3-mosaic.vrt ${OutputDir}/Output_Chunk_* || exit 1
    GDAL_CACHE_MAX=${RAMUse} gdalwarp ${GDALOptions} ${OutputDir}/breaks-order3-mosaic.vrt ${OutputDir}/breaks-order3-mosaic.tif || exit 1
fi

# Clean all intermediary files (not input, those are cleaned out by executors)
rm ${OutputDir}/Log_Chunk_* ${OutputDir}/Output_Chunk_*
