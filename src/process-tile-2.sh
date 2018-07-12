#!/bin/sh
# Process a tile from start to finish in 3 steps

# Tile to process
VI=$1
Tile=$2

# Step 1: chunk input into ~1400 chunks locally
Rscript ./detect-breaks.R -t $Tile -v $VI --crop-only -m none -o 2 || exit 1

# Step 2: run the processing on the cluster
./spark-submit.sh detect-breaks.R -t $Tile -v $VI -o 2 || exit 1

# Step 3: postprocess the result locally
Rscript ./postprocess-breaks.r -i /data/users/Public/greatemerald/modis/breaks/$VI/$Tile/breaks-order2.tif || exit 1
