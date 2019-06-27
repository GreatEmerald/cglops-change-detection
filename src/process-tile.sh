#!/bin/sh
# Process a tile from start to finish in 3 steps

# Tile to process
VI=$1
Tile=$2
Order=$3

# Step 1: chunk input into ~1400 chunks locally
Rscript ./detect-breaks.R -t $Tile -v $VI --crop-only -m none -o $Order || exit 1

# Step 2: run the processing on the cluster
./spark-submit.sh detect-breaks.R 3 -t $Tile -v $VI -o $Order || exit 1

# Step 3: postprocess the result locally
#Rscript ./postprocess-breaks.r -i /data/users/Public/greatemerald/modis/breaks-bfast-2018/$VI/$Tile/breaks-order${Order}.tif || exit 1
