#!/bin/bash
# Build VRTs of input data
outdir=/data/users/Public/greatemerald/modis-utm/input-vrt
indir=/data/mep_cg1/MOD_S10/MOD_UTM/deroob
# Input is the index type (e.g. NDMI_16d_Int16)
for zone in $indir/*
do
  zonename=$(basename $zone)
  mkdir -p $outdir/$1
  gdalbuildvrt $outdir/$1/$zonename.vrt $zone/*$1.tif
done
