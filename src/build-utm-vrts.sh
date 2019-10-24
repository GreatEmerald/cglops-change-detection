#!/bin/bash
for dir in $1/*
do gdalbuildvrt $dir.vrt $dir/*.tif
done
