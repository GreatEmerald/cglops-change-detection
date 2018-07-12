#!/bin/bash
# Process a list of tiles one by one

for i in {18..23}; do
    tile=X${i}Y06
    ./process-tile-2.sh EVI $tile || exit 1
done
