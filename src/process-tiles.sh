#!/bin/bash
# Process a list of tiles one by one

for i in {20..23}; do
    tile=X${i}Y05
    ./process-tile.sh EVI $tile 2 || exit 1
done
