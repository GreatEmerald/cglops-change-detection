#!/bin/bash
# Process a list of tiles one by one

# Stopped at X17Y03
tilesy=({05..10})
for i in ${tilesy[@]}; do
  for n in {15..23}; do
    tile=X${n}Y${i}
    echo $tile
    ./process-tile.sh EVI $tile 3 || exit 1
  done
done
