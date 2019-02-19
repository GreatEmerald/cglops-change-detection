#!/bin/bash
# Process a list of tiles one by one

tilesy=({03..10})
for i in ${tilesy[@]}; do
  for n in {15..23}; do
    tile=X${n}Y${i}
    echo $tile
    ./process-tile-monitor.sh NDMI $tile 3 || exit 1
  done
done
