#!/bin/bash
# Process a list of tiles one by one

# Missing: X20/21/22/23Y07
tilesy=(07)
for i in ${tilesy[@]}; do
  for n in {20..23}; do
    tile=X${n}Y${i}
    echo $tile
    ./process-tile.sh NDMI $tile 3 || exit 1
  done
done
