#!/bin/bash
# Process a list of tiles one by one

tilesy=(04 {07..10})
for i in ${tilesy[@]}; do
  for n in {19..23}; do
    tile=X${n}Y${i}
    echo $tile
    ./process-tile.sh EVI $tile 3 || exit 1
  done
done
