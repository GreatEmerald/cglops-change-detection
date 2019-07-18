#!/bin/bash
# Process a list of tiles one by one

vi=NIRV
tilesy=({03..10})
for i in ${tilesy[@]}; do
  for n in {15..23}; do
    tile=X${n}Y${i}
    if [[ ! -d /data/users/Public/greatemerald/modis/breaks-bfast-2018/$vi/$tile ]]; then
      echo $tile
      ./process-tile.sh $vi $tile 3 || exit 1
    fi
  done
done
