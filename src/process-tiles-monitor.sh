#!/bin/bash
# Process a list of tiles one by one

# X15 complete, X16 still missing some for Y03?
# Processed X17-X23 for Y03-Y10, continuing with X15-X23 and Y06-Y10
# X19Y07 stopped due to kerberos, continuing from X15Y08
tilesy=({08..10})
for i in ${tilesy[@]}; do
  for n in {15..23}; do
    tile=X${n}Y${i}
    echo $tile
    ./process-tile-monitor.sh NDMI $tile 3 || exit 1
  done
done
