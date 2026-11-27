#!/bin/bash
# usage: bash run_h2_hunt.sh KLO KHI KEEP   (ranges k = KLO..KHI in turn; times.txt)
cd "$(dirname "$0")"; export MSS_REPO=${MSS_REPO:-../../..}
for k in $(seq $1 $2); do
  s=$(date +%s); python3 h2_hunt.py $k $3 --write 2> err_hunt_$k.txt; echo "h2_hunt k=$k keep $3 exit=$? $(( $(date +%s) - s )) s" >> times.txt
done
