#!/bin/bash
# usage: bash run_h2_hunt_round.sh R KLO KHI KEEP MAXS   (round R > 0 of the hunt: seed 4000 + k + 1000 R; times.txt)
cd "$(dirname "$0")"; export MSS_REPO=${MSS_REPO:-../../..}; export H2_ROUND=$1 H2_MAXS=$5
for k in $(seq $2 $3); do
  s=$(date +%s); python3 h2_hunt.py $k $4 --write 2> err_hunt_${k}_r$1.txt; echo "h2_hunt round $1 k=$k keep $4 maxs $5 exit=$? $(( $(date +%s) - s )) s" >> times.txt
done
