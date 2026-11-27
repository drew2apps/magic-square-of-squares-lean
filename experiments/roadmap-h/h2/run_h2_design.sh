#!/bin/bash
# usage: bash run_h2_design.sh KLO KHI TAKE   (ranges k = KLO..KHI in turn; times.txt)
cd "$(dirname "$0")"; export MSS_REPO=${MSS_REPO:-../../..}
for k in $(seq $1 $2); do
  s=$(date +%s); python3 h2_design.py $k $3 --write 2> err_design_$k.txt; echo "h2_design k=$k take $3 exit=$? $(( $(date +%s) - s )) s" >> times.txt
done
