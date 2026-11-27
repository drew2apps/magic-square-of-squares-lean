#!/bin/bash
# Ranges stopped by hand at the end of the run (27 Nov, about 10:04 UTC): each is rerun from its .part file with its
# sample limit (H2_MAXS) or take count set to what was done, so it replays the recorded ratios and writes its output.
# usage: bash run_h2_finalize.sh   (times.txt)
cd "$(dirname "$0")"; export MSS_REPO=${MSS_REPO:-../../..}
fin() {  # k keep maxs
  s=$(date +%s); H2_MAXS=$3 python3 h2_hunt.py $1 $2 --write 2>> err_hunt_$1.txt; echo "h2_hunt k=$1 keep $2 maxs $3 (finalized from .part) exit=$? $(( $(date +%s) - s )) s" >> times.txt
}
fin 25 120 1040
fin 30 15 265
fin 33 10 248
fin 36 10 94
s=$(date +%s); python3 h2_design.py 28 249 --write 2>> err_design_28.txt; echo "h2_design k=28 take 249 (finalized from .part) exit=$? $(( $(date +%s) - s )) s" >> times.txt
