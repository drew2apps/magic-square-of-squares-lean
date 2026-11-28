#!/bin/bash
# I3: 50,014 new degenerate ratios (seeds 5000 + k, 2942 per range k = 13..29), one core.
cd "$(dirname "$0")"
export G1_DIR=${G1_DIR:-$MSS_REPO/experiments/roadmap-g/g1}
for k in ${KS:-$(seq 13 29)}; do
  [ -f i3_new_k$k.txt ] && continue
  t=$SECONDS; python3 i3_scan.py $k --new --write 2> >(grep -v 'stack size' > err_new_k$k.txt); echo "i3_scan k=$k --new exit=$? $((SECONDS - t)) s" >> times.txt
done
