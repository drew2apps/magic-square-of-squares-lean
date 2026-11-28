#!/bin/bash
# I3 full scan: H1's ranges k = 13..29 in sequence (one core). Paths: G1_DIR, else the repository's roadmap-g/g1.
cd "$(dirname "$0")"
export G1_DIR=${G1_DIR:-$MSS_REPO/experiments/roadmap-g/g1}
for k in ${KS:-$(seq 13 29)}; do
  [ -f i3_scan_k$k.txt ] && continue
  t=$SECONDS; python3 i3_scan.py $k --write 2> >(grep -v 'stack size' > err_k$k.txt); echo "i3_scan k=$k exit=$? $((SECONDS - t)) s" >> times.txt
done
