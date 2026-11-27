#!/bin/bash
# usage: bash run_h1_checks.sh   (reruns without --write and records whether each output matches; times.txt)
cd "$(dirname "$0")"; export MSS_REPO=${MSS_REPO:-../../..}
for c in "h1_pairs.py --negative" "h1_scan.py --negative" "h1_rank0_exact.py --negative"; do
  s=$(date +%s); r=$(python3 $c 2>&1 | grep -v '^PARI stack size' | tail -1); echo "check $c: $r $(( $(date +%s) - s )) s" >> times.txt
done
