#!/bin/bash
# G1 full run: one worker per argument list of ranges; 20000 random ratios and 200 controls per range
cd "$(dirname "$0")"
for k in "$@"; do
  s=$(date +%s)
  nice python3 g1_stress.py $k 20000 200 --write > log_k$k.txt 2> err_k$k.txt
  echo "k=$k exit=$? $(( $(date +%s) - s )) s" >> times.txt
done
