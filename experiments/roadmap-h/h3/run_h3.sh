#!/bin/bash
# usage: bash run_h3.sh   (reruns every mode without --write, compares with the committed outputs; times.txt)
cd "$(dirname "$0")"; export MSS_REPO=${MSS_REPO:-../../..}
for a in "" "--negative" "--check 2000"; do
  s=$(date +%s); r=$(python3 h3_tables.py $a 2>/dev/null | tail -1); echo "h3_tables $a: $r $(( $(date +%s) - s )) s" >> times.txt
done
