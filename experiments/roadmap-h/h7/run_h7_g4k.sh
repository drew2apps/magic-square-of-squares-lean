#!/bin/bash
# H7: F4's g4.K route (experiments/census-wide/g4k/dm_g4k.py, unchanged) on G3's open ratios, one ratio at a time.
# Output h7_g4k.jsonl (one JSON record per ratio, as dm_g4k.py prints it); times.txt.  MSS_REPO as usual.
cd "$(dirname "$0")"; export MSS_REPO=${MSS_REPO:-../../..}
: > h7_g4k.jsonl.part
for r in "$@"; do
  s=$(date +%s); timeout 2400 python3 $MSS_REPO/experiments/census-wide/g4k/dm_g4k.py $r 2> err_g4k_${r/,/_}.txt | grep -v '^PARI stack' >> h7_g4k.jsonl.part
  echo "dm_g4k $r exit=${PIPESTATUS[0]} $(( $(date +%s) - s )) s" >> times.txt
done
mv h7_g4k.jsonl.part h7_g4k.jsonl
