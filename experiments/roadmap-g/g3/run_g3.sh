#!/bin/bash
# G3: second descent (wide_pari.py, unchanged) on the degenerate-type ratios, by chunks of 1000 in m, two processes each.
cd "$(dirname "$0")"
for lo in $(seq 11000 1000 29000); do
  hi=$((lo+1000)); [ -s sd_${lo}_1.jsonl ] && grep -q '"pair"' sd_${lo}_1.jsonl && [ -f done_$lo ] && continue
  python3 -c "import json; P=json.load(open('fam_11000_30000.json'))['open']; json.dump([p for p in P if $lo < p[0] <= $hi], open('chunk_$lo.json','w'))"
  s=$(date +%s)
  for k in 0 1; do nice -n 5 python3 ../../census-wide/wide_pari.py chunk_$lo.json $k 2 > sd_${lo}_$k.jsonl 2> sd_${lo}_$k.err & done; wait
  echo "chunk $lo-$hi: $(( $(date +%s)-s )) s" >> times.txt; touch done_$lo
done
