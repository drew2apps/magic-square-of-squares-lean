#!/bin/bash
# I4: G3's pipeline (unchanged: census-wide/wide_pari.py, then pipe_g3.py) on 14000 < m <= 20000, chunks of 1000 in m,
# two processes per chunk, CHUNKS chunks at a time (default 1).  Usage: run_i4.sh LO1 LO2 ...  (e.g. 14000 15000)
cd "$(dirname "$0")"
export MSS_REPO=${MSS_REPO:-../../..}
for lo in "$@"; do
  hi=$((lo+1000)); [ -f done_$lo ] && continue
  python3 -c "import json; P=json.load(open('fam_11000_30000.json'))['open']; json.dump([p for p in P if $lo < p[0] <= $hi], open('chunk_$lo.json','w'))"
  s=$(date +%s)
  for k in 0 1; do nice -n 5 python3 ../../census-wide/wide_pari.py chunk_$lo.json $k 2 > sd_${lo}_$k.jsonl 2> sd_${lo}_$k.err & done; wait
  echo "chunk $lo-$hi: $(( $(date +%s)-s )) s" >> times.txt; touch done_$lo
  s=$(date +%s); python3 pipe_g3.py $lo >> stages.log 2>&1; echo "stages $lo: exit $? $(( $(date +%s)-s )) s" >> times.txt
done
