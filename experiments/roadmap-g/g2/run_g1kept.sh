#!/bin/bash
# G2 kill sets for G1's 63 kept ratios, one process per ratio under an 8 GB limit on virtual memory (run together in one process,
# they need more than 8 GB of real memory).  A ratio that runs out of memory is written as SKIPPED.  Writes g2_g1_kept_0.txt (in the order of g1_kept.txt; delete it to start over).
cd "$(dirname "$0")"
grep -h "class dim" ../g1/g1_k*.txt | grep -v degenerate | sed 's/^(\([0-9]*\),\([0-9]*\)).*/\1,\2/' > g1_kept.txt
touch g2_g1_kept_0.txt; d=$(mktemp -d)    # resumes: ratios already written are kept
while read r; do
  grep -q "^$r " g2_g1_kept_0.txt && continue
  echo "$r" > $d/one.txt
  line=$( (ulimit -v 8000000; nice python3 g2_killsets.py $d/one.txt 0 1 2>/dev/null) | grep bad-alone)
  [ -n "$line" ] || line="$r  SKIPPED (out of memory at 8 GB)"
  echo "$line" >> g2_g1_kept_0.txt
done < g1_kept.txt
rm -r $d
