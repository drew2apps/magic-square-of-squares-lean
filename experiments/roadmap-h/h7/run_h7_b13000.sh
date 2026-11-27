#!/bin/bash
# H7: G3 block 13000 < m <= 14000 left two open ratios, (13348,3795) and (13899,1180).  On both: F4's g4.K route
# (dm_g4k.py), F4's exact chain (g4k_exact.py) and H7's exact rank step (h7_g4k_rank1.py), all unchanged, print-only.
# Outputs h7_b13000_g4k.jsonl, h7_b13000_g4k_exact.txt, h7_b13000_rank1.txt; times.txt.  MSS_REPO, H1_DIR as usual.
cd "$(dirname "$0")"; export MSS_REPO=${MSS_REPO:-../../..}; G=$MSS_REPO/experiments/census-wide/g4k
: > h7_b13000_g4k.jsonl; : > h7_b13000_g4k_exact.txt
for r in 13348,3795 13899,1180; do
  s=$(date +%s); timeout 2400 python3 $G/dm_g4k.py $r 2> err_b13000_${r/,/_}.txt | grep -v '^PARI stack' >> h7_b13000_g4k.jsonl
  echo "dm_g4k $r exit=${PIPESTATUS[0]} $(( $(date +%s) - s )) s" >> times.txt
  s=$(date +%s); timeout 3000 python3 $G/g4k_exact.py $r 2>> err_b13000_${r/,/_}.txt | grep -v '^PARI stack' >> h7_b13000_g4k_exact.txt
  echo "g4k_exact $r exit=${PIPESTATUS[0]} $(( $(date +%s) - s )) s" >> times.txt
done
s=$(date +%s); python3 h7_g4k_rank1.py 13348,3795 13899,1180 2>> err_b13000_rank1.txt | grep -v '^PARI stack' > h7_b13000_rank1.txt
echo "h7_g4k_rank1 13348,3795 13899,1180 exit=${PIPESTATUS[0]} $(( $(date +%s) - s )) s" >> times.txt
