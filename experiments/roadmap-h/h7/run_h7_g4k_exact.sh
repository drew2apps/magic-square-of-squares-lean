#!/bin/bash
# H7: F4's exact g4.K chain (g4k_exact.py, g4k_rank_exact.py, unchanged, print-only mode) on (12201,7360),
# the one open G3 ratio where PARI's rank of g4.K is [1,1].  Outputs h7_g4k_exact.txt, h7_g4k_rank_exact.txt.
cd "$(dirname "$0")"; export MSS_REPO=${MSS_REPO:-../../..}; G=$MSS_REPO/experiments/census-wide/g4k
s=$(date +%s); timeout 3000 python3 $G/g4k_exact.py 12201,7360 2> err_g4k_exact.txt | grep -v '^PARI stack' > h7_g4k_exact.txt
echo "g4k_exact 12201,7360 exit=${PIPESTATUS[0]} $(( $(date +%s) - s )) s" >> times.txt
s=$(date +%s); timeout 3000 python3 $G/g4k_rank_exact.py 12201,7360 2> err_g4k_rank_exact.txt | grep -v '^PARI stack' > h7_g4k_rank_exact.txt
echo "g4k_rank_exact 12201,7360 exit=${PIPESTATUS[0]} $(( $(date +%s) - s )) s" >> times.txt
