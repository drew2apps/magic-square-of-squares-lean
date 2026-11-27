#!/bin/bash
# G4 second pass on E+ and E-: worker i of N.  Writes pass2_<i>_<N>.jsonl (resumes it).
cd "$(dirname "$0")"
export MSS_REPO=${MSS_REPO:-../../..}
i=$1; N=$2; s=$(date +%s)
G4_I=$i G4_N=$N nice python3 g4_pass2.py g4_0_3.json g4_1_3.json g4_2_3.json --out pass2_${i}_${N}.jsonl --curves E+,E- 2> err_pass2_${i}_${N}.txt
echo "pass2 worker $i/$N exit=$? $(( $(date +%s) - s )) s" >> times.txt
