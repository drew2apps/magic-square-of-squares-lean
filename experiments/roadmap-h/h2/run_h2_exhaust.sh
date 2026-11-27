#!/bin/bash
# usage: bash run_h2_exhaust.sh LO HI i N   (worker i of N; times.txt)
cd "$(dirname "$0")"; export MSS_REPO=${MSS_REPO:-../../..}
s=$(date +%s); python3 h2_exhaust.py $1 $2 $3 $4 --write 2> err_ex_$1_$3.txt; echo "h2_exhaust $1 $2 worker $3 of $4 exit=$? $(( $(date +%s) - s )) s" >> times.txt
