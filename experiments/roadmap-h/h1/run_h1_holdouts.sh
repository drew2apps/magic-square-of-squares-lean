#!/bin/bash
# usage: bash run_h1_holdouts.sh   (times.txt)
cd "$(dirname "$0")"; export MSS_REPO=${MSS_REPO:-../../..}
s=$(date +%s); python3 h1_holdouts.py --write 2> err_holdouts.txt; echo "h1_holdouts exit=$? $(( $(date +%s) - s )) s" >> times.txt
