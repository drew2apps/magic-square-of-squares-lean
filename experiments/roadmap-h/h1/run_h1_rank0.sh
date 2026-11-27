#!/bin/bash
# usage: bash run_h1_rank0.sh   (times.txt)
cd "$(dirname "$0")"; export MSS_REPO=${MSS_REPO:-../../..}
s=$(date +%s); python3 h1_rank0.py --write 2> err_rank0.txt; echo "h1_rank0 exit=$? $(( $(date +%s) - s )) s" >> times.txt
