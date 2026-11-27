#!/bin/bash
# usage: bash run_h1_allphantoms.sh   (times.txt)
cd "$(dirname "$0")"; export MSS_REPO=${MSS_REPO:-../../..}
s=$(date +%s); python3 h1_allphantoms.py --write 2> err_allphantoms.txt; echo "h1_allphantoms exit=$? $(( $(date +%s) - s )) s" >> times.txt
