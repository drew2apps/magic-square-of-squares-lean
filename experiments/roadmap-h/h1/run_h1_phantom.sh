#!/bin/bash
# usage: bash run_h1_phantom.sh   (times.txt)
cd "$(dirname "$0")"; export MSS_REPO=${MSS_REPO:-../../..}
s=$(date +%s); python3 h1_phantom.py --write 2> err_phantom.txt; echo "h1_phantom exit=$? $(( $(date +%s) - s )) s" >> times.txt
