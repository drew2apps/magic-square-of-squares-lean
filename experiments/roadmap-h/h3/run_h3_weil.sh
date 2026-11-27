#!/bin/bash
# usage: bash run_h3_weil.sh   (writes h3_weil.txt, h3_weil.json, h3_weil_negative.txt; times.txt)
cd "$(dirname "$0")"; export MSS_REPO=${MSS_REPO:-../../..}
s=$(date +%s); python3 h3_weil.py --write 2> err_weil.txt; echo "h3_weil --write exit=$? $(( $(date +%s) - s )) s" >> times.txt
s=$(date +%s); python3 h3_weil.py --negative --write 2>> err_weil.txt; echo "h3_weil --negative --write exit=$? $(( $(date +%s) - s )) s" >> times.txt
