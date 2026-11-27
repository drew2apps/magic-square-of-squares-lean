#!/bin/bash
# usage: bash run_h1_partners.sh   (times.txt)
cd "$(dirname "$0")"; export MSS_REPO=${MSS_REPO:-../../..}
s=$(date +%s); python3 h1_partners.py --write 2> err_partners.txt; echo "h1_partners exit=$? $(( $(date +%s) - s )) s" >> times.txt
