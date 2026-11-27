#!/bin/bash
# usage: bash run_h6.sh m n   (splitting, then ranks; times.txt)
cd "$(dirname "$0")"
s=$(date +%s); python3 h6_split.py $1 $2 --write; echo "h6_split $1 $2 exit=$? $(( $(date +%s) - s )) s" >> times.txt
s=$(date +%s); python3 h6_ranks.py $1 $2 --write 2> err_ranks_$1_$2.txt; echo "h6_ranks $1 $2 exit=$? $(( $(date +%s) - s )) s" >> times.txt
