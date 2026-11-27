#!/bin/bash
cd "$(dirname "$0")"
export MSS_REPO=${MSS_REPO:-../../..}
i=$1; N=$2; s=$(date +%s)
nice python3 -u g4_search.py square_20000.json $i $N 8 --write > log_${i}_${N}.jsonl 2> err_${i}_${N}.txt
echo "worker $i/$N exit=$? $(( $(date +%s) - s )) s" >> times.txt
