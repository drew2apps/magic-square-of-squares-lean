#!/bin/bash
# usage: bash run_h1_scan.sh i N   (worker i of N over k = 13..29); MSS_REPO points at the repository
cd "$(dirname "$0")"
i=$1; N=$2; j=0
for k in $(seq 29 -1 13); do
  if [ $((j % N)) -eq $i ]; then
    s=$(date +%s); python3 h1_scan.py $k 1177 --write 2> err_scan_$k.txt; e=$?
    echo "scan k=$k exit=$e $(( $(date +%s) - s )) s" >> times.txt
  fi
  j=$((j+1))
done
