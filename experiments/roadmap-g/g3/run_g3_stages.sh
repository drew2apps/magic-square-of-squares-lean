#!/bin/bash
# waits for each chunk's second descent, then runs pipe_g3.py on it.  MSS_REPO is the repository (default ../../..).
cd "$(dirname "$0")"
export MSS_REPO=${MSS_REPO:-../../..}
for lo in $(seq 11000 1000 29000); do
  until [ -f done_$lo ]; do sleep 60; done
  [ -f stage_$lo.json ] || { s=$(date +%s); python3 pipe_g3.py $lo >> stages.log 2>&1; echo "stages $lo: exit $? $(( $(date +%s)-s )) s" >> times.txt; }
done
