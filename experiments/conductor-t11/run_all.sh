#!/bin/bash
# runs t11_theta.py on each curve with --write and records wall times (runlogs/times.txt); not part of the outputs
cd "$(dirname "$0")"
mkdir -p runlogs
for lab in "$@"; do
  t0=$(date +%s)
  timeout 1500 python3 -I t11_theta.py "$lab" --write > "runlogs/$(echo "$lab" | tr -c 'A-Za-z0-9\n' '_').log" 2>&1
  echo "$lab exit $? $(( $(date +%s) - t0 )) s" >> runlogs/times.txt
done
