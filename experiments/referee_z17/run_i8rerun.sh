#!/bin/bash
# z17: t10's analysis from scratch on the 3 "neither" ratios and 15 of the 360 (seed 17), one fresh process each
cd "$(dirname "$0")"
ulimit -v 5000000
for r in $(cat z17_i8_sample.txt); do
  grep -q "\"m\": ${r%,*}, \"n\": ${r#*,}" z17_i8rerun.jsonl 2>/dev/null && continue
  s=$(date +%s)
  python3 -I z17_i8rerun.py $r 2>>z17_i8rerun.err | grep -v 'PARI stack' >> z17_i8rerun.jsonl
  echo "$r exit=${PIPESTATUS[0]} $(( $(date +%s) - s )) s" >> z17_i8rerun_times.txt
done
