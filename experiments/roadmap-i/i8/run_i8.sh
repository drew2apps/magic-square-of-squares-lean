#!/bin/bash
# I8: N workers over i2_needgood_list.txt, each a fresh process every 40 ratios (PARI's memory grows); then the merge.
# MSS_REPO: a checkout of main with experiments/conductor-t10 (default ../../..).
cd "$(dirname "$0")"; N=${N:-3}
ulimit -v 7000000
for i in $(seq 0 $((N-1))); do
  ( s=$(date +%s); while true; do python3 i8_pairct.py $i $N --max 40 --write 2>>err_i8_$i.txt | grep -v 'PARI stack'; r=${PIPESTATUS[0]}; [ $r -eq 3 ] || break; done
    echo "i8 worker $i of $N exit=$r $(( $(date +%s) - s )) s" >> times.txt ) &
done
wait
python3 i8_pairct.py --merge $N --write 2>/dev/null | grep -v 'PARI stack'
