#!/bin/bash
# I2 on the ratios that need a good prime (one worker per argument N; default 1)
cd "$(dirname "$0")"
export G1_DIR=${G1_DIR:-$MSS_REPO/experiments/roadmap-g/g1} G2_DIR=${G2_DIR:-$MSS_REPO/experiments/roadmap-g/g2} H4_DIR=${H4_DIR:-$([ -d ../../roadmap-h/h4 ] && echo ../../roadmap-h/h4 || echo ../h4)}
N=${N:-1}
for i in $(seq 0 $((N-1))); do
  ( ulimit -v ${MEMKB:-5000000}; t=$SECONDS; python3 i2_shape.py ${LIST:-i2_needgood_list.txt} $i $N --write 2> >(grep -v 'stack size' > err_i2_${LIST:-n}_$i.txt); echo "i2_shape needgood worker $i of $N exit=$? $((SECONDS - t)) s" >> times.txt ) &
done
wait
