#!/bin/bash
# I1 designed hunt: three seeds in parallel. Paths: the repository at MSS_REPO; roadmap H's folders at ../../roadmap-h (the repository's layout), or next to this one.
cd "$(dirname "$0")"
export G1_DIR=${G1_DIR:-$MSS_REPO/experiments/roadmap-g/g1} G2_DIR=${G2_DIR:-$MSS_REPO/experiments/roadmap-g/g2} H4_DIR=${H4_DIR:-$([ -d ../../roadmap-h/h4 ] && echo ../../roadmap-h/h4 || echo ../h4)}
for s in 1 2 3; do
  ( t=$SECONDS; python3 i1_hunt.py $s ${TAKE:-300} --write 2> >(grep -v 'stack size' > err_hunt_${I1_MODE:-A}_s$s.txt); echo "i1_hunt mode ${I1_MODE:-forced} seed $s take ${TAKE:-300} exit=$? $((SECONDS - t)) s" >> times.txt ) &
done
wait
