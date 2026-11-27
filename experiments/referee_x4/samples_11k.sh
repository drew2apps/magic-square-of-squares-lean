#!/bin/bash
# Review x4 (run 6): the item-2 reruns, with the repository's programs, into OUT (about 25 min on 2 cores; not for CI).
# Written by an AI model; not read by a person.  usage: bash samples_11k.sh OUT   then copy OUT/* to samples_11k/ to compare.
# The repository programs name the session's folders (experiments/README.md); PYTHONPATH supplies them instead.
set -eu
OUT=$(mkdir -p "$1" && cd "$1" && pwd); HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../.." && pwd); S=$HERE/samples_11k
export PYTHONPATH=$ROOT/scripts/gauss_descent:$ROOT/scripts/second_descent:$ROOT/experiments/census-wide
cd "$ROOT/experiments/census-wide"
for k in 0 1; do python3 wide_pari.py "$S/sd_in.json" $k 2 > "$OUT/sd_$k.jsonl" & done; wait
python3 "$ROOT/scripts/g1_full/global.py" "$S/g_in.json" "$OUT/g_out.json" > "$OUT/g.log" 2>&1
python3 "$ROOT/experiments/census2000/dm_pari.py" "$S/dm_in.json" "$OUT/dm_out.jsonl" > "$OUT/dm.log" 2>&1
cd route
for d in 220962238 -220962238; do
  python3 route_pari.py 10512,145 $d 30000 > "$OUT/route_10512_$d.txt" 2>&1
  python3 rank_hard.py 10512,145,$d 8 2>/dev/null | grep -v 'PARI stack' > "$OUT/rh_10512_$d.txt"
done
cd "$ROOT/scripts/cover_descent/structure"
python3 class_selmer_linalg.py 11000 0 --from 10000 --list "$OUT/list_10000_11000.txt" > "$OUT/linalg_10000_11000.txt"
cd .. && python3 tuple_descent.py batch magic "$OUT/list_10000_11000.txt" 1 > "$OUT/tuple_10000_11000.txt"
