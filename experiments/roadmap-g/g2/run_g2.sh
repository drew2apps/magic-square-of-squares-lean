#!/bin/bash
# G2 kill sets, worker i of N, on the repository's two lists (the bundle ran on byte-identical copies of them, left out here).
# Writes g2_tuple_selmer_997_<i>.txt and g2_class_selmer_15000_list_<i>.txt; g2_997_all.txt and g2_15000_all.txt are the
# workers' outputs concatenated (the committed per-worker files are left out; sorted, they equal the _all files).
cd "$(dirname "$0")"
i=$1; N=$2; R=${MSS_REPO:-../../..}/scripts/cover_descent
for L in $R/tuple_selmer_997.txt $R/structure/class_selmer_15000_list.txt; do
  b=$(basename $L .txt)
  s=$(date +%s); nice python3 g2_killsets.py $L $i $N --write > log_${b}_$i.txt 2> err_${b}_$i.txt
  echo "$(basename $L) worker $i/$N exit=$? $(( $(date +%s) - s )) s" >> times.txt
done
