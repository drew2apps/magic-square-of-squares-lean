#!/bin/bash
# F0: the class-criterion step of the census, rerun with the repository's own programs at commit 24d6ef6 (unchanged).
# The session ran it through a PARI stand-in for sympy (not committed); this copy runs the programs with sympy itself
# (packaging edit, run 6).  Written by an AI model (Claude); not read by a person.  Experiment, not counted.  Each step is compared with the committed file by diff; the log records the result
# and the time.  usage: bash class_step.sh REPO OUTDIR
set -u -o pipefail
REPO=$1; OUT=$2; mkdir -p "$OUT"
cd "$REPO" && echo "commit $(git rev-parse HEAD)" > "$OUT/commit.txt"
cd "$REPO/scripts/cover_descent"
t0=$(date +%s)
# (1) the four-class descent on every ratio 1000 < m <= 5000 that avoids the 484 classes (CI job scan-5000, first line)
python3 tuple_descent.py batch magic ../gauss_descent/scan_certificates.txt 1 > "$OUT/tuple_scan_5000.txt"
diff "$OUT/tuple_scan_5000.txt" tuple_scan_5000.txt > "$OUT/tuple_scan_5000.diff"; echo "(1) tuple_scan_5000: diff exit $?, $(( $(date +%s) - t0 )) s" >> "$OUT/log.txt"
t0=$(date +%s); cd structure
# (2) the criterion by linear algebra, 5000 < m <= 10000, and its list (CI job scan-15000, first half)
python3 class_selmer_linalg.py 10000 0 --from 5000 --list "$OUT/list_5000_10000.txt" > "$OUT/linalg_5000_10000.txt"
diff "$OUT/linalg_5000_10000.txt" <(head -n 4 class_selmer_linalg_15000.txt) > "$OUT/linalg_5000_10000.diff"; e1=$?
diff "$OUT/list_5000_10000.txt" <(awk -F'[ ,]' '$1 <= 10000' class_selmer_15000_list.txt) > "$OUT/list_5000_10000.diff"; e2=$?
echo "(2) criterion 5000 < m <= 10000: summary diff exit $e1, list diff exit $e2, $(( $(date +%s) - t0 )) s" >> "$OUT/log.txt"
t0=$(date +%s)
# (3) the criterion, 10000 < m <= 11000 (a new range for the summary; the list is compared with the committed list)
python3 class_selmer_linalg.py 11000 0 --from 10000 --list "$OUT/list_10000_11000.txt" > "$OUT/linalg_10000_11000.txt"
diff "$OUT/list_10000_11000.txt" <(awk -F'[ ,]' '$1 > 10000 && $1 <= 11000' class_selmer_15000_list.txt) > "$OUT/list_10000_11000.diff"
echo "(3) criterion 10000 < m <= 11000: list diff exit $?, $(( $(date +%s) - t0 )) s" >> "$OUT/log.txt"
t0=$(date +%s); cd ..
# (4) the four-class descent on the listed ratios with m <= 11000 (CI job scan-15000, second part, restricted)
cat "$OUT/list_5000_10000.txt" "$OUT/list_10000_11000.txt" > "$OUT/list_5000_11000.txt"
python3 tuple_descent.py batch magic "$OUT/list_5000_11000.txt" 1 > "$OUT/tuple_selmer_11000.txt"
diff "$OUT/tuple_selmer_11000.txt" <(awk -F'[(),]' '$2 <= 11000' tuple_selmer_15000.txt) > "$OUT/tuple_selmer_11000.diff"
echo "(4) tuple descent on the list, 5000 < m <= 11000: diff exit $?, $(grep -c 'nothing is left$' "$OUT/tuple_selmer_11000.txt") of $(wc -l < "$OUT/list_5000_11000.txt") with nothing left, $(( $(date +%s) - t0 )) s" >> "$OUT/log.txt"
echo done >> "$OUT/log.txt"
