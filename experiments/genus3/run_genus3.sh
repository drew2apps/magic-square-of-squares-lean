#!/bin/bash
# usage: bash run_genus3.sh   (every script of this folder compared with its committed output; times.txt)
# Written by an AI model (build agent t12); not read by a person.
cd "$(dirname "$0")"
: > times.txt
r() { local name="$1" out="$2"; shift 2; local s=$(date +%s); python3 -I "$@" 2>/dev/null | grep -v 'PARI stack' | diff -q - "$out" > /dev/null; echo "$name exit=$? $(( $(date +%s) - s )) s" >> times.txt; }
r genus3_split genus3_split.txt genus3_split.py
r genus3_points genus3_points.txt genus3_points.py
r "genus3_points --negative" genus3_points_negative.txt genus3_points.py --negative
r genus3_isog genus3_isog_627_364.txt genus3_isog.py
r genus3_mult genus3_mult.txt genus3_mult.py
r kcurve_exact kcurve_exact.txt kcurve_exact.py
r "kcurve_exact --negative" kcurve_exact_negative.txt kcurve_exact.py --negative
r genus3_rank_upper genus3_rank_upper.txt genus3_rank_upper.py
cd session
for s in verify_swap_and_genus3 genus3_ranks kcurve_not_qcurve; do r "session/$s" $s.txt $s.py; done
r session/isogeny_classes isogeny_classes_627_364.txt isogeny_classes.py
mv times.txt ../times_session.txt 2>/dev/null; cd ..; cat times_session.txt >> times.txt; rm -f times_session.txt
