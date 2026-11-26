#!/bin/bash
# run the bundle's eleven CI lines; record pass/fail and wall time
R=/home/claude/conductor/x3
L=/tmp/x3/logs
F() { grep -v -e '^PARI stack size' -e '^# matches'; }
run() { # name dir cmd
  local name=$1 dir=$2; shift 2
  cd $dir; local t0=$(date +%s.%N)
  bash -c "$*" > $L/$name.out 2> $L/$name.err; local rc=$?
  local t1=$(date +%s.%N)
  printf "%-4s rc=%d %6.1fs  %s\n" $name $rc $(echo "$t1-$t0"|bc) "$*" >> $L/summary.txt
}
export -f F
G=$R/experiments/census-wide/g4k; K=$R/experiments/census-wide/rank1
FF="grep -v -e '^PARI stack size' -e '^# matches'"
run c01 $G "python3 -I g4k_exact.py | grep -q '^# matches'"
run c03 $G "python3 -I g4k_rank_exact.py | grep -q '^# matches'"
run c05 $G "python3 -I dm_g4k.py | grep -q '^# matches'"
run c07 $K "python3 -I run_shim.py ../../../scripts/g1_full/rank_indep.py | diff - ../../../scripts/g1_full/rank_indep.txt"
run c08 $K "python3 -I f2_general.py 3 | $FF | diff - f2_general_3.log"
run c09 $K "python3 -I f2_general_neg.py | $FF | diff - f2_general_negative.txt"
run c10 $K "python3 -I f2_d1_local2.py | $FF | diff - f2_d1_local2.txt"
run c04 $G "python3 -I g4k_rank_exact.py --negative | $FF | diff - g4k_rank_exact_negative.txt"
run c06 $G "python3 -I dm_g4k.py --negative | $FF | diff - dm_g4k_negative.txt"
run c11 $K "python3 -I f2_d1_local2.py --negative | $FF | diff - f2_d1_local2_negative.txt"
run c02 $G "python3 -I g4k_exact.py --negative | $FF | diff - g4k_exact_negative.txt"
echo DONE >> $L/summary.txt
