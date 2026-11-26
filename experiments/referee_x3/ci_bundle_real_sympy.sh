#!/bin/bash
R=/home/claude/conductor/x3; L=/tmp/x3/logs
run() { local name=$1 dir=$2; shift 2; cd $dir; local t0=$(date +%s.%N)
  bash -c "$*" > $L/$name.out 2> $L/$name.err; local rc=$?; local t1=$(date +%s.%N)
  printf "%-6s rc=%d %6.1fs  %s\n" $name $rc $(echo "$t1-$t0"|bc) "$*" >> $L/summary_rs.txt; }
G=$R/experiments/census-wide/g4k; K=$R/experiments/census-wide/rank1
FF="grep -v -e '^PARI stack size' -e '^# matches'"
W="python3 -I $R/experiments/referee_x3/run_real_sympy.py"
run rs01 $G "$W g4k_exact.py | grep -q '^# matches'"
run rs03 $G "$W g4k_rank_exact.py | grep -q '^# matches'"
run rs07 $K "python3 -I ../../../scripts/g1_full/rank_indep.py | diff - ../../../scripts/g1_full/rank_indep.txt"
run rs08 $K "$W f2_general.py 3 | $FF | diff - f2_general_3.log"
run rs09 $K "$W f2_general_neg.py | $FF | diff - f2_general_negative.txt"
run rs10 $K "$W f2_d1_local2.py | $FF | diff - f2_d1_local2.txt"
run rs04 $G "$W g4k_rank_exact.py --negative | $FF | diff - g4k_rank_exact_negative.txt"
run rs11 $K "$W f2_d1_local2.py --negative | $FF | diff - f2_d1_local2_negative.txt"
run rs02 $G "$W g4k_exact.py --negative | $FF | diff - g4k_exact_negative.txt"
echo DONE >> $L/summary_rs.txt
