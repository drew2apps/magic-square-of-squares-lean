#!/bin/bash
# F5 (6557,324) w*f0*f8: the streamed L-value; each stage logs its exit and wall time
cd "$(dirname "$0")/.."
export MSS_PAIR=6557,324 MSS_CURVE='w*f0*f8'
D=r6557
st() { s=$(date +%s); "$@"; echo "$* : exit $? $(( $(date +%s)-s )) s" >> $D/times.txt; }
st python3 lvalue_stream.py --primes $D 2> $D/primes.err
for J in 0 1 2 3; do st nice -n 5 python3 lvalue_stream.py --ap 4 $J $D 2> $D/ap$J.err & done; wait
for g in A B C; do st nice -n 5 python3 lvalue_stream.py --series $g $D 2> $D/series_$g.err & done; wait
st python3 lvalue_stream.py --merge $D --write > $D/merge.out 2> $D/merge.err
st python3 lvalue_stream.py --negative $D --write > $D/negative.out 2> $D/negative.err
st python3 lvalue_stream.py --apcheck $D --write > $D/apcheck.out 2> $D/apcheck.err
