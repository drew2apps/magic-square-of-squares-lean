#!/bin/bash
# z17: run one proposed CI line from the worktree root; record exit code and wall time (machine shared, load ~11 on 2 cores).
# usage: runci.sh NAME cmd...
name=$1; shift
out=experiments/referee_z17/ci/$name
s=$(date +%s.%N)
"$@" > $out.out 2> $out.err; r=$?
e=$(date +%s.%N)
printf '%s exit=%d wall=%.1f s  (%s)\n' "$name" $r $(echo "$e - $s" | bc) "$*" | tee -a experiments/referee_z17/ci/times.txt
