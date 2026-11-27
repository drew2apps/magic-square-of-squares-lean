#!/bin/bash
# Review z9 of roadmap H (run 7), written by an AI model; not read by a person.  Reruns a sample of G3's block
# 13000 < m <= 14000 with the repository's programs, as referee_x6/samples/samples.sh did for run 6: census-wide's
# sd_pari.py and wide_pari.py name the session's folders, so patched copies are made in a scratch folder first; nothing
# in the repository is changed.  Inputs in g3_samples_in/: second descent on (13001,9679), (13348,3795), (13447,1303);
# rank-0 quotients (scripts/g1_full/global.py) on (13015,4872), (13049,5543), (13348,3795); Demjanenko-Manin
# (census2000/dm_pari.py) on (13439,7440), (13447,1303) with their partners.  Then g3_samples_compare.py.
# usage: bash g3_samples.sh [SCRATCH]   (about 3 minutes of CPU; not for CI)
set -e
H=$(cd "$(dirname "$0")" && pwd); R=$(cd "$H/../.." && pwd); S=${1:-/tmp/z9_g3_samples}; mkdir -p $S/cw
O='/tmp/claude-0/-home-user/9ec7a8d3-42dc-53c8-8518-1a663ac809f1/scratchpad/repo/scripts'
sed "s#$O/#$R/scripts/#" $R/experiments/census-wide/sd_pari.py > $S/cw/sd_pari.py
sed "s#/mnt/project-files/mss-experiments/census-wide#$S/cw#" $R/experiments/census-wide/wide_pari.py > $S/cw/wide_pari.py
cd $S/cw
nice python3 -I wide_pari.py $H/g3_samples_in/sample.json 0 1 2>/dev/null | grep -v 'PARI stack' > $S/sample_sd.jsonl
nice python3 -I $R/scripts/g1_full/global.py $H/g3_samples_in/r0_in.json $S/r0_out.json > $S/r0.log 2>&1
nice python3 -I $R/experiments/census2000/dm_pari.py $H/g3_samples_in/dm_in.json $S/dm_out.jsonl > $S/dm.log 2>&1
cd $H && python3 -I g3_samples_compare.py $S
