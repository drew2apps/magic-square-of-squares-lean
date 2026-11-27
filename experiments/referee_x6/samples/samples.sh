#!/bin/bash
# Review x6 (run 6), written by an AI model; not read by a person.  The reruns of roadmap-g/g3 and g4 stages with the
# repository's programs under PARI 2.17.4, as they were made (about 5 minutes).  census-wide's sd_pari.py, wide_pari.py and
# route_pari.py name the session's folders (experiments/README.md, "Rerunning"), so patched copies are made in a scratch
# folder first; nothing in the repository is changed.  Writes into the scratch folder only; compare.py compares the
# copies committed here with the committed records.   usage: bash samples.sh [SCRATCH]
set -e
H=$(cd "$(dirname "$0")" && pwd); R=$(cd "$H/../../.." && pwd); S=${1:-/tmp/x6_samples}; mkdir -p $S/cw $S/route
O='/tmp/claude-0/-home-user/9ec7a8d3-42dc-53c8-8518-1a663ac809f1/scratchpad/repo/scripts'
sed "s#$O/#$R/scripts/#" $R/experiments/census-wide/sd_pari.py > $S/cw/sd_pari.py
sed "s#/mnt/project-files/mss-experiments/census-wide#$S/cw#" $R/experiments/census-wide/wide_pari.py > $S/cw/wide_pari.py
sed "s#$O/second_descent#$R/scripts/second_descent#; s#sys.path.insert(0, '..')#sys.path.insert(0, '$S/cw')#" $R/experiments/census-wide/route/route_pari.py > $S/route/route_pari.py
cd $S/cw
python3 -I wide_pari.py $H/sample.json 0 1 2>/dev/null | grep -v 'PARI stack' > sample_sd.jsonl
python3 -I $R/scripts/g1_full/global.py $H/r0_in.json r0_out.json > r0.log 2>&1
python3 -I $R/experiments/census2000/dm_pari.py $H/dm_in.json dm_out.jsonl > dm.log 2>&1
cd $S/route && python3 -I route_pari.py 12416,5313 251858174 > log_12416.txt 2>&1
cd $S && python3 -I $R/experiments/roadmap-g/g4/g4_search.py $H/g4_in.json 0 1 2>/dev/null | grep -v 'PARI stack' > g4_rerun.jsonl
echo "outputs in $S"
