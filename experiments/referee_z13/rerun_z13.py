"""Reviewer z13: rerun z13_indep.run from scratch on a fixed subset of the jobs (every STEP-th line of jobs.txt, plus
the first curve of each category) and compare every field but 'time' with the committed records in
indep_a.jsonl, indep_b.jsonl.  Written by an AI model (reviewer z13), not read by a person.
usage: rerun_z13.py [STEP]   (default 25)"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
STEP = int(sys.argv[1]) if len(sys.argv) > 1 else 25
sys.argv = sys.argv[:1]
sys.path.insert(0, HERE)
import z13_indep as Z
jobs = [l.split() for l in open(os.path.join(HERE, 'jobs.txt')) if l.strip()]
pick = []; seen = set()
for i, j in enumerate(jobs):
    if i % STEP == 0 or j[3] not in seen: pick.append(j)
    seen.add(j[3])
old = {}
for fn in ('indep_a.jsonl', 'indep_b.jsonl'):
    for l in open(os.path.join(HERE, fn)):
        d = json.loads(l); old[(d['m'], d['n'], d['cv'])] = d
same = diff = 0
for m, n, cv, cat in pick:
    r = json.loads(json.dumps(Z.run(int(m), int(n), cv, cat)))
    o = old[(int(m), int(n), cv)]
    keys = sorted((set(r) | set(o)) - {'time'})
    bad = [k for k in keys if r.get(k) != o.get(k)]
    if bad: diff += 1; print('DIFFERS', m, n, cv, bad)
    else: same += 1
    print(m, n, cv, cat, 'pairs', r['cmp']['agree'] + r['cmp']['disagree'], 'disagree', r['cmp']['disagree'], 'same as committed' if not bad else 'DIFFERENT')
print('rerun', len(pick), 'curves: identical', same, 'different', diff)
