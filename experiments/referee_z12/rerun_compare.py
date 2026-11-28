"""Run 10 (the conductor's, for CI): rerun review z12's full test (fulltest.run: rank_second's complete 2-descent and
Cassels-Tate pairing against ct08's closed formula) on a few curves from scratch and compare every row with the
committed records in full0.jsonl / full1.jsonl.  Prints one line per curve and a summary."""
import os, sys, json, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
CURVES = [(5623, 1017, '08'), (3320, 2303, '08'), (771, 251, '08'), (1323, 803, '08'), (511, 260, 'LU'), (988, 525, 'LU')] + [
    tuple(int(x) for x in a.split(',')[:2]) + (a.split(',')[2],) for a in sys.argv[1:]]
stored = {}
for f in ('full0.jsonl', 'full1.jsonl'):
    for l in open(os.path.join(HERE, f)):
        r = json.loads(l)
        if 'error' not in r: stored[(r['m'], r['n'], r['cv'])] = r
import fulltest
key = lambda row: (tuple(row['a']), row['inSel'], tuple(row['rs']), json.dumps(row['ct'], sort_keys=True),
                   json.dumps(row.get('ct_all', {}), sort_keys=True))
bad = 0
with tempfile.TemporaryDirectory() as d:
    out = os.path.join(d, 'rerun.jsonl')
    for (m, n, cv) in CURVES:
        r = json.loads(json.dumps(fulltest.run(m, n, cv, out, seeds=(1, 7))))   # the committed form (string keys, lists)
        s = stored[(m, n, cv)]
        same = r['sel'] == s['sel'] and sorted(map(key, r['rows'])) == sorted(map(key, s['rows']))
        nrs = sum(len(x['rs']) for x in r['rows'])
        agree = all(x['rs'] == x['ct']['1'] for x in r['rows'] if x['inSel'] and x['ct'])
        print((m, n, cv), 'dim Sel2', r['dim'], 'rows', len(r['rows']), 'pairings', nrs,
              'same as committed' if same else 'DIFFERS from committed', '; ct08 = rank_second' if agree else '; DISAGREE')
        bad += (not same) + (not agree)
print('all same, all agree' if bad == 0 else 'FAILURES: %d' % bad)
sys.exit(1 if bad else 0)
