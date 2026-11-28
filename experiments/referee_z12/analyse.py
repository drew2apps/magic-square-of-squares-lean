import os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__))  # repository copy: paths relative to this folder
"""Count agreement of ct08.pairing with rank_second (DFull.pairing) over distinct (ratio, curve, a, b) triples."""
import json, sys, glob, collections
from fractions import Fraction as Fr
recs = []
for f in sorted(glob.glob(_os.path.join(HERE, 'full[01].jsonl'))):
    for l in open(f): recs.append(json.loads(l))
cnt = collections.Counter(); dis = []; seeddep = []; lin = collections.Counter(); notsel = []; errs = [r for r in recs if 'error' in r]
ratios = collections.defaultdict(set); none = collections.Counter(); kinds = collections.Counter()
for r in recs:
    if 'error' in r: continue
    ratios[r['cv']].add((r['m'], r['n']))
    sel = [tuple(b) for b in r['sel']]
    for row in r['rows']:
        a = tuple(row['a']); kind = 'tuple' if any(w[0] == 'tuple' for w in row['why']) else 'sel'
        if not row['inSel']: notsel.append((r['m'], r['n'], r['cv'], a, kind)); continue
        rs = row['rs']
        if not row['ct']: cnt[(r['cv'], 'no ct (a not (1,d,d))')] += len(sel); continue
        c1 = row['ct']['1']; c7 = row['ct']['7']
        for b, h, v1, v7 in zip(sel, rs, c1, c7):
            if h is None: none[(r['cv'], 'rs None')] += 1; continue
            if v1 is None or v7 is None: none[(r['cv'], 'ct None')] += 1
            if v1 != v7: seeddep.append((r['m'], r['n'], r['cv'], a, b, v1, v7))
            ok = v1 == h
            cnt[(r['cv'], kind, 'agree' if ok else 'DISAGREE')] += 1
            kinds[(r['cv'], kind, h)] += 1
            if not ok: dis.append((r['m'], r['n'], r['cv'], a, b, h, v1))
        # linearity: ct_all element i = product of sel[k] over bits of i; expected value = product of rs over those bits
        if row['ct_all']:
            for i, v in enumerate(row['ct_all']['v']):
                e = 1
                for k in range(len(sel)):
                    if i >> k & 1: e *= rs[k] if rs[k] is not None else 0
                lin[(r['cv'], 'ct over all Sel2 = linear extension of rs' if v == e else 'LINEARITY/AGREEMENT FAIL')] += 1
print('curves done:', {k: len(v) for k, v in ratios.items()}, 'errors:', len(errs), errs[:3])
print('distinct (ratio,curve,a,b) triples, b over the Sel2 basis:')
for k, v in sorted(cnt.items()): print('  ', k, v)
print('values of rank_second on those triples:', dict(sorted(kinds.items())))
print('None counts:', dict(none))
print('ct08 seed 1 vs seed 7 differences:', len(seeddep), seeddep[:3])
print('ct08 over all of Sel2 vs linear extension of rank_second:', dict(lin))
print('tuple classes / a not in Sel2:', len(notsel), notsel[:5])
print('disagreements:', len(dis))
for d in sorted(dis, key=lambda d: (d[0] * d[1], d[0]))[:10]: print('  ', d)
