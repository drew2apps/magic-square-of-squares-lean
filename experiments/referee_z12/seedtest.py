import os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__))  # repository copy: paths relative to this folder
import json, sys, collections, os
sys.path.insert(0, _os.path.join(HERE, '..', 'research-pairing')); os.chdir(_os.path.join(HERE, '..', 'conductor-t10'))
import ct08
from parse import sqf
c = collections.Counter(); ex = []
for f in sys.argv[1:]:
    for l in open(f):
        r = json.loads(l)
        m, n = r['m'], r['n']
        A, C = (m * m - n * n, m * m + n * n) if r['cv'] == '08' else (8 * m * n, 4 * (m * m + n * n))
        for row in r['rows'][:6]:
            if not row['ct']: continue
            dl = row['ap'][1]
            for b, h in zip(r['sel'], row['rs']):
                vals = []; t2 = []
                for sd in (2, 3, 4, 5, 11, 12345):
                    v, loc = ct08.pairing(A, C, dl, b[0], sqf(b[0] * b[1]), seed=sd)
                    vals.append(v); t2.append(loc.get(2))
                c[(r['cv'], 'global const' if len(set(vals)) == 1 else 'GLOBAL SEED-DEP', 'agree rs' if set(vals) == {h} else 'NO', '2-adic term const' if len(set(t2)) == 1 else 'T2 SEED-DEP', 'None' if None in vals else '')] += 1
                if len(set(t2)) > 1 and len(ex) < 3: ex.append((m, n, r['cv'], dl, b, vals, t2))
print(dict(c)); print(ex)
