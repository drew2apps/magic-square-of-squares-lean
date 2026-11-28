import os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__))  # repository copy: paths relative to this folder
import json, glob, collections, sys, os
sys.path.insert(0, _os.path.join(HERE, '..', 'research-pairing')); os.chdir(_os.path.join(HERE, '..', 'conductor-t10'))
import ct08
from parse import sqf
from sympy import factorint
c = collections.Counter(); ex = []
for f in sorted(glob.glob(_os.path.join(HERE, 'full[01].jsonl'))):
    for l in open(f):
        r = json.loads(l)
        if 'error' in r: continue
        m, n, cv = r['m'], r['n'], r['cv']
        A, C = (m * m - n * n, m * m + n * n) if cv == '08' else (8 * m * n, 4 * (m * m + n * n))
        for row in r['rows']:
            if not row['ct']: continue
            dl = row['ap'][1]
            for b, h in zip(r['sel'], row['rs']):
                if b[0] != 1: continue
                beta = b[1]
                v, loc = ct08.pairing(A, C, dl, 1, beta)
                t = loc[2]
                for p in factorint(abs(C + A)):
                    if p != 2: t *= ct08.hil(dl, beta, p)
                for p in factorint(abs(A)):
                    if p != 2 and not ct08.issq(ct08.Fr(dl), p): t *= ct08.hil(A * C, beta, p)
                c[(cv, t == h)] += 1
                if t != h and len(ex) < 3: ex.append((m, n, cv, dl, b, h, t))
print('b_r = 1 product formula vs rank_second (distinct (ratio, a, b), b a basis element with b_r = 1):', dict(c), ex)
