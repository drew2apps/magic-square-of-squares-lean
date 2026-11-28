import os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__))  # repository copy: paths relative to this folder
import sys, os, collections
sys.path.insert(0, _os.path.join(HERE, '..', 'research-pairing')); os.chdir(_os.path.join(HERE, '..', 'conductor-t10'))
from parse import load, sqf
from ct08 import pairing
R = load(); idx = {'0': 0, '8': 1, 'L': 2, 'U': 3}
for which in ('08', 'LU'):
    c = collections.Counter(); dist = collections.Counter(); seen = set()
    for cv in R:
        if cv['cv'] != which: continue
        m, n = cv['ratio']; S = m * m + n * n
        for l in cv['lines']:
            if not l['b']: continue
            D = l['tup'][idx[which[0]]] * l['tup'][idx[which[1]]]
            if which == '08': A, C = m * m - n * n, S; dl = sqf(m * n * D)
            else: A, C = 8 * m * n, 4 * S; dl = sqf((m * m - n * n) * D)
            b_r, b_m = l['b']; b_s = sqf(b_r * b_m)
            v1, _ = pairing(A, C, dl, b_r, b_s, conv=1)
            c[('conv1 value', v1)] += 1; c[('b_r=1', b_r == 1)] += 1
            key = ((m, n), l['cls'], l['b'])
            if key not in seen:
                seen.add(key); dist[('conv1 value', v1)] += 1; dist[('b_r=1', b_r == 1)] += 1
    print(which, 'per kill line:', dict(c), ' per distinct (ratio,class,b):', dict(dist))
