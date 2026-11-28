import os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__))  # repository copy: paths relative to this folder
import json, glob, collections
from sympy import factorint
def sqf(n):
    s = -1 if n < 0 else 1; r = 1
    for q, e in factorint(abs(n)).items():
        if e % 2: r *= q
    return s * r
idx = {'0': 0, '8': 1, 'L': 2, 'U': 3}
c = collections.Counter(); kinds = collections.Counter()
for f in sorted(glob.glob(_os.path.join(HERE, 'full[01].jsonl'))):
    for l in open(f):
        r = json.loads(l)
        if 'error' in r: continue
        m, n, cv = r['m'], r['n'], r['cv']
        for row in r['rows']:
            for w in row['why']:
                if w[0] != 'tuple': continue
                t = eval(w[1]); D = sqf(t[idx[cv[0]]] * t[idx[cv[1]]])
                a = tuple(row['a'])
                if cv == '08': ok = a[0] == 1 and abs(a[1]) == abs(sqf(m * n * D))
                else: ok = a[0] == sqf(2 * (m * m + n * n)) and abs(a[1]) == abs(sqf(2 * D)) and abs(row['ap'][1]) == abs(sqf((m * m - n * n) * D))
                c[(cv, ok)] += 1
                kinds[(cv, tuple(v for k, v in sorted(zip(['08','0L','0U','8L','8U','LU'], [None]*6))))] += 0
print('class formula of section 1 on every residual tuple of the sampled ratios (both signs):', dict(c))
