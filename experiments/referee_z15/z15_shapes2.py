"""z15: (a) the tuples that survive E_08 and E_LU, by shape; (b) per ratio, does the one E_08 class of the S1 tuples
(resp. the one E_LU class of the S2 tuples) survive that curve; (c) how many ghost tuples per ghost ratio.
From the session's out_*.jsonl (copied to session_data/).  Written by an AI model (reviewer z15), not read by a person."""
import os, json, collections
import cypari2
P = cypari2.Pari(); HERE = os.path.dirname(os.path.abspath(__file__))
def sqf(x): return int(P.core(abs(int(x))))
CV = ['08', 'LU', '0L', '0U', '8L', '8U']
rows = [json.loads(l) for f in ('out_997.jsonl', 'out_15000_0.jsonl', 'out_15000_1.jsonl') for l in open(os.path.join(HERE, 'session_data', f))]
def shape(m, n, t):
    A = m * m - n * n; a = sqf(m * n * t[0] * t[1]); b = sqf(A * t[2] * t[3])
    if (a, b) == (sqf(2 * m * n), 1): return 'S1'
    if (a, b) == (1, sqf(2 * A)): return 'S2'
    return 'a=b=1' if a == b == 1 else 'a=1' if a == 1 else 'b=1' if b == 1 else 'other'
c08lu = collections.Counter(); ratio_s = collections.Counter(); per_ratio = collections.Counter()
for r in rows:
    m, n = r['m'], r['n']; g = collections.defaultdict(list)
    for t, d in r['tuples']:
        s = shape(m, n, t); k = [c for c in CV if d[c] == 1]
        if '08' not in k and 'LU' not in k: c08lu[(s, 'killed on a mixed curve' if k else 'SURVIVES')] += 1
        if s in ('S1', 'S2'): g[s].append(d)
    for s, ds in g.items():
        cur = '08' if s == 'S1' else 'LU'
        ratio_s[(s, 'class survives E_%s' % cur if ds[0][cur] == 0 else 'class killed on E_%s' % cur,
                 'all its tuples killed' if all(any(d[c] == 1 for c in CV) for d in ds) else 'a tuple survives')] += 1
    if g: per_ratio[sum(len(v) for v in g.values())] += 1
print('z15_shapes2.py: written by an AI model (reviewer z15), not read by a person.')
print('(a) tuples surviving E_08 and E_LU, by shape:')
for k, v in sorted(c08lu.items()): print('   ', k, v)
print('(b) per ratio and shape:')
for k, v in sorted(ratio_s.items()): print('   ', k, v)
print('(c) ghost tuples per ghost ratio:', sorted(per_ratio.items()))
