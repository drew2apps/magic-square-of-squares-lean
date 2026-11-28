"""z15: is "every survivor has the ghost shape" informative?  Per ratio: residual size, number of ghost-shape tuples,
survivors.  Written by an AI model (reviewer z15), not read by a person.  From session_data/out_*.jsonl."""
import os, json, collections
import cypari2
P = cypari2.Pari(); HERE = os.path.dirname(os.path.abspath(__file__))
def sqf(x): return int(P.core(abs(int(x))))
CV = ['08', 'LU', '0L', '0U', '8L', '8U']
rows = [json.loads(l) for f in ('out_997.jsonl', 'out_15000_0.jsonl', 'out_15000_1.jsonl') for l in open(os.path.join(HERE, 'session_data', f))]
def ghost(m, n, t):
    A = m * m - n * n; a = sqf(m * n * t[0] * t[1]); b = sqf(A * t[2] * t[3])
    return (a, b) == (sqf(2 * m * n), 1) or (a, b) == (1, sqf(2 * A))
kind = collections.Counter(); kind_s = collections.Counter(); detail = []
for r in rows:
    m, n = r['m'], r['n']; N = len(r['tuples'])
    g = sum(ghost(m, n, t) for t, d in r['tuples'])
    s = sum(1 for t, d in r['tuples'] if not any(d[c] == 1 for c in CV))
    k = 'all ghost' if g == N else 'no ghost' if g == 0 else 'mixed'
    kind[k] += 1; kind_s[k] += s > 0
    if s: detail.append((m, n, N, g, s))
print('z15_shapes3.py: written by an AI model (reviewer z15), not read by a person.')
print('ratios by composition of the residual set:', dict(kind), '; of them with a survivor:', dict(kind_s))
print('survivor ratios (m, n, residual size, ghost tuples, survivors):', detail)
print('residual size 2, both ghost:', sum(1 for r in rows if len(r['tuples']) == 2 and all(ghost(r['m'], r['n'], t) for t, d in r['tuples'])))
# ratio-level: which tuples reach the mixed curves (survive E_08 and E_LU)?
reach = collections.Counter(); reach_s = collections.Counter()
for r in rows:
    m, n = r['m'], r['n']
    R = [(ghost(m, n, t), d) for t, d in r['tuples'] if d['08'] != 1 and d['LU'] != 1]
    if not R: k = 'no tuple reaches the mixed curves'
    else:
        gh = any(g for g, d in R); ng = any(not g for g, d in R)
        k = 'ghost and non-ghost reach' if gh and ng else 'only ghost reach' if gh else 'only non-ghost reach'
    reach[k] += 1
    reach_s[k] += any(not any(d[c] == 1 for c in CV) for t, d in r['tuples'])
print('ratios by which tuples survive E_08 and E_LU:', dict(reach), '; with a survivor:', dict(reach_s))
