"""z15: the ghost shape, cosets, set sizes and failure rates of EXPLORATION_4 sections 2-3, recomputed with my own code
from (i) the session's out_*.jsonl (copied to session_data/) and (ii) the residual tuples recomputed with the
repository's code (z15_residual_all.jsonl).  Squarefree parts by PARI's core.  Written by an AI model (reviewer z15),
not read by a person.  usage: python3 z15_shapes.py > z15_shapes.txt"""
import os, json, collections
import cypari2
P = cypari2.Pari()
HERE = os.path.dirname(os.path.abspath(__file__))
def sqf(x): return int(P.core(abs(int(x))))
CV = ['08', 'LU', '0L', '0U', '8L', '8U']
rows = []
for f in ('out_997.jsonl', 'out_15000_0.jsonl', 'out_15000_1.jsonl'):
    for l in open(os.path.join(HERE, 'session_data', f)):
        r = json.loads(l); r['src'] = f; rows.append(r)
mine = {}
p = os.path.join(HERE, 'z15_residual_all.jsonl')
if os.path.exists(p):
    for l in open(p):
        r = json.loads(l); mine[(r['m'], r['n'])] = sorted(tuple(t) for t in r['tuples'])
print('z15_shapes.py: written by an AI model (reviewer z15), not read by a person.')
print('session ratios', len(rows), '; tuples', sum(len(r['tuples']) for r in rows),
      '; errors', sum('error' in r for r in rows))
if mine:
    same = sum(sorted(tuple(t) for t, d in r['tuples']) == mine.get((r['m'], r['n'])) for r in rows)
    print('residual sets equal to the repository-code recomputation (z15_residual_all.jsonl):', same, 'of', len(rows))

def killed(d): return [c for c in CV if d[c] == 1 or d[c] == 'outside']
def shape(m, n, t):
    A = m * m - n * n
    a = sqf(m * n * t[0] * t[1]); b = sqf(A * t[2] * t[3])
    if (a, b) == (sqf(2 * m * n), 1): return 'S1', a, b
    if (a, b) == (1, sqf(2 * A)): return 'S2', a, b
    if a == 1 and b == 1: return 'a=b=1', a, b
    if a == 1: return 'a=1', a, b
    if b == 1: return 'b=1', a, b
    return 'other', a, b

cnt = collections.Counter(); ghost_r = set(); ghost_n = 0; surv = []
breakdown = collections.Counter(); incons = []; torsion_viol = []; per_curve = collections.Counter()
outside = 0; two_trivial_mixed = collections.Counter()
for r in rows:
    m, n = r['m'], r['n']
    s1v = set(); s2v = set()
    for t, d in r['tuples']:
        s, a, b = shape(m, n, t)
        k = killed(d)
        for c in k: per_curve[c] += 1
        outside += any(d[c] == 'outside' for c in CV)
        cnt[(s, bool(k))] += 1
        # a = 1: the class on E_08 is trivial; b = 1: the class on E_LU is a torsion class; such a class is never killed there
        if a == 1 and d['08'] != 0: torsion_viol.append((m, n, t, '08'))
        if b == 1 and d['LU'] != 0: torsion_viol.append((m, n, t, 'LU'))
        if s in ('S1', 'S2'):
            ghost_n += 1; ghost_r.add((m, n))
            if not k: surv.append((m, n, t, s))
            else: breakdown['on E_08 or E_LU' if ('08' in k or 'LU' in k) else 'only on mixed curves'] += 1
            (s1v if s == 'S1' else s2v).add(d['08'] if s == 'S1' else d['LU'])
        elif not k: surv.append((m, n, t, s))
        # tuples that survive both E_08 and E_LU: what kills them
        if '08' not in k and 'LU' not in k: two_trivial_mixed[(s if s in ('S1', 'S2') else 'not ghost', bool(k))] += 1
    if len(s1v) > 1 or len(s2v) > 1: incons.append((m, n))
print()
print('shape counts (shape, killed on some curve):')
for k, v in sorted(cnt.items(), key=str): print('   ', k, v)
print('ghost-shape tuples (S1 or S2):', ghost_n, 'at', len(ghost_r), 'ratios; killed:', sum(breakdown.values()), dict(breakdown))
print('tuples killed on no curve:', len(surv), 'at', len(set((m, n) for m, n, t, s in surv)), 'ratios:')
for x in surv: print('   ', x)
print('ratios where the S1 tuples (resp. S2) do not all have the same E_08 (resp. E_LU) verdict (they have one class there):', incons)
print('a = 1 killed on E_08, or b = 1 killed on E_LU (impossible for a torsion class):', torsion_viol[:5], len(torsion_viol))
print('classes outside Sel2:', outside)
print('killed per curve:', {c: per_curve[c] for c in CV})
print('tuples surviving E_08 and E_LU, by (shape, killed on a mixed curve):', dict(two_trivial_mixed))

# 3. cosets: V is a coset of an F_2-subspace iff 2^rank(V - v0) = |V| (exact for every size)
def vec(t):
    v = set()
    for i, d in enumerate(t):
        f = P.factor(d)
        for q, e in zip(f[0], f[1]):
            if int(e) % 2: v.add((i, int(q)))
    return frozenset(v)
def f2rank(vs):
    idx = {}; rows_ = []
    for v in vs:
        x = 0
        for e in v:
            if e not in idx: idx[e] = len(idx)
            x |= 1 << idx[e]
        rows_.append(x)
    basis = []
    for x in rows_:
        for b in basis:
            if x ^ b < x: x ^= b
        if x:
            basis.append(x); basis.sort(reverse=True)
    return len(basis)
aff = collections.Counter(); sizes = collections.Counter(); odd = []; stab = collections.Counter(); bysrc = collections.Counter()
for r in rows:
    m, n = r['m'], r['n']
    T = [tuple(t) for t, d in r['tuples']]
    V = [vec(t) for t in T]; Vs = set(V)
    assert len(Vs) == len(V), 'two tuples with one vector'
    v0 = V[0]; rk = f2rank([v ^ v0 for v in V])
    isc = 2 ** rk == len(V)
    aff[isc] += 1; sizes[len(V)] += 1; bysrc[(r['src'] == 'out_997.jsonl', isc)] += 1
    if len(V) % 2: odd.append((m, n))
    # translation stabiliser of V (a non-trivial one forces |V| even)
    st = [w for w in set(v ^ v0 for v in V) if w and all((v ^ w) in Vs for v in V)]
    A = m * m - n * n; c1 = sqf(2 * A); c2 = sqf(2 * m * n)
    w1 = vec((1, 1, c1, c1)); w2 = vec((c2, c2, 1, 1))
    stab[('stabiliser dim>0', bool(st), 'contains (1,1,c,c) c=sqf(2(m2-n2))', w1 in st, 'contains (c,c,1,1) c=sqf(2mn)', w2 in st)] += 1
print()
print('coset of an F2-subspace (exact rank test, all sizes):', dict(aff), '; by list (m<=5000?, coset?):', dict(bysrc))
print('sizes:', sorted(sizes.items()))
print('odd sizes:', odd)
print('translation stabilisers:')
for k, v in sorted(stab.items(), key=str): print('   ', k, v)
# failure rates
n997 = sum(1 for r in rows if r['src'] == 'out_997.jsonl'); n15 = len(rows) - n997
s997 = len(set((m, n) for m, n, t, s in surv if m <= 5000)); s15 = len(set((m, n) for m, n, t, s in surv if m > 5000))
print('ratios with tuples after bad primes: m<=5000 %d, 5000<m<=15000 %d; with a tuple after the six pairings: %d, %d'
      % (n997, n15, s997, s15))
print('rates: %d/997 = %.2f%%, %d/8739 = %.2f%%, %d/997 = %.3f%%, %d/8739 = %.3f%%'
      % (n997, 100 * n997 / 997, n15, 100 * n15 / 8739, s997, 100 * s997 / 997, s15, 100 * s15 / 8739))
