#!/usr/bin/env python3
"""Checks that need no Sage or PARI; run in CI.

1. Unit tests of classify_t, the exact test applied to every candidate t in dm.py, dm2.py, dm3.py.
   The function is taken from the source of dm.py without importing Sage.
2. Consistency of dm_evidence.txt, dm_closures.txt and open_pairs_1000.txt: 157 pairs, 119 closed
   (91 E04, 22 D1, 6 E0L), 38 open, no overlap; every recorded t is degenerate; the listed bound
   amax is at least the bound that follows from the recorded heights.
This does not check the heights, generators or ranks themselves."""
import ast, math, os, re, sys
from fractions import Fraction as Fr
here = os.path.dirname(os.path.abspath(__file__))

src = open(os.path.join(here, 'dm.py'), encoding='utf-8').read()
ns = {'math': math, 'Fr': Fr}
for node in ast.parse(src).body:
    if isinstance(node, ast.FunctionDef) and node.name in ('issq_int', 'issq', 'sqrt_fr', 'classify_t'):
        exec(compile(ast.Module([node], []), 'dm.py', 'exec'), ns)
classify_t = ns['classify_t']; issq = ns['issq']

# --- 1. classify_t
assert issq(Fr(9, 4)) and not issq(Fr(2)) and not issq(Fr(-1)) and issq(0)
# degenerate squares: t = +-1 when m^2 + n^2 is a square, t = 0 when (m^2 + n^2)/2 is a square
for m, n in [(4, 3), (12, 5), (15, 8), (21, 20), (72, 65)]:
    assert classify_t(m, n, Fr(1)) == 'degenerate' and classify_t(m, n, Fr(-1)) == 'degenerate'
    assert classify_t(m, n, Fr(0)) == 'no'
for m, n in [(7, 1), (17, 7), (23, 7), (41, 1)]:
    assert classify_t(m, n, Fr(0)) == 'degenerate' and classify_t(m, n, Fr(1)) == 'no'
# ordinary values: never all five squares
cnt = 0
for m, n in [(2, 1), (3, 2), (15, 8), (231, 16)]:
    for p in range(-30, 31):
        for q in range(1, 12):
            if math.gcd(p, q) == 1 and classify_t(m, n, Fr(p, q)) != 'no': cnt += 1
assert cnt == 2, cnt                     # only (15, 8) at t = 1 and t = -1
# a value where four of the five forms are squares must still be 'no': (4,3) at t = 0 has F0, F8 squares
assert classify_t(4, 3, Fr(0)) == 'no'

# --- 2. the tables
def pairs(line): return [tuple(int(v) for v in s.split(',')) for s in line.split()]
lines = open(os.path.join(here, 'open_pairs_1000.txt')).read().split('\n')
all157 = set()
for l in lines:
    if l and not l.startswith('#'): all157 |= set(pairs(l))
assert len(all157) == 157
ev = [l.split() for l in open(os.path.join(here, 'dm_evidence.txt')) if l.strip() and not l.startswith('#')]
assert len(ev) == 119
closed = set(); meth = {'E04': 0, 'D1': 0, 'E0L': 0}
def partner(m, n): return (m + n, m - n) if (m + n) % 2 else ((m + n) // 2, (m - n) // 2)
for r in ev:
    of = pairs(r[0])[0]; run = pairs(r[1])[0]; f = dict(x.split('=', 1) for x in r[3:] if '=' in x)
    assert of not in closed; closed.add(of); meth[r[2]] += 1
    assert run == of or run == partner(*of) or partner(*run) == of, r[:2]
    hG = float(f['hG']); B = float(f['B']); amax = int(f['amax'])
    assert re.fullmatch(r'\[((-?\d+(/\d+)?:deg)(;-?\d+(/\d+)?:deg)*)?\]', f['t']), f['t']
    for item in f['t'][1:-1].split(';'):          # every recorded t is re-tested, on the pair the curve was run for
        if item: assert classify_t(run[0], run[1], Fr(item.split(':')[0])) == 'degenerate', (r[:2], item)
    if r[2] == 'E04':
        m, n = run; kappa = 2 * (m * m + n * n) * n * n
        assert amax >= (2 * B + 2 * math.log(kappa)) / hG, r[:3]
        x = Fr(f['G']); S = m * m + n * n         # the recorded generator is a rational point of E04 with y != 0
        assert issq(x * (x + 2 * S * m * m) * (x + 2 * S * n * n)) and x * (x + 2 * S * m * m) * (x + 2 * S * n * n) != 0, r[:3]
    elif r[2] == 'D1':
        assert amax >= 2 * B / hG, r[:3]
        m, n = run; S = m * m + n * n; D2 = (m * m - n * n) ** 2; x = Fr(f['G'])
        c3 = 4 * S * S * m * m * n * n; c2 = -16 * S * m ** 3 * n ** 3; c1 = 4 * S * S * D2; c0 = -16 * S * m * n * D2
        v = x ** 3 + c1 * x * x + c0 * c2 * x + c0 * c0 * c3     # the recorded generator is a rational point of D1 with y != 0
        assert issq(v) and v != 0, r[:3]
    else:
        K = float(f['K']); m, n = run
        W = (2 * K + 4 * (math.log(max(m, n)) + math.log(2))) / (2 * hG)
        assert amax >= W + abs(int(f['rho'])) / 2 and float(f['worst']) < K, r[:3]
assert meth == {'E04': 91, 'D1': 22, 'E0L': 6}, meth
assert closed <= all157
cl = open(os.path.join(here, 'dm_closures.txt')).read().split('\n')
still = set()
for i, l in enumerate(cl):
    if l.startswith('# Still open'): still |= set(pairs(cl[i + 1]))
assert len(still) == 38 and not (still & closed) and (still | closed) == all157
listed = set(pairs(l.split()[0])[0] for l in cl if l and not l.startswith('#') and len(l.split()) == 4)
assert listed == closed
print('check_evidence OK: classify_t tests pass; 157 = 119 closed (91 E04, 22 D1, 6 E0L) + 38 open')
