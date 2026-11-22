#!/usr/bin/env python3
"""Checks that need no Sage or PARI; run in CI.

1. Unit tests of classify_t, the exact test applied to every candidate t in dm.py, dm2.py, dm3.py.
   The function is taken from the source of dm.py without importing Sage.
2. Consistency of dm_evidence.txt, dm_closures.txt and open_pairs_1000.txt: 157 pairs, 122 closed
   (91 E04, 25 D1, 6 E0L), 35 open, no overlap; every recorded t is re-tested and is degenerate; the listed bound
   amax is at least the bound that follows from the recorded heights.
3. Every recorded E04 and D1 generator is a rational point of its curve with y != 0.
   Every recorded E0L model is isomorphic over Q to the Weierstrass curve of G1.e0l_weierstrass for the
   pair, and the recorded generator is a rational point of that model, not of order 2.
4. Every line has an answer 'rank <= 1' in magma_check.txt (required first of all for the lines marked upper=magma).
5. dm_log.jsonl (JSON records assembled from the scripts' outputs, one per line; a partner line repeats the record of
   its run pair) agrees with dm_evidence.txt in the pair, the run pair, the curve name, the x-coordinate of the
   generator, the bound amax, and the height to 0.001.
This does not check the heights, the ranks, or that the points generate."""
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
assert len(ev) == 122
closed = set(); meth = {'E04': 0, 'D1': 0, 'E0L': 0}
def cinv(a1, a2, a3, a4, a6):
    b2 = a1 * a1 + 4 * a2; b4 = 2 * a4 + a1 * a3; b6 = a3 * a3 + 4 * a6
    return b2 * b2 - 24 * b4, -b2 ** 3 + 36 * b2 * b4 - 216 * b6
magma = set(tuple(l.split()[:3]) for l in open(os.path.join(here, 'magma_check.txt')) if l.rstrip().endswith('rank <= 1'))
def partner(m, n): return (m + n, m - n) if (m + n) % 2 else ((m + n) // 2, (m - n) // 2)
for r in ev:
    of = pairs(r[0])[0]; run = pairs(r[1])[0]; f = dict(x.split('=', 1) for x in r[3:] if '=' in x)
    assert of not in closed; closed.add(of); meth[r[2]] += 1
    if f.get('upper') == 'magma': assert tuple(r[:3]) in magma, r[:3]
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
        mdl = [Fr(z) for z in f['model'].split(',')]; x = Fr(f['G']); mm, nn = Fr(m), Fr(n); S = mm * mm + nn * nn
        qa = 2 * nn ** 2 * S; qb = 8 * mm ** 3 * nn; qc = 2 * (7 * mm ** 4 - 4 * mm ** 2 * nn ** 2 + nn ** 4)
        qd = 4 * mm * (mm * mm - nn * nn) * (3 * mm * mm - nn * nn) / nn; qq = 2 * mm * (mm * mm - nn * nn) / nn
        a2w = qc - qd ** 2 / (4 * qq ** 2)
        c4w, c6w = cinv(qd / qq, a2w, 2 * qq * qb, -4 * qq ** 2 * qa, a2w * (-4 * qq ** 2 * qa)); c4m, c6m = cinv(*mdl)
        assert c4m != 0 and c6m != 0, r[:3]
        u2 = (c6w / c6m) / (c4w / c4m)                 # c4w = u^4 c4m, c6w = u^6 c6m with u rational: isomorphic over Q
        assert issq(u2) and u2 * u2 == c4w / c4m, r[:3]
        disc = (mdl[0] * x + mdl[2]) ** 2 + 4 * (x ** 3 + mdl[1] * x * x + mdl[3] * x + mdl[4])
        assert issq(disc) and disc != 0, r[:3]         # G is a rational point of the model, not of order 2
assert meth == {'E04': 91, 'D1': 25, 'E0L': 6}, meth
assert all(tuple(r[:3]) in magma for r in ev) and len(magma) == 122   # Magma's bound 'rank <= 1' is recorded for every line, and for nothing else
import json
log = [json.loads(l) for l in open(os.path.join(here, 'dm_log.jsonl'))]   # the scripts' JSON output: one record per line of dm_evidence.txt, same order
assert len(log) == 122
for r, j in zip(ev, log):
    f = dict(x.split('=', 1) for x in r[3:] if '=' in x); g = j['G'][0] if isinstance(j['G'], list) else j['G']
    assert [j['of'], j['run'], j['method']] == r[:3] and str(g) == f['G'] and int(j['amax']) == int(f['amax']) and j['status'] == 'closed' and '%d,%d' % tuple(j['pair']) == r[1], r[:3]
    assert abs(float(j['hG']) - float(f['hG'])) < 1e-3 and (('upper' in j) == (f.get('upper') == 'magma')), r[:3]
assert closed <= all157
cl = open(os.path.join(here, 'dm_closures.txt')).read().split('\n')
still = set()
for i, l in enumerate(cl):
    if l.startswith('# Still open'): still |= set(pairs(cl[i + 1]))
assert len(still) == 35 and not (still & closed) and (still | closed) == all157
listed = set(pairs(l.split()[0])[0] for l in cl if l and not l.startswith('#') and len(l.split()) == 4)
assert listed == closed
print('check_evidence OK: classify_t tests pass; 157 = 122 closed (91 E04, 25 D1, 6 E0L) + 35 open')
