#!/usr/bin/env python3
"""Checks that need no Sage or PARI; run in CI.

1. Unit tests of classify_t, the exact test applied to every candidate t in dm.py, dm2.py, dm3.py.
   The function is taken from the source of dm.py without importing Sage.
2. Consistency of dm_evidence.txt, dm_closures.txt and open_pairs_1000.txt: 157 pairs, 131 closed
   (91 E04, 25 D1, 6 E0L, and 9 with four points on E0L, section 6), 26 open, no overlap; every recorded t is re-tested and is degenerate; the listed bound
   amax is at least the bound that follows from the recorded heights.
3. Every recorded E04 and D1 generator is a rational point of its curve with y != 0.
   Every recorded E0L model is isomorphic over Q to the Weierstrass curve of G1.e0l_weierstrass for the
   pair, and the recorded generator is a rational point of that model, not of order 2.
4. Every line has an answer 'rank <= 1' in magma_check.txt (required first of all for the lines marked upper=magma).
5. dm_log.jsonl (JSON records assembled from the scripts' outputs, one per line; a partner line repeats the record of
   its run pair) agrees with dm_evidence.txt in the pair, the run pair, the curve name, the x-coordinate of the
   generator, the bound amax, and the height to 0.001.
6. dm4_log.jsonl (the JSON records printed by dm4.py, the argument with four points on E0L for rank 2 and 3): the pair
   is one of the 157 and is not closed otherwise; every recorded t is degenerate; the recorded model is the
   Weierstrass curve of the pair and the recorded generators are rational points of the minimal model, with a
   positive definite height matrix; the elementary constants lo4, up4, wpm, winv are recomputed from m and n; eL, eU,
   kappa and cstar are recomputed from them and from the recorded BUP and U; Hstar satisfies the condition of the
   argument; the saturation bound covers the index bound that follows from the recorded regulator and height
   lower bound. dm4_cross.jsonl: the same checks for 43 pairs of dm_evidence.txt that dm4.py closes a second time.
This does not check the heights, the ranks, or that the points generate."""
import ast, math, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import satcert
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
# --- 6. the four-point closures
def res2(A, B):                                  # resultant of two quadratics, coefficients from the leading one down
    a2, a1, a0 = A; b2, b1, b0 = B
    return (a2 * b0 - a0 * b2) ** 2 - (a2 * b1 - a1 * b2) * (a1 * b0 - a0 * b1)
def e0l_c4c6(m, n):                              # c4, c6 of the Weierstrass curve of G1.e0l_weierstrass
    mm, nn = Fr(m), Fr(n); S = mm * mm + nn * nn
    qa = 2 * nn ** 2 * S; qb = 8 * mm ** 3 * nn; qc = 2 * (7 * mm ** 4 - 4 * mm ** 2 * nn ** 2 + nn ** 4)
    qd = 4 * mm * (mm * mm - nn * nn) * (3 * mm * mm - nn * nn) / nn; qq = 2 * mm * (mm * mm - nn * nn) / nn
    a2w = qc - qd ** 2 / (4 * qq ** 2)
    return cinv(qd / qq, a2w, 2 * qq * qb, -4 * qq ** 2 * qa, a2w * (-4 * qq ** 2 * qa))
def same_curve(c46a, c46b):
    (c4a, c6a), (c4b, c6b) = c46a, c46b
    if c4a == 0 or c6a == 0 or c4b == 0 or c6b == 0: return False
    u2 = (c6a / c6b) / (c4a / c4b)
    return issq(u2) and u2 * u2 == c4a / c4b
def isqrt_up(q):
    """An integer >= sqrt(q), q a non-negative fraction."""
    v = math.isqrt(q.numerator // q.denominator + 1) + 1
    return v
def count_vectors(M, thr, s):
    """The number of integer vectors c with u^T M u <= thr, u = 2c - s (M a positive definite integer matrix).
    Exact: M = L^T D L by rational elimination, u^T M u = sum_i d_i (u_i + sum_{k>i} l_ik u_k)^2, coordinates from the last."""
    r = len(M); A = [[Fr(v) for v in row] for row in M]; d = [None] * r; l = [[Fr(0)] * r for _ in range(r)]
    for i in range(r):
        d[i] = A[i][i]; assert d[i] > 0
        for k in range(i + 1, r): l[i][k] = A[i][k] / d[i]
        for a in range(i + 1, r):
            for b in range(i + 1, r): A[a][b] -= d[i] * l[i][a] * l[i][b]
    def rec(i, u, rem):
        if i < 0: return 1
        off = sum(l[i][k] * u[k] for k in range(i + 1, r)); w = isqrt_up(rem / d[i]); tot = 0
        lo = math.floor((-off - w + s[i]) / 2) - 1; hi = math.ceil((-off + w + s[i]) / 2) + 1
        for c in range(lo, hi + 1):
            u[i] = 2 * c - s[i]; t = d[i] * (u[i] + off) ** 2
            if t <= rem: tot += rec(i - 1, u, rem - t)
        return tot
    return rec(r - 1, [0] * r, Fr(thr))
def check4(j, new, fresh=False):
    """One record of dm4.py. new = True: a closure of the census (rank 2 or 3, pair not closed otherwise);
    new = False: a second closure of a pair that dm_evidence.txt already closes (rank 1, 2 or 3)."""
    m, n = j['pair']; of = (m, n); S = m * m + n * n; r = j['rank']
    assert j['status'] == 'closed' and of in all157 and j['rank_bounds'][1] == r and r in ((2, 3) if new else (1, 2, 3)), of
    assert fresh or ((of not in closed) if new else (of in closed)), of
    for tt, verdict in j['t'].items():
        assert verdict == 'degenerate' and classify_t(m, n, Fr(tt)) == 'degenerate', (of, tt)
    mdl = [Fr(z) for z in j['model']]; mn = [Fr(z) for z in j['min_model']]
    assert same_curve(e0l_c4c6(m, n), cinv(*mdl)) and same_curve(cinv(*mdl), cinv(*mn)), of
    assert mdl[0] == 0 and mdl[2] == 0 and mdl[4] == 0 and all(z.denominator == 1 for z in mdl), of      # Y^2 = x (x - al)(x - be), integral
    for x, y in j['G'] + [j['R0']]:
        x, y = Fr(x), Fr(y)
        assert y * y + mn[0] * x * y + mn[2] * y == x ** 3 + mn[1] * x * x + mn[3] * x + mn[4], of
    g = j['gram']; assert len(g) == r and len(j['G']) == r
    det = g[0][0] if r == 1 else g[0][0] * g[1][1] - g[0][1] ** 2 if r == 2 else (g[0][0] * (g[1][1] * g[2][2] - g[1][2] ** 2) - g[0][1] * (g[0][1] * g[2][2] - g[1][2] * g[0][2]) + g[0][2] * (g[0][1] * g[1][2] - g[1][1] * g[0][2]))
    assert g[0][0] > 0 and (r == 1 or g[0][0] * g[1][1] - g[0][1] ** 2 > 0) and det > 0 and abs(det - j['reg']) < 1e-6 * det, of
    assert j['kmax'] >= math.sqrt(j['reg0'] * {1: 1.0, 2: 4.0 / 3, 3: 2.0}[r] / j['lam'] ** r) and j['lam'] > 0 and j['reg0'] >= j['reg'], of
    Pi = (2 * m * n ** 3, m ** 4 - n ** 4, -2 * m ** 3 * n); Psi0 = (n * n, 0, m * m); Psi = (S, -4 * m * n, S)
    r0, r1 = abs(res2(Pi, Psi0)), abs(res2(Pi, Psi)); assert [str(r0), str(r1)] == j['res'] and r0 > 0 and r1 > 0, of
    octic = m ** 8 + 14 * m ** 4 * n ** 4 + n ** 8
    assert r0 == m * m * n * n * octic and r1 == (m * m - n * n) ** 2 * octic, of      # the forms of G1.fourPoint_res_psi0_dvd, G1.fourPoint_res_psi_dvd
    lo4 = math.log(n * n * (m - n) ** 2) - math.log(2) - 2 * math.log(r0) - 2 * math.log(r1)
    up4 = math.log(2) + 2 * math.log(2 * m ** 3 * n + m ** 4 - n ** 4)
    wpm = math.log(32 * m * m * n * n * S) + math.log(2 * m * m)
    winv = math.log(4 * (m * m - n * n) ** 2 * (m + n) ** 2) + math.log(8 * (m + n) ** 2)
    BUP, U = j['BUP'], j['U']; eL = (BUP - lo4) / 4; eU = (up4 + U) / 4; kap = (BUP + U + max(wpm, winv)) / 4
    for name, v in (('lo4', lo4), ('up4', up4), ('wpm', wpm), ('winv', winv), ('eL', eL), ('eU', eU), ('kappa', kap)):
        assert abs(v - j[name]) < 1e-9 * max(1, abs(v)), (of, name)
    H = j['Hstar']; L = H - eL; ee = eL + eU
    # the bounds of G1.fourPoint_rank2_bound and G1.fourPoint_rank3_bound (FOUR_POINTS.md, Propositions 2 and 3)
    if r <= 2:
        assert L > 0, of
        cs = (ee * ee / 2 + 4 * kap * kap) / L
    else:
        assert L > 2 * kap, of
        cs = 4 * kap + (ee * ee + 16 * kap * ee) / (4 * (L - kap))
    assert abs(cs - j['cstar']) < 1e-9 * cs and j['solver_selftest'].startswith('ok'), of
    gc = j['gram_check']; assert 0 < gc['mu'] and 0 <= gc['relative'] <= 0.01 and j['arithmetic'].startswith('exact'), of
    assert abs(gc['relative'] - r * gc['max_deviation'] / gc['mu']) <= 1e-6 * gc['relative'] + 1e-15, of
    assert -eL <= j['observed N(a) - H'][0] and j['observed N(a) - H'][1] <= eU, of
    # -- the basis: bound on the index, and the saturation certificate (satcert.py; nothing from PARI or Sage is used here)
    herm = {1: 1.0, 2: 4.0 / 3, 3: 2.0}[r]
    assert 0 < j['lam_used'] <= 0.99 * j['lam'] * (1 + 1e-12) and j['kmax'] ** 2 * j['lam_used'] ** r >= j['reg0'] * herm * (1 - 1e-12), of
    al, be = (int(v) for v in j['roots']); assert mdl == [0, -(al + be), 0, al * be, 0] and al != be and al * be != 0, of
    gT = [(Fr(a), Fr(b)) for a, b in j['G_T']]; tT = [(Fr(a), Fr(b)) for a, b in j['tors_T']]
    assert len(gT) == r and len(tT) == j['ntors'] - 1 and len(set(tT)) == len(tT) and all(satcert.on_curve(P, al, be) for P in gT + tT + [(Fr(j['R0_T'][0]), Fr(j['R0_T'][1]))]), of
    assert all((Fr(e), Fr(0)) in tT for e in (0, al, be)), of
    g = 0                                                # #E(Q)_tors <= 16 (Mazur), so for a good prime q >= 17 its order is prime
    for q in satcert.primes_upto(400):                   # to q, and torsion of order prime to q injects into E(F_q)
        if q >= 17 and satcert.good(al, be, q): g = math.gcd(g, satcert.count(al, be, q))
    assert g == j['ntors'], (of, g)                      # so tors_T is the whole torsion subgroup (without the origin)
    assert all(satcert.is_prime(int(p)) and all(satcert.good(al, be, q) for q in qs) for p, qs in j['sat'].items()), of
    assert satcert.verify(j['sat'], al, be, gT, tT, j['kmax']), of
    # -- the two lists: the integer vectors, counted again by a different enumeration (Fincke-Pohst, exact)
    G0i = [[int(v) for v in row] for row in j['G0_2e60']]
    assert all(abs(G0i[a][b] / 2.0 ** 60 - j['gram'][a][b]) <= j['gram_check']['max_deviation'] + 1e-9 for a in range(r) for b in range(r)), of
    for name, bound in (('A', H + eU), ('B', cs)):
        L4 = j['lists'][name]; thr = int(L4['threshold'])
        assert abs(L4['bound'] - bound) < 1e-9 * bound and thr >= 4 * bound / (1 - gc['relative']) * 2.0 ** 60 * (1 - 1e-9), (of, name)
        assert L4['centre2'] == (j['rho'] if name == 'A' else [0] * r), (of, name)
        assert count_vectors(G0i, thr, L4['centre2']) == L4['vectors'], (of, name)
    assert j['nA'] == j['lists']['A']['vectors'] * j['ntors'] and j['nz'] == j['lists']['B']['vectors'] * j['ntors'], of
    return of, r
log4 = [json.loads(l) for l in open(os.path.join(here, 'dm4_log.jsonl'))]
closed4 = set(); ranks4 = {2: 0, 3: 0}
for j in log4:
    of, r = check4(j, True); assert of not in closed4; closed4.add(of); ranks4[r] += 1
cross4 = [json.loads(l) for l in open(os.path.join(here, 'dm4_cross.jsonl'))]      # second closures of pairs closed on a curve of rank 1
crossed = set()
for j in cross4:
    of, r = check4(j, False); assert of not in crossed; crossed.add(of)
assert len(crossed) == 43
closed |= closed4
assert closed <= all157
cl = open(os.path.join(here, 'dm_closures.txt')).read().split('\n')
still = set()
for i, l in enumerate(cl):
    if l.startswith('# Still open'): still |= set(pairs(cl[i + 1]))
assert len(log4) == 9 and ranks4 == {2: 4, 3: 5}
assert len(still) == 26 and not (still & closed) and (still | closed) == all157
listed = set(pairs(l.split()[0])[0] for l in cl if l and not l.startswith('#') and len(l.split()) == 4)
assert listed == closed
assert set(pairs(l.split()[0])[0] for l in cl if l and not l.startswith('#') and len(l.split()) == 4 and l.split()[2] == 'E0L4') == closed4
if '--negative' in sys.argv:
    # Negative tests: a record that is altered in one place must be rejected by check4.
    import copy
    def must_fail(j, change, label):
        k = copy.deepcopy(j); change(k)
        try: check4(k, tuple(k['pair']) in closed4, fresh=True)
        except AssertionError: return
        raise SystemExit('negative test not caught: %s for %s' % (label, j['pair']))
    def double_first(k):
        al, be = (int(v) for v in k['roots']); x, y = Fr(k['G_T'][0][0]), Fr(k['G_T'][0][1]); a2, a4 = -(al + be), al * be
        l = (3 * x * x + 2 * a2 * x + a4) / (2 * y); x3 = l * l - a2 - 2 * x; k['G_T'][0] = [str(x3), str(l * (x - x3) - y)]
    def composite(k):
        p = sorted(k['sat'], key=int)[-1]; k['sat'][p] = [q * 3 for q in k['sat'][p]]
    ntest = 0
    for j in log4 + cross4[:6]:
        check4(j, tuple(j['pair']) in closed4, fresh=True)
        must_fail(j, double_first, 'first basis point doubled (index 2)')
        must_fail(j, lambda k: k['sat'].__setitem__('2', k['sat']['2'][:1]), 'certificate for 2 shortened')
        must_fail(j, lambda k: k['sat'].pop(sorted(k['sat'], key=int)[-1]), 'certificate for the largest prime removed')
        must_fail(j, composite, 'auxiliary moduli not prime')
        must_fail(j, lambda k: k.__setitem__('kmax', k['kmax'] - 2), 'index bound lowered')
        must_fail(j, lambda k: k.__setitem__('tors_T', k['tors_T'][:-1]), 'a torsion point missing')
        must_fail(j, lambda k: k['lists']['A'].__setitem__('vectors', k['lists']['A']['vectors'] - 1), 'vector count of list A')
        must_fail(j, lambda k: k['lists']['B'].__setitem__('threshold', str(int(k['lists']['B']['threshold']) // 2)), 'threshold of list B halved')
        must_fail(j, lambda k: k.__setitem__('lam_used', k['lam']), 'height lower bound without the margin')
        must_fail(j, lambda k: k.__setitem__('cstar', k['cstar'] * 0.99), 'cstar lowered')
        ntest += 10
    print('negative tests: %d altered records, all rejected' % ntest)
if '--dm4' in sys.argv:
    # CI: records just produced by dm4.py (a file of JSON lines) get the same checks as the stored ones and must agree
    # with them in everything that does not depend on the choice of basis.
    fresh = [json.loads(l) for l in open(sys.argv[sys.argv.index('--dm4') + 1]) if l.strip()]
    assert fresh, 'no fresh records'
    stored = dict((tuple(j['pair']), j) for j in log4 + cross4)
    for j in fresh:
        of = tuple(j['pair']); assert of in stored, of
        check4(j, of in closed4, fresh=True); k = stored[of]
        for name in ('status', 'rank', 'rank_bounds', 'model', 'min_model', 'res', 't', 'ntors'):
            assert j[name] == k[name], (of, name)
        for name in ('BUP', 'U', 'eL', 'eU', 'kappa'):
            assert abs(j[name] - k[name]) < 1e-9 * abs(k[name]), (of, name)
        assert abs(j['reg'] - k['reg']) < 1e-6 * k['reg'], of
    print('fresh dm4 records agree with the stored ones:', ' '.join('(%d,%d)' % tuple(j['pair']) for j in fresh))
print('check_evidence OK: classify_t tests pass; 157 = 122 closed on curves of rank 1 (91 E04, 25 D1, 6 E0L) + 9 closed with four points on E0L (a separate, unreviewed tier) + 26 open')
