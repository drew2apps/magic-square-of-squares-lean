"""z14 task 1: independent test, on ratios NOT in the pilot, of EXPLORATION_3 §1 on the mixed curve
    E_mix: y^2 = (x + 2C^2)(x + 2AB)(x - 2AB),   (A, B, C) = (m^2 - n^2, 2mn, m^2 + n^2),
the curve E_0L = E_0U = E_8L = E_8U of t10_pairct.py (checked here with t10_pairct.Pair on every ratio).
Reference: the repository's own complete 2-descent and pairing (rank_second.full_selmer2, t10_pairct.DFullX).
  (a) explore3's cassels.py (copied as e3_cassels.py, paths only changed) against DFull on every (a, b), a over Sel^2
      (NA_ALL elements drawn with a fixed seed, the basis always included), b over the basis; DFull's own value on
      basis x basis (upper triangle, symmetry checked on 3 entries), extended bilinearly; DFull direct on 2 non-basis a.
  (b) the Redei product prod_{p | W odd exp} ((x0 + y0 s)/p) against DFull, for pairs of classes (1,d,d), (1,beta,beta)
      of Sel^2 with x^2 - d y^2 = beta z^2 soluble; several primitive solutions and both signs of s per prime.
  (c) V = {d : (1,d,d) in Sel^2} against N = {d : (d, W)_v = 1 for all v}; the note's conic
      R+ s^2 - R- D^2 = 2AB d w^2 soluble (PARI qfsolve, a finder; solutions checked) against d in N.
Written by an AI model (reviewer z14); not read by a person.  Experiment.
usage: mixtest.py OUT.jsonl i0 i1    (ratios i0 <= i < i1 of ratios.txt)"""
import sys, os, json, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
ARGS = sys.argv[1:]; sys.argv = sys.argv[:1]
sys.path.insert(0, os.path.join(HERE, '..', 'conductor-t10'))
import t10_pairct as T
RS = T.RS; RI = T.RI
sys.path.insert(0, HERE)
import e3_cassels as E3
from fractions import Fraction as Fr
from math import gcd
PARI = T.PARI
PARI.allocatemem(2 * 10 ** 8, 8 * 10 ** 8)          # z14: smaller PARI stack (e3_cassels asks for 10^9; memory is shared)
NA_ALL = 24           # a: the basis and random elements of Sel^2, NA_ALL in all (all of Sel^2 if 2^dim - 1 <= NA_ALL)
NREDEI = 150          # pairs (d, beta) per curve for (b)
NSOL = 8              # primitive solutions per Redei pair

def sqf(x):
    x = Fr(x); return RI.sqfree_part(x.numerator * x.denominator)

def factor(n):
    n = abs(int(n))
    if n < 2: return {}
    f = PARI.factor(n); out = {int(q): int(e) for q, e in zip(f[0], f[1])}
    for q in out: assert PARI.isprime(q)
    return out

def triple(c):
    """(c1, c2, c1 c2 mod squares) for squarefree c1, c2 (cheap; explore3's triple() factors c1 c2 with PARI)"""
    a1, a2 = int(c[0]), int(c[1]); g = gcd(a1, a2)
    return (a1, a2, (a1 // g) * (a2 // g))

def hilb(a, b, p):
    return RS.hilbert(a, b, p)

def vec(c, gens):
    """(c1, c2) -> bit vector over gens (both components); None if a prime outside gens."""
    v = 0; L = len(gens)
    for jj, d in enumerate(c):
        d = sqf(d); rest = abs(d)
        for i, g in enumerate(gens):
            if (g == -1 and d < 0) or (g != -1 and d % g == 0): v |= 1 << (i + jj * L)
            if g != -1 and rest % g == 0: rest //= g
        if rest != 1: return None
    return v

def coords(v, basis_vecs):
    """coordinates of v in the F2-span of basis_vecs (list of ints); None if not in the span."""
    rows = [(b, 1 << i) for i, b in enumerate(basis_vecs)]
    piv = []
    for r, c in rows:
        for pr, pc in piv:
            if r ^ pr < r: r ^= pr; c ^= pc
        if r: piv.append((r, c)); piv.sort(reverse=True)
        else: raise AssertionError('basis not independent')
    c = 0
    for pr, pc in piv:
        if v ^ pr < v: v ^= pr; c ^= pc
    return c if v == 0 else None

def predict(M, ca, cb, k):
    s = 0
    for i in range(k):
        if ca >> i & 1:
            for j in range(k):
                if cb >> j & 1 and M[i][j] == -1: s ^= 1
    return -1 if s else 1

def primitive(v):
    den = 1
    for x in v: den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
    v = [int(Fr(x) * den) for x in v]; g = 0
    for x in v: g = gcd(g, x)
    return [x // g for x in v]

def conic_solutions(d, beta, nsol, rng):
    """primitive solutions of x^2 - d y^2 = beta z^2: PARI qfsolve (a finder, checked), then the second intersections
    with lines through it and random small points (all checked)."""
    s = PARI.qfsolve(PARI.matdiagonal([1, -d, -beta]))
    if str(PARI.type(s)) != 't_COL': return None
    P = primitive([Fr(str(x)) for x in s])
    Q = lambda v: v[0] * v[0] - d * v[1] * v[1] - beta * v[2] * v[2]
    Bf = lambda u, w: u[0] * w[0] - d * u[1] * w[1] - beta * u[2] * w[2]
    assert Q(P) == 0
    sols = [P]
    tries = 0
    while len(sols) < nsol and tries < 200:
        tries += 1
        R = [rng.randint(-9, 9) for _ in range(3)]
        if Q(R) == 0 or not any(R): continue
        t = Fr(-2 * Bf(P, R), Q(R))
        P2 = primitive([P[i] + t * R[i] for i in range(3)])
        assert Q(P2) == 0
        if not any(P2) or P2[2] == 0: continue
        if any(P2 == x or [-y for y in P2] == x for x in sols): continue
        sols.append(P2)
    return sols

def redei_value(sol, d, Wodd):
    """(value, zero-flag, per-prime sign-independence) for one primitive solution."""
    x0, y0, z0 = sol
    tot = 1; zero = False; signdep = False; hack = 1
    for p in Wodd:
        assert PARI.kronecker(d, p) == 1, 'd not a square mod p | W'
        s = int(PARI.lift(PARI.sqrt(PARI.Mod(d % p, p))))
        assert (s * s - d) % p == 0
        v1 = int(PARI.kronecker((x0 + y0 * s) % p, p)); v2 = int(PARI.kronecker((x0 - y0 * s) % p, p))
        if v1 == 0 or v2 == 0: zero = True; hv = v1 or v2
        else:
            hv = v1
            if v1 != v2: signdep = True
        tot *= v1; hack *= hv
    return (None if zero else tot), hack, zero, signdep

def run(m, n, out):
    t0 = time.time()
    rng = random.Random(1000003 * m + n)
    A, B, C = m * m - n * n, 2 * m * n, m * m + n * n
    Rp, Rm = C * C + A * B, C * C - A * B; W = Rp * Rm
    es = (-2 * C * C, -2 * A * B, 2 * A * B)
    rec = dict(m=m, n=n)
    # the model is the mixed curve of t10_pairct (translation only)
    rec['model_ok'] = []
    for j, k in (('0', 'L'), ('0', 'U'), ('8', 'L'), ('8', 'U')):
        Pr = T.Pair(m, n, j, k)
        rec['model_ok'].append(Pr.u == 1 and (Pr.es[1] - Pr.es[0], Pr.es[2] - Pr.es[0]) == (es[1] - es[0], es[2] - es[0]))
    sh = es[1]; esh = tuple(e - sh for e in es)
    bad, loc, sel, kum, enc, places = RS.full_selmer2(esh, random.Random(1))
    k = len(sel); gens = [-1] + bad
    bv = [vec(c, gens) for c in sel]
    rec['dim'] = k; rec['sel'] = sel
    # DFull on basis x basis
    D = [T.DFullX(esh, a) for a in sel]
    M = [[None] * k for _ in range(k)]
    nfail = 0
    for i in range(k):
        for j in range(i, k):
            h, per = T.safe_pair(D[i], sel[j])
            M[i][j] = M[j][i] = h
            if h is None: nfail += 1
    symchk = []
    for _ in range(3 if k > 1 else 0):
        i, j = sorted(rng.sample(range(k), 2))
        symchk.append(T.safe_pair(D[j], sel[i])[0] == M[i][j])
    rec['M_symcheck'] = symchk
    rec['M'] = M; rec['M_fail'] = nfail
    sys.stderr.write('  M done %.1f s (k=%d)\n' % (time.time() - t0, k))
    rec['M_sym'] = all(M[i][j] == M[j][i] for i in range(k) for j in range(k))
    rec['M_diag'] = all(M[i][i] == 1 for i in range(k))
    # all elements of Sel^2 (index = coordinate mask)
    def elt(mask):
        c1 = c2 = 1
        for i in range(k):
            if mask >> i & 1: c1 *= sel[i][0]; c2 *= sel[i][1]
        return (sqf(c1), sqf(c2))
    if 2 ** k - 1 <= NA_ALL: masks = list(range(1, 2 ** k))
    else:
        masks = set(1 << i for i in range(k))
        while len(masks) < NA_ALL: masks.add(rng.randrange(1, 2 ** k))
        masks = sorted(masks)
    rec['na'] = len(masks)
    # (a) cassels.py
    cnt = {}; dis = []
    for mask in masks:
        a = elt(mask)
        if a[0] < 0: cnt['a not real'] = cnt.get('a not real', 0) + 1; continue
        try: forms = E3.tangent_forms(es, triple(a))
        except Exception as ex: cnt['forms fail'] = cnt.get('forms fail', 0) + 1; continue
        for j in range(k):
            b = sel[j]
            if b[0] < 0 and False: pass
            want = predict(M, mask, 1 << j, k)
            try: v, info = E3.cassels(es, triple(a), triple(b), seed=1, forms=forms)
            except Exception as ex: v, info = None, {'fail': repr(ex)}
            key = 'agree' if v == want else ('fail' if v is None else 'DISAGREE')
            cnt[key] = cnt.get(key, 0) + 1
            if key == 'DISAGREE' and len(dis) < 5: dis.append((a, b, v, want))
            if key == 'fail' and len(dis) < 5: dis.append((a, b, 'fail', str(info)[:200]))
    rec['cassels'] = cnt; rec['cassels_bad'] = dis
    sys.stderr.write('  cassels done %.1f s\n' % (time.time() - t0))
    # cassels.py with a second seed on 40 (a, b)
    cnt2 = {}
    for _ in range(12):
        mask = rng.choice(masks); a = elt(mask); j = rng.randrange(k)
        if a[0] < 0: continue
        v1 = E3.cassels(es, triple(a), triple(sel[j]), seed=1)[0]
        v2 = E3.cassels(es, triple(a), triple(sel[j]), seed=77)[0]
        cnt2[str((v1 == v2, v2 == predict(M, mask, 1 << j, k)))] = cnt2.get(str((v1 == v2, v2 == predict(M, mask, 1 << j, k))), 0) + 1
    rec['cassels_seed77'] = cnt2
    # DFull direct on 4 non-basis a (bilinearity of DFull's evaluation)
    nb = [x for x in masks if bin(x).count('1') > 1]
    cntd = {}
    for mask in (rng.sample(nb, min(2, len(nb))) if nb else []):
        a = elt(mask); Da = T.DFullX(esh, a)
        for j in range(k):
            h, per = T.safe_pair(Da, sel[j])
            key = 'agree' if h == predict(M, mask, 1 << j, k) else ('fail' if h is None else 'DISAGREE')
            cntd[key] = cntd.get(key, 0) + 1
    rec['dfull_bilinear'] = cntd
    sys.stderr.write('  dfull bilinear done %.1f s\n' % (time.time() - t0))
    # (b), (c): V = {d : (1, d, d) in Sel^2}
    Vmask = []
    for mask in range(2 ** k):
        if mask == 0 or elt(mask)[0] == 1: Vmask.append(mask)
    V = [elt(x)[1] for x in Vmask]
    rec['dimV'] = len(Vmask).bit_length() - 1
    assert 2 ** rec['dimV'] == len(Vmask)
    fW = factor(W)
    Wodd = sorted(p for p, e in fW.items() if e % 2 and p != 2)
    rec['W_odd_primes_mod8'] = sorted(set(p % 24 for p in Wodd))
    rec['W_has_2_odd'] = fW.get(2, 0) % 2 == 1
    PW = sorted(set([2]) | set(fW) | set(factor(2 * A * B)))
    def isnorm(d): return all(hilb(d, W, p) == 1 for p in PW) and hilb(d, W, 'inf') == 1
    rec['V_in_N'] = all(isnorm(d) for d in V)
    # N over <-1, primes of 2AB>
    gAB = [-1] + sorted(factor(2 * A * B))
    Nmasks = 0; Ntot = 2 ** len(gAB); conic = {}
    dl = []
    for mask in range(Ntot):
        d = 1
        for i, g in enumerate(gAB):
            if mask >> i & 1: d *= g
        nd = isnorm(d)
        Nmasks += nd
        dl.append((d, nd))
    rec['dimN_AB'] = Nmasks.bit_length() - 1; assert 2 ** rec['dimN_AB'] == Nmasks
    for d, nd in rng.sample(dl, min(24, len(dl))):
        s = PARI.qfsolve(PARI.matdiagonal([Rp, -Rm, -2 * A * B * d]))
        sol = str(PARI.type(s)) == 't_COL'
        if sol:
            P = primitive([Fr(str(x)) for x in s]); assert Rp * P[0] ** 2 - Rm * P[1] ** 2 == 2 * A * B * d * P[2] ** 2
        key = str((nd, sol)); conic[key] = conic.get(key, 0) + 1
    rec['conic_vs_norm'] = conic
    sys.stderr.write('  norm done %.1f s\n' % (time.time() - t0))
    # V inside the group generated by -1 and the primes of 2AB?
    rec['V_supported_on_2AB'] = all(all(q in gAB for q in factor(d)) for d in V)
    # (b) Redei
    pairs = [(x, y) for x in Vmask for y in Vmask if x and y]
    if len(pairs) > NREDEI: pairs = rng.sample(pairs, NREDEI)
    rb = {}; rbad = []
    for x, y in sorted(pairs):
        d, beta = elt(x)[1], elt(y)[1]
        want = predict(M, x, y, k)
        sols = conic_solutions(d, beta, NSOL, rng)
        if sols is None:
            key = 'insoluble'; rb[key] = rb.get(key, 0) + 1
            rb['insoluble, pairing %d' % want] = rb.get('insoluble, pairing %d' % want, 0) + 1
            continue
        vals = [redei_value(s_, d, Wodd) for s_ in sols]
        clean = set(v[0] for v in vals if v[0] is not None)
        hacks = set(v[1] for v in vals)
        nz = sum(v[2] for v in vals); sd = any(v[3] for v in vals)
        rb['soluble'] = rb.get('soluble', 0) + 1
        rb['solutions'] = rb.get('solutions', 0) + len(sols)
        if nz: rb['solutions with p | x0 + y0 s'] = rb.get('solutions with p | x0 + y0 s', 0) + nz
        if sd: rb['sign of s matters'] = rb.get('sign of s matters', 0) + 1
        if len(clean) > 1: rb['depends on the solution'] = rb.get('depends on the solution', 0) + 1
        if len(hacks) > 1: rb['zero-case rule depends on the solution'] = rb.get('zero-case rule depends on the solution', 0) + 1
        ok = clean == {want} and hacks == {want}
        rb['agree' if ok else 'DISAGREE'] = rb.get('agree' if ok else 'DISAGREE', 0) + 1
        rb['pairing %d' % want] = rb.get('pairing %d' % want, 0) + 1
        if not ok and len(rbad) < 5: rbad.append((d, beta, want, sorted(clean), sorted(hacks), sols[:3]))
    rec['redei'] = rb; rec['redei_bad'] = rbad
    # DFull directly on 3 Redei pairs (not through bilinearity)
    cdr = {}
    for x, y in (rng.sample(pairs, min(3, len(pairs))) if pairs else []):
        h, per = T.safe_pair(T.DFullX(esh, elt(x)), elt(y))
        key = 'agree' if h == predict(M, x, y, k) else 'DISAGREE'
        cdr[key] = cdr.get(key, 0) + 1
    rec['redei_direct'] = cdr
    rec['time'] = round(time.time() - t0, 1)
    with open(out, 'a') as f: f.write(json.dumps(rec, sort_keys=True) + '\n')
    return rec

if __name__ == '__main__':
    out, i0, i1 = ARGS[0], int(ARGS[1]), int(ARGS[2])
    lines = [l.split() for l in open(os.path.join(HERE, 'ratios.txt')) if l.strip() and l[0] != '#']
    done = set()
    if os.path.exists(out):
        for l in open(out): r = json.loads(l); done.add((r['m'], r['n']))
    for l in lines[i0:i1]:
        m, n = map(int, l[0].split(','))
        if (m, n) in done: continue
        r = run(m, n, out)
        sys.stderr.write('%d,%d dim %d dimV %d cassels %s redei %s time %.1f\n' % (m, n, r['dim'], r['dimV'], r['cassels'], r['redei'], r['time']))
