#!/usr/bin/env python3
"""F4: Demjanenko-Manin on C: Y^2 = g4(s) K(s) in exact arithmetic, given rank C(Q) = 1.
Written by an AI model (Claude); not read by a person.  Experiment, not counted.

    python3 g4k_exact.py                 the three ratios; checks the output equals g4k_exact.txt
    python3 g4k_exact.py --write         write g4k_exact.txt
    python3 g4k_exact.py --negative      controls that must be refused (with --write: g4k_exact_negative.txt)
    python3 g4k_exact.py m,n [m,n ...]   print only

PARI is used only for hints, each checked exactly here: the point G (ellrank, ellsaturation) and the change of
coordinates to a minimal model (ellminimalmodel).  Everything that the conclusion rests on is exact (Fractions and
integers), with the repository's own routines at 24d6ef6 where they exist:
  (1) torsion of order 2 (exact_last_pair.torsion: point counts at good primes);
  (2) an enclosure lo <= hhat(G) <= hi (exact_last_pair.enclose: Silverman's bounds, 7 exact doublings);
  (3) saturation of G at every prime <= 199 (exact_last_pair.saturation);
  (4) every point P of infinite order has hhat(P) >= hi / 210^2: on the minimal model (checked isomorphic and minimal
      at every bad prime), Q = m P lies in E0(Q_p) at every p with m = lcm(12, n_p over the multiplicative primes I_n);
      for such Q, h(x(Q)) - hhat(Q) <= log+|j|/4 + log+|b2/12| + log 2* + 1.946 (Cremona-Prickett-Siksek, Theorems 2
      and 15, as quoted in exact_last_pair.lam_exact; there for models with a1 = a3 = 0, here with b2 = a1^2 + 4 a2 of
      the minimal model), and a search of every x = a/e^2 below that height finds no point other than 2-torsion.
      With (3), IF rank C(Q) = 1 then G generates C(Q) modulo torsion;
  (5) the Demjanenko-Manin bound: a solution gives P, P' on C with x(P) x(P') = c3^2 (see dm_g4k.py), so
      |hhat(P) - hhat(P')| <= B + U + 2 log|c3| (B, U of (2), on the model the relation is written in).  With P = aG + T0
      and P' = bG + T1 (a, b >= 0), |a| != |b| gives 2 max(a, b) - 1 <= (B + U + 2 log|c3|) / lo; |a| = |b| gives
      P' = +-P + T with T torsion, and then x(P) is a root of x^2 = c3^2 or of e x^2 + (beta - e^2 - c3^2) x + c3^2 e = 0
      (beta = f'(e), from x(P + T) = e + beta / (x - e)).  Every aG + T0 with a up to the bound and every rational root
      is tested in the five forms (classify_t of global.py, exact).
So the conclusion is: IF rank C(Q) = 1, every solution t is degenerate.  Rank 1 is PARI's ellrank upper bound (it uses
the Cassels-Tate pairing); the exact bound of g4k_rank_exact.py is 2, 3 and 5 for the three curves."""
import sys, os, io, contextlib, math
from fractions import Fraction as Fr
from math import gcd, isqrt
HERE = os.path.dirname(os.path.abspath(__file__))
_IN_REPO = os.path.join(HERE, '..', '..', '..')
REPO = os.environ.get('MSS_REPO', _IN_REPO)
sys.path.insert(0, HERE)
import dm_g4k
sys.path.insert(0, os.path.join(REPO, 'scripts', 'g1_full')); sys.path.insert(0, os.path.join(REPO, 'scripts', 'second_descent'))
with contextlib.redirect_stdout(io.StringIO()):
    import exact_last_pair as X
    import sigma_exact as SE
    import rank_indep as RI
import cypari2
pari = cypari2.Pari()
RATIOS = ['7735,4968', '6557,324', '6881,6233']
LMAX = 199

def fr(z): return Fr(int(pari.numerator(z)), int(pari.denominator(z)))

def classify_t(m, n, t):
    """global.py's classify_t (exact)."""
    s = Fr(m * m + n * n, 2)
    F0 = m * m + n * n * t * t; F8 = n * n + m * m * t * t; F4 = s * (1 + t * t); L = F4 - 2 * m * n * t; U = F4 + 2 * m * n * t
    sq = lambda q: (lambda q: q >= 0 and isqrt(q.numerator) ** 2 == q.numerator and isqrt(q.denominator) ** 2 == q.denominator)(Fr(q))
    if not all(sq(v) for v in (F0, F8, F4, L, U)): return 'no'
    cells = [(m * t + n) ** 2, (m * t - n) ** 2, (n * t + m) ** 2, (n * t - m) ** 2, F0, F8, F4, L, U]
    return 'SOLUTION' if len(set(cells)) == 9 and all(c != 0 for c in cells) else 'degenerate'

def to_t(X_, c3):
    w = Fr(X_) / c3
    if w < 0 or isqrt(w.numerator) ** 2 != w.numerator or isqrt(w.denominator) ** 2 != w.denominator: return []
    r = Fr(isqrt(w.numerator), isqrt(w.denominator)); return [r, -r] if r else [r]

class GenEll:
    """A general integral Weierstrass model [a1, a2, a3, a4, a6] (invariants only)."""
    def __init__(s, a):
        s.a1, s.a2, s.a3, s.a4, s.a6 = a
        a1, a2, a3, a4, a6 = a
        s.b2 = a1 * a1 + 4 * a2; s.b4 = 2 * a4 + a1 * a3; s.b6 = a3 * a3 + 4 * a6
        s.b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
        s.c4 = s.b2 ** 2 - 24 * s.b4; s.c6 = -s.b2 ** 3 + 36 * s.b2 * s.b4 - 216 * s.b6
        s.D = -s.b2 ** 2 * s.b8 - 8 * s.b4 ** 3 - 27 * s.b6 ** 2 + 9 * s.b2 * s.b4 * s.b6

def change(a, u, r, s_, t):
    """The a-invariants after x = u^2 x' + r, y = u^3 y' + s u^2 x' + t (Silverman III.1)."""
    a1, a2, a3, a4, a6 = a
    return [Fr(a1 + 2 * s_, u), Fr(a2 - s_ * a1 + 3 * r - s_ * s_, u ** 2), Fr(a3 + r * a1 + 2 * t, u ** 3),
            Fr(a4 - s_ * a3 + 2 * r * a2 - (t + r * s_) * a1 + 3 * r * r - 2 * s_ * t, u ** 4),
            Fr(a6 + r * a4 + r * r * a2 + r ** 3 - t * a3 - t * t - r * t * a1, u ** 6)]

def minimal_at(E, p):
    """Sufficient for minimality at p: v(c4) < 4, v(c6) < 6 or v(D) < 12 (a change with p | u would divide c4 by p^4,
    c6 by p^6 and D by p^12).  Otherwise, for p = 2 or 3 only, no change with u = p and r, s, t modulo p^6 gives an
    integral model (r, s, t of such a change are p-integral: as in exact_last_pair.minimal_at, written there for
    a1 = a3 = 0 and p = 2)."""
    v = lambda x: RI.vp(x, p) if x else 99
    if v(E.c4) < 4 or v(E.c6) < 6 or v(E.D) < 12: return True
    assert p in (2, 3), ('minimality at %d needs another test' % p)
    a = [E.a1, E.a2, E.a3, E.a4, E.a6]; M = p ** 6
    for r in range(M):
        for s_ in range(p):
            for t in range(M):
                if all(c.denominator == 1 for c in change(a, p, r, s_, t)): return False
    return True

def lower_height(Em, need):
    """Every point of infinite order has hhat >= need, on the minimal model Em (see (4) above).  Returns the data."""
    bad = sorted(RI.factor_int(abs(Em.D)))
    for p in bad: assert minimal_at(Em, p), ('not minimal at', p)
    m = 12
    for p in bad:
        if Em.c4 % p: n_ = RI.vp(Em.D, p); m = m * n_ // gcd(m, n_)
    j = Fr(Em.c4 ** 3, Em.D)
    Binf = SE.lnplus_hi(abs(j)) / 4 + SE.lnplus_hi(Fr(abs(Em.b2), 12)) + (SE.LN2[1] if Em.b2 else 0) + Fr(1946, 1000)
    hb = m * m * need + Binf
    assert SE.ln_bounds(X.XMAX)[0] > hb, 'height search too long: refused'
    X0 = 1
    while SE.ln_bounds(X0 + 1)[0] <= hb: X0 += 1
    found, tors = [], []
    e = 1
    while e * e <= X0:
        e2 = e * e; e4 = e2 * e2; e6 = e4 * e2
        for a in range(-X0, X0 + 1):
            if gcd(a, e) != 1: continue
            # x = a/e^2 is the x of a point iff 4x^3 + b2 x^2 + 2 b4 x + b6 is a square; times e^6:
            v = 4 * a ** 3 + Em.b2 * a * a * e2 + 2 * Em.b4 * a * e4 + Em.b6 * e6
            if v == 0: tors.append(Fr(a, e2))
            elif v > 0 and isqrt(v) ** 2 == v: found.append(Fr(a, e2))
        e += 1
    assert not found, ('points of small height found', found)
    return dict(m=m, Binf=Binf, hb=hb, X0=X0, tors=tors, bad=bad)

def torsion2(Ed):
    """C(Q)_tors = Z/2, exactly: (i) the gcd of #C(F_q) over good q < 400 is a power of 2, so there is no odd torsion;
    (ii) the cubic has one rational root e (the quotient quadratic has a non-square discriminant), so C(Q)[2] = Z/2 and
    the torsion is cyclic of 2-power order; (iii) T = (e, 0) is not 2P: in y^2 = X(X^2 + a X + b), X = x - e,
    X(2P) = (X^2 - b)^2 / (4 y^2) = 0 needs X(P)^2 = b, and neither X = sqrt(b) nor -sqrt(b) gives a square f."""
    g = 0; used = []
    for q in range(3, 400):
        if RI.is_prime_td(q) and Ed.D % q: g = gcd(g, X.count_mod(Ed, q)); used.append(q)
    while g % 2 == 0: g //= 2
    assert g == 1, ('odd torsion not excluded', g)
    rts = [int(v) for v in pari.nfroots(None, pari('x^3+(%d)*x^2+(%d)*x+(%d)' % (Ed.a2, Ed.a4, Ed.a6)))]   # hint
    e = rts[0]; assert e ** 3 + Ed.a2 * e * e + Ed.a4 * e + Ed.a6 == 0
    a, b = 3 * e + Ed.a2, 3 * e * e + 2 * Ed.a2 * e + Ed.a4
    dq = a * a - 4 * b
    assert dq < 0 or isqrt(dq) ** 2 != dq, 'three rational 2-torsion points'
    if b >= 0 and isqrt(b) ** 2 == b:
        d = isqrt(b)
        for Xp in (d, -d):
            v = Xp * (Xp * Xp + a * Xp + b)
            assert not (v >= 0 and isqrt(v) ** 2 == v), 'T is 2P: a point of order 4'
    return (Fr(e), Fr(0)), used

def run(mn, G_override=None, say=print):
    m, n = mn
    ainv = dm_g4k.g4k_curve(m, n)[0]; _, a2, _, a4, a6 = ainv
    S = m * m + n * n; g4 = [2 * S, 2 * S]; K = [S * S, 2 * S * S - 16 * m * m * n * n, S * S]
    c3 = g4[1] * K[2]
    Ed = X.Ell(a2, a4, a6)
    # hints from PARI, checked exactly
    E = pari.ellinit(ainv)
    pts = pari.ellrank(E, 2)[3]
    if len(pts) == 0: pts = pari.ellrank(E, 4)[3]
    G = pari.ellsaturation(E, pts, 1000)[0]
    G = (fr(G[0]), fr(G[1])) if G_override is None else G_override
    assert Ed.on(G), 'G is not on the curve'
    T, used = torsion2(Ed)
    lo, hi, B, U = X.enclose(Ed, G, 7)
    assert lo > 0, 'hhat(G) not shown positive'
    with contextlib.redirect_stdout(io.StringIO()):
        X.saturation(Ed, G, T, LMAX)
    nxt = 211
    need = hi / (nxt - 1) ** 2
    mm_ = pari('(E) -> my(v, M = ellminimalmodel(E, &v)); [M[1..5], v]')(E)
    Emin_hint, uv = mm_[0], mm_[1]
    u, r, s_, t = (fr(v) for v in uv)
    amin = change([0, a2, 0, a4, a6], u, r, s_, t)
    assert all(c.denominator == 1 for c in amin), 'the minimal-model hint is not integral'
    assert [int(c) for c in amin] == [int(v) for v in Emin_hint[:5]]
    Em = GenEll([int(c) for c in amin])
    L = lower_height(Em, need)
    say('(%d,%d) C: y^2 = x^3 + c2 x^2 + c1 c3 x + c0 c3^2 (x = c3 s)' % (m, n))
    say('  torsion Z/2 (no odd part by point counts at %d good primes; one rational root; T not in 2C(Q)); hhat(G) in [%.6f, %.6f] (exact rational enclosure)'
        % (len(used), float(lo), float(hi)))
    say('  saturation of G at every prime <= %d (exact)' % LMAX)
    say('  minimal model [%s] (from the hint, checked: integral, minimal at each of %d bad primes)'
        % (', '.join(str(int(c)) for c in amin), len(L['bad'])))
    say('  m = %d, height search to %d: no point of infinite order of small height, 2-torsion x %s'
        % (L['m'], L['X0'], [str(v) for v in L['tors']]))
    say('  so every point of infinite order has hhat >= hi/210^2; IF rank C(Q) = 1, G generates C(Q) modulo torsion')
    # (5) the Demjanenko-Manin bound, exact
    C = 2 * SE.ln_bounds(abs(c3))[1]
    lhs = (B + U + C) / lo
    amax = int((lhs + 1) / 2)                   # 2a - 1 <= lhs
    cand = {}
    P = None
    for a in range(0, amax + 1):
        for Q in (P, Ed.add(P, T)):
            if Q is not None: cand.setdefault(Q[0], 'a = %d' % a)
        P = Ed.add(P, G)
    e = T[0]; beta = 3 * e * e + 2 * a2 * e + a4
    roots = [Fr(c3), Fr(-c3)]
    qa, qb, qc = e, beta - e * e - c3 * c3, c3 * c3 * e
    if qa == 0: roots.append(Fr(-qc, qb)) if qb else None
    else:
        dsc = qb * qb - 4 * qa * qc
        if dsc >= 0 and isqrt(dsc) ** 2 == dsc:
            roots += [Fr(-qb + isqrt(dsc), 2 * qa), Fr(-qb - isqrt(dsc), 2 * qa)]
    for x0 in roots: cand.setdefault(x0, 'P\' = +-P + T')
    verdicts = {}
    for X_ in cand:
        if X_ == e: continue
        for tt in to_t(X_, c3): verdicts[tt] = classify_t(m, n, tt)
    verdicts[Fr(0)] = classify_t(m, n, Fr(0))     # t = 0 (x = 0, P' = O) is tested directly
    hits = {str(k): v for k, v in verdicts.items() if v != 'no'}
    status = 'ALERT' if 'SOLUTION' in verdicts.values() else 'closed'
    say('  Demjanenko-Manin: B + U + 2 log|c3| <= %.3f, so a <= %d; %d candidate x, %d values of t; not "no": %s; %s'
        % (float(B + U + C), amax, len(cand), len(verdicts), hits, status))
    say('  conclusion: IF rank C(Q) = 1, every solution t of (%d,%d) is degenerate' % (m, n) if status == 'closed'
        else '  conclusion: NOT closed')
    return dict(status=status, amax=amax, hits=hits, G=G, Ed=Ed, T=T)

def negative(say):
    ok = []
    def refused(nm, fn):
        try:
            with contextlib.redirect_stdout(io.StringIO()): fn()
            say('  NOT REFUSED: ' + nm); ok.append(False)
        except AssertionError as ex:
            say('  refused: %s (%s)' % (nm, str(ex).split(',')[0].strip("('"))); ok.append(True)
    say('controls on (7735,4968):')
    base = run((7735, 4968), say=lambda s: None)
    G, Ed, T = base['G'], base['Ed'], base['T']
    G3 = Ed.mul(3, G)
    refused('3G in place of G (not saturated at 3)', lambda: run((7735, 4968), G_override=G3, say=lambda s: None))
    refused('the height search on y^2 = x^3 - 2, which has the point (3, 5) of small height (the search must see it)',
            lambda: lower_height(GenEll([0, 0, 0, 0, -2]), Fr(1, 1000)))
    refused('a point that is not on the curve', lambda: run((7735, 4968), G_override=(G[0] + 1, G[1]), say=lambda s: None))
    global classify_t
    keep = classify_t
    classify_t = lambda mm, nn, t: 'SOLUTION' if t == 1 else keep(mm, nn, t)
    try:
        r = run((7735, 4968), say=lambda s: None)
        good = r['status'] == 'ALERT'; ok.append(good)
        say('  classifier forged to call t = 1 a solution: status %s (must be ALERT)' % r['status'])
    finally: classify_t = keep
    say('  unchanged run: status %s, t not "no": %s (t = 1, -1 must be found)' % (base['status'], sorted(base['hits'])))
    ok.append(base['status'] == 'closed' and {'1', '-1'} <= set(base['hits']))
    say('negative checks: %s' % ('all as required (%d of %d)' % (sum(ok), len(ok)) if all(ok) else 'FAILED'))
    return all(ok)

if __name__ == '__main__':
    OUT = []
    def say(s): OUT.append(s); print(s, flush=True)
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if '--negative' in sys.argv:
        negative(say); name = 'g4k_exact_negative.txt'
    else:
        for a in (args or RATIOS): run(tuple(int(v) for v in a.split(',')), say=say)
        name = 'g4k_exact.txt'
    txt = '\n'.join(OUT) + '\n'
    if '--write' in sys.argv: open(os.path.join(HERE, name), 'w').write(txt)
    elif not args and os.path.exists(os.path.join(HERE, name)):
        assert txt == open(os.path.join(HERE, name)).read(), 'output differs from ' + name
        print('# matches ' + name)
