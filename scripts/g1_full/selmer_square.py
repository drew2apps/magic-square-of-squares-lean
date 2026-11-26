#!/usr/bin/env python3
"""Exact 2-isogeny Selmer groups of E+ and E- at the square-type ratios (m^2 + n^2 a square).  Experiment (item C1 of
AGENT_TASKS.md): nothing here is counted, nothing goes into the census.

    python3 selmer_square.py [--max M] [--jobs J] [--write]     table for coprime m > n >= 1, m <= M (default 2500)
    python3 selmer_square.py --negative                           forged inputs that must be refused
    python3 selmer_square.py --max M --dump-tsv | --dump-json     the rows of the .tsv / the .json on stdout (for CI)

Written by an AI model (Claude); not read by a person.

THE CURVES.  E+ at (m, n) is sigma_exact.model(m, n, 2), E- is model(m, n, -2) = model(m, -n, 2) (checked on each
curve).  The root of the cubic is e = -2 S^3 (S = m^2 + n^2), the other two roots are not real (the quadratic factor
has discriminant -64 m^2 n^2 (m - n)^2 (m + n)^6 < 0), so (e, 0) is the only rational point of order 2.  With
x = X - e the curve is  E: y^2 = x(x^2 + a x + b)  and the 2-isogenous curve is  E': Y^2 = X(X^2 + a' X + b'), with
(checked exactly on every curve against the model)
    a  = -4(m^6 - m^4 n^2 - 8 m^3 n^3 - m^2 n^4 + n^6),      b  = 4 S^2 H(m, n) H(m, -n),
    a' = -2a,                                               b' = a^2 - 4b = -64 m^2 n^2 (m - n)^2 (m + n)^6,
H(m, n) = m^4 - 2m^3 n + 2m^2 n^2 + 2m n^3 + n^4 (UNIFORM_RANK.md).  For E- replace n by -n.

THE DESCENT (Silverman-Tate III.4-III.6, as in rank_indep.py, whose functions are used unchanged).  Sel(E) is the set
of squarefree d | b (any sign) such that C_d: N^2 = d M^4 + a M^2 e^2 + (b/d) e^4 has a point over R and over Q_p
for every p | 2 b b'; Sel(E') the same with (a', b').  (For p not dividing 2 b b', C_d has good reduction and a point:
rank_indep.py's docstring.)  The images of E(Q), E'(Q) lie in them, and 2^rank = |im alpha| |im alpha'| / 4, so
    rank E(Q) <= s := log2 |Sel(E)| + log2 |Sel(E')| - 2,
the "Selmer bound" of the table.  Known: alpha(T) = b, alpha(R) = 2 S H(m, n) (R of UNIFORM_RANK.md, 2R = P) and
alpha'(T') = b' = -1 mod squares, so |Sel(E)| >= 4 and |Sel(E')| >= 2 (Theorem 1 of UNIFORM_RANK.md: R has infinite
order, rank >= 1).  s = 1 exactly when Sel(E) = <b, 2SH> and Sel(E') = {1, -1}.

LOCAL IMAGES.  For each bad prime p and each class delta of Q_p*/Q_p*^2 the solubility of C_delta over Q_p is decided
by rank_indep.solvable (an exhaustive search over residue discs that certifies True or False, or says None).  The
table is accepted only if every local answer is certified (no None) and the local images satisfy the product formula
    |im_p(E)| |im_p(E')| = 4 for odd p,  8 for p = 2
(|E(Q_p)/phi^(E'(Q_p))| |E'(Q_p)/phi(E(Q_p))| = |E[phi](Q_p)| |E'[phi^](Q_p)| |2|_p^-1, a consequence of
Schaefer's formula |E'(Q_p)/phi E(Q_p)| = |E[phi](Q_p)| |phi'(0)|_p^-1 c_p(E')/c_p(E), since phi^ phi = 2), and each
local image is a subgroup containing the classes of the known points.  This formula is used as a check of the solver,
not in the bound.  The real place: real_ok (exact).

FACTORISATIONS.  b and b' are factored through their factors m, n, m - n, m + n, sqrt(S), H(m, n), H(m, -n) (each by
rank_indep.factor_int: found by PARI or sympy, checked by multiplying back and proving every factor prime).

OUTPUT (with --write): selmer_square_<M>.tsv, one line per curve; selmer_square_<M>.json, the local images and the
Selmer groups (read by selmer_square_rule.py); selmer_square_<M>.txt, the summary printed on stdout (progress and time
go to stderr, so stdout can be compared with the file).  M = 5000 takes about 26 minutes on two cores (--jobs 2).
--negative runs the forged inputs of selmer_square_rule.negative (with --write: selmer_square_negative.txt).
"""
import os, sys, json, time
from fractions import Fraction as Fr
from math import gcd, isqrt

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import rank_indep as RI
import sigma_exact as SE

WRITE = '--write' in sys.argv

# Factorisations are found by sympy (one Pari instance per call, as in rank_indep._find_factor, uses too much memory
# over two thousand curves); rank_indep.factor_int still checks each one (product, and every factor proved prime).
def _find_factor_sympy(N):
    import sympy
    return {int(p): int(k) for p, k in sympy.factorint(abs(int(N))).items()}
RI._find_factor = lambda N: (
    {int(p): int(k) for p, k in RI.HINTS['factor'][str(N)]} if str(N) in RI.HINTS['factor'] else _find_factor_sympy(N))

_vp_orig = RI.vp
def _vp_fast(x, p):
    """Same values as rank_indep.vp (checked on random inputs at import), faster for p = 2."""
    assert x != 0
    if p == 2: return ((x & -x).bit_length() - 1) if x > 0 else (((-x) & x).bit_length() - 1)
    return _vp_orig(x, p)
import random as _r
_rg = _r.Random(7)
for _ in range(2000):
    _x = _rg.choice([-1, 1]) * _rg.randrange(1, 1 << 200) << _rg.randrange(0, 40)
    for _p in (2, 3, 5):
        assert _vp_fast(_x, _p) == _vp_orig(_x, _p)
RI.vp = _vp_fast

def H(m, n): return m**4 - 2*m**3*n + 2*m*m*n*n + 2*m*n**3 + n**4

def square_ratios(M):
    out = []
    for m in range(2, M + 1):
        for n in range(1, m):
            if gcd(m, n) == 1:
                S = m*m + n*n
                if isqrt(S)**2 == S: out.append((m, n))
    return out

def fac_str(F):
    if not F: return '1'
    return '*'.join(('%d' % p) if k == 1 else ('%d^%d' % (p, k)) for p, k in sorted(F.items()))

def merge(*Fs):
    out = {}
    for F in Fs:
        for p, k in F.items(): out[p] = out.get(p, 0) + k
    return out

def curve_data(m, n, s):
    """The shifted model of E+ (s = 2) or E- (s = -2) at (m, n), with every closed form checked against model()."""
    M = SE.model(m, n, s)
    nn = n if s == 2 else -n
    M2 = SE.model(m, nn, 2)
    assert (M['a2'], M['a4'], M['a6']) == (M2['a2'], M2['a4'], M2['a6']), 'E- is not E+ at (m, -n)'
    a2, a4, a6 = M['a2'], M['a4'], M['a6']
    S = m*m + n*n
    e = -2 * S**3
    assert e**3 + a2*e*e + a4*e + a6 == 0, 'e is not a root'
    c1 = a2 + e; c0 = a4 + e*c1                         # cubic = (X - e)(X^2 + c1 X + c0)
    assert a6 == -e*c0, 'not a factorisation'
    assert c1*c1 - 4*c0 < 0, 'the quadratic factor has real roots'   # so e is the only rational root
    a = 3*e + a2; b = 3*e*e + 2*a2*e + a4
    assert a == -4*(m**6 - m**4*nn**2 - 8*m**3*nn**3 - m*m*nn**4 + nn**6)
    assert b == 4 * S*S * H(m, nn) * H(m, -nn)
    ap, bp = -2*a, a*a - 4*b
    assert bp == -64 * m*m * nn*nn * (m - nn)**2 * (m + nn)**6
    return dict(nn=nn, a=a, b=b, ap=ap, bp=bp, S=S)

def descend(m, n, s, depth_boost=0):
    D = curve_data(m, n, s)
    nn, a, b, ap, bp, S = D['nn'], D['a'], D['b'], D['ap'], D['bp'], D['S']
    k = isqrt(S); assert k*k == S
    parts = dict(m=m, n=abs(nn), mmn=abs(m - nn), mpn=abs(m + nn), k=k, H=H(m, nn), Hm=H(m, -nn))
    F = {key: RI.factor_int(v) for key, v in parts.items()}
    Fb = merge({2: 2}, {p: 4*e for p, e in F['k'].items()}, F['H'], F['Hm'])
    Fbp = merge({2: 6}, {p: 2*e for p, e in F['m'].items()}, {p: 2*e for p, e in F['n'].items()},
                {p: 2*e for p, e in F['mmn'].items()}, {p: 6*e for p, e in F['mpn'].items()})
    RI.HINTS['factor'][str(abs(b))] = sorted(Fb.items())
    RI.HINTS['factor'][str(abs(bp))] = sorted(Fbp.items())
    # depths of the disc search: rank_indep's defaults (14 at 2, 8 at odd p), raised with the valuations of b b'
    v2 = RI.vp(b * bp, 2); vodd = max([e for p, e in merge(Fb, Fbp).items() if p != 2] + [0])
    RI.DEPTH[2] = max(14, v2 + 4) + depth_boost; RI.DEPTH_ODD = max(8, vodd + 2) + depth_boost
    try:
        R = RI.first_descent(a, b)        # factor_int(b), factor_int(b') check the factorisations above
    finally:
        RI.DEPTH[2] = 14; RI.DEPTH_ODD = 8
    del RI.HINTS['factor'][str(abs(b))]; del RI.HINTS['factor'][str(abs(bp))]
    undecided = []
    for side, loc in (('E', R['locS']), ("E'", R['locSp'])):
        for p, (poss, cert) in loc.items():
            if poss != cert: undecided.append((side, p))
    # product formula and subgroup property of the local images (checks of the local solver)
    for p in R['bad']:
        nE, nEp = len(R['locS'][p][0]), len(R['locSp'][p][0])
        assert nE * nEp == (8 if p == 2 else 4) or undecided, ('product formula fails', m, nn, p, nE, nEp)
    # known classes
    sq = RI.sqfree_part
    known_E = [sq(b), sq(2 * S * H(m, nn)), sq(sq(b) * sq(2 * S * H(m, nn)))]
    known_Ep = [sq(bp)]
    assert sq(bp) == -1
    for d in known_E: assert d in R['S'], ('a known class of E is not in Sel(E)', d)
    for d in known_Ep: assert d in R['Sp'], ("a known class of E' is not in Sel(E')", d)
    # Selmer sets are groups
    for X in (R['S'], R['Sp']):
        Xs = set(X)
        for d1 in X:                       # d1, d2 squarefree: the squarefree part of d1 d2 is d1 d2 / gcd(d1, d2)^2
            for d2 in X: assert (d1 * d2) // gcd(d1, d2)**2 in Xs, 'Selmer set not closed under products'
    nS, nSp = len(R['S']), len(R['Sp'])
    assert nS & (nS - 1) == 0 and nSp & (nSp - 1) == 0
    bound = nS.bit_length() - 1 + nSp.bit_length() - 1 - 2
    assert bound >= 1
    return dict(m=m, n=nn, F=F, bad=R['bad'], S=sorted(R['S'], key=abs), Sp=sorted(R['Sp'], key=abs),
                bound=bound, undecided=undecided,
                locS={p: sorted(map(list, R['locS'][p][0])) for p in R['bad']},
                locSp={p: sorted(map(list, R['locSp'][p][0])) for p in R['bad']})

def verify_record(m, n, s, S_claim, Sp_claim, bound_claim=None):
    """A claimed pair of Selmer groups is accepted only if it equals the one recomputed here."""
    r = run_one((m, n, s))
    assert not r['undecided']
    assert set(S_claim) == set(r['S']), ('Sel(E) differs', sorted(set(S_claim) ^ set(r['S'])))
    assert set(Sp_claim) == set(r['Sp']), ("Sel(E') differs", sorted(set(Sp_claim) ^ set(r['Sp'])))
    if bound_claim is not None: assert bound_claim == r['bound'], 'bound differs'
    return True

def obstruction(m, n, s, d, side='E'):
    """The first place (R or a bad prime) where C_d is certified insoluble, for a squarefree d dividing b (side 'E')
    or b' (side "E'"); None if C_d is soluble everywhere (certified).  Used to show why a forged class is refused."""
    D = curve_data(m, n, s)
    A, B = (D['a'], D['b']) if side == 'E' else (D['ap'], D['bp'])
    assert B % d == 0
    if not RI.real_ok(d, A, B): return 'R'
    bad = sorted(set([2]) | set(RI.factor_int(D['b'])) | set(RI.factor_int(D['bp'])))
    for p in bad:
        r = RI.solvable([RI.quartic(d, A, B)], p)
        if r is False: return p
        assert r is True, ('undecided', p)
    return None

def run_one(args):
    m, n, s = args
    t0 = time.time()
    r = descend(m, n, s)
    boost = 0
    while r['undecided'] and boost < 8:
        boost += 4
        r = descend(m, n, s, depth_boost=boost)
    r['boost'] = boost; r['sec'] = time.time() - t0
    return r

def row(r):
    F = r['F']
    return '\t'.join([str(r['m']), str(r['n']), 'E+' if r['n'] > 0 else 'E-', str(r['bound']),
                      str(len(r['S'])), str(len(r['Sp'])),
                      ','.join(map(str, r['S'])), ','.join(map(str, r['Sp'])),
                      fac_str(F['m']), fac_str(F['n']), fac_str(F['mmn']), fac_str(F['mpn']), fac_str(F['k']),
                      fac_str(F['H']), fac_str(F['Hm']), str(len(r['undecided']))])

HEADER = '\t'.join(['m', 'n_signed', 'curve', 'selmer_bound', '|Sel(E)|', "|Sel(E')|", 'Sel(E)', "Sel(E')",
                    'fac m', 'fac |n|', 'fac |m-n|', 'fac |m+n|', 'fac sqrt(S)', 'fac H(m,n)', 'fac H(m,-n)', 'undecided'])

def main():
    M = 2500
    if '--max' in sys.argv: M = int(sys.argv[sys.argv.index('--max') + 1])
    J = 1
    if '--jobs' in sys.argv: J = int(sys.argv[sys.argv.index('--jobs') + 1])
    rats = square_ratios(M)
    tasks = [(m, n, s) for (m, n) in rats for s in (2, -2)]
    print('square-type coprime ratios with m <= %d: %d; curves: %d' % (M, len(rats), len(tasks)), file=sys.stderr, flush=True)
    t0 = time.time()
    if J > 1:
        from multiprocessing import Pool
        res = []
        with Pool(J) as pool:
            for i, r in enumerate(pool.imap(run_one, tasks, chunksize=2)):
                res.append(r)
                if (i + 1) % 200 == 0: print('  %d curves done (m = %d), %.0f s' % (i + 1, r['m'], time.time() - t0), file=sys.stderr, flush=True)
    else:
        res = [run_one(t) for t in tasks]
    und = [r for r in res if r['undecided']]
    assert not und, ('undecided local tests', [(r['m'], r['n'], r['undecided']) for r in und][:5])
    from collections import Counter
    cnt = Counter(r['bound'] for r in res)
    cntE = Counter((('E+' if r['n'] > 0 else 'E-'), r['bound']) for r in res)
    lines = ['selmer_square.py: exact 2-isogeny Selmer bounds, square-type coprime ratios m <= %d (%d ratios, %d curves)'
             % (M, len(rats), len(tasks)),
             'every local test certified; product formula |im_p(E)||im_p(E\')| = 4 (odd p), 8 (p = 2) holds at every bad prime',
             'Selmer bound distribution: ' + ', '.join('%d: %d' % (k, cnt[k]) for k in sorted(cnt)),
             'by curve: ' + ', '.join('%s %d: %d' % (c, k, v) for (c, k), v in sorted(cntE.items())),
             'both curves of a ratio with bound 1: %d of %d ratios'
             % (sum(1 for (m, n) in rats if all(r['bound'] == 1 for r in res if r['m'] == m and abs(r['n']) == n)), len(rats))]
    recs = [dict(m=r['m'], n=r['n'], bound=r['bound'], S=r['S'], Sp=r['Sp'], bad=r['bad'],
                 locS={str(p): v for p, v in r['locS'].items()},
                 locSp={str(p): v for p, v in r['locSp'].items()},
                 F={k: sorted(v.items()) for k, v in r['F'].items()}) for r in res]
    if '--dump-tsv' in sys.argv or '--dump-json' in sys.argv:
        # for CI: the rows of the .tsv, or the .json, on stdout (nothing is written); with m <= M they are the first
        # rows of selmer_square_5000.tsv / .json, which are ordered by m
        if '--dump-tsv' in sys.argv:
            print(HEADER)
            for r in res: print(row(r))
        else:
            json.dump(recs, sys.stdout); print()
        return
    for l in lines: print(l)
    print('time %.0f s' % (time.time() - t0), file=sys.stderr)
    if WRITE:
        with open(os.path.join(here, 'selmer_square_%d.tsv' % M), 'w') as f:
            f.write(HEADER + '\n')
            for r in res: f.write(row(r) + '\n')
        with open(os.path.join(here, 'selmer_square_%d.json' % M), 'w') as f:
            json.dump(recs, f)
        with open(os.path.join(here, 'selmer_square_%d.txt' % M), 'w') as f:
            f.write('\n'.join(lines) + '\n')

if __name__ == '__main__':
    if '--negative' in sys.argv:
        import selmer_square_rule as SR
        ok, res = SR.negative()
        SR.say('negative tests: %d of %d refused' % (sum(res), len(res)))
        if WRITE: open(os.path.join(here, 'selmer_square_negative.txt'), 'w').write('\n'.join(SR.out) + '\n')
        sys.exit(0 if ok else 1)
    main()
