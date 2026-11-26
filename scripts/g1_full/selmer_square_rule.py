#!/usr/bin/env python3
"""The Selmer groups of selmer_square.py predicted from Legendre symbols alone, and the lemmas behind the prediction.
Experiment (item C1 of AGENT_TASKS.md): nothing counted, nothing in the census.

    python3 selmer_square_rule.py [--max M] [--write]   identities (sympy), congruence lemmas (finite checks), then the
                                                         rule against every curve of selmer_square_<M>.json
    python3 selmer_square_rule.py --negative            forged records / forged rules that must be refused
    python3 selmer_square_rule.py --max 5000 --tsv-to PATH   also the .tsv, to PATH (for CI; nothing is written here)

Written by an AI model (Claude); not read by a person.  The lemmas and their proofs are in SELMER_SQUARE.md; this script
checks each computational step: the polynomial identities in general m, n (sympy), the congruences (finite checks over
residues), the hypotheses on every curve (prime congruences, coprimality, non-squares), and then the predicted local
images and Selmer groups against the exact descent of selmer_square.py (rank_indep.solvable, no prediction used there).

THE RULE (E = E+ at (m, n), coprime, m^2 + n^2 = k^2, mn != 0; E- is n -> -n).  a = 32 m^3 n^3 - 4 (m^2 - n^2)^2 S.
  P_b  = odd primes of k H(m, n) H(m, -n),  P_b' = odd primes of m n (m - n)(m + n)  (disjoint).
  p in P_b, j = v_p(b):  im_p(E) = all four classes unless (j even and (a/p) = -1), then the two unit classes;
                         im_p(E') = {1} unless (j even and (a/p) = -1), then the two unit classes.
  p in P_b' (split iff (-2a/p) = 1):  im_p(E) = {1} if split, the two unit classes if not;
                                       im_p(E') = all four if split, the two unit classes if not.
  p = 2:  im_2(E) = {1, 2},  im_2(E') = {1, -1, 2, -2}  (mod squares).      R:  im(E) = {1},  im(E') = {1, -1}.
  (-2a/p) = (2/p) for p | mn,  = 1 for p | m + n,  = (-1/p) for p | m - n.
So Sel(E) = {d > 0 squarefree, d | 2 prod P_b : no p in P_b with (j even, (a/p) = -1) divides d, odd part of d = 1 mod 8,
(d/p) = 1 for every split p in P_b'} and Sel(E') = {d squarefree, d | 2 prod P_b' : no non-split p divides d, odd part
of d = +-1 mod 8, (d/p) = 1 for every p in P_b except those with (j even, (a/p) = -1)}.
"""
import os, sys, json
from math import gcd, isqrt
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import rank_indep as RI

WRITE = '--write' in sys.argv
out = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); out.append(s)

def H(m, n): return m**4 - 2*m**3*n + 2*m*m*n*n + 2*m*n**3 + n**4
def a_of(m, n): return 32*m**3*n**3 - 4*(m*m - n*n)**2*(m*m + n*n)

# ------------------------------------------------------------------ 1. identities in general m, n
def identities():
    import sympy as sp
    m, n, x = sp.symbols('m n x')
    S = m*m + n*n; u = m*m - n*n
    Hs = H(m, n); Hm = H(m, -n)
    a = -4*(m**6 - m**4*n**2 - 8*m**3*n**3 - m**2*n**4 + n**6)        # selmer_square.curve_data (checked against model)
    b = 4*S**2*Hs*Hm
    bp = -64*m**2*n**2*(m - n)**2*(m + n)**6
    kap = 4*m*n*(m - n)*(m + n)**3
    Z = lambda e: sp.expand(e) == 0
    chk = [
        ('I1  a = 32 m^3 n^3 - 4 (m^2 - n^2)^2 S', Z(a - a_of(m, n))),
        ("I2  b' = a^2 - 4b = -64 m^2 n^2 (m - n)^2 (m + n)^6 = -(2 kappa)^2, kappa = 4mn(m - n)(m + n)^3", Z(a*a - 4*b - bp) and Z(bp + 4*kap**2)),
        ('I3  x^2 + a x + b = (x + a/2)^2 + kappa^2', Z(x*x + a*x + b - ((x + a/2)**2 + kap**2))),
        ('I4  a - 32 m^3 n^3 = -4 (m^2 - n^2)^2 S  (so a = 32 m^3 n^3 mod S)', Z(a - 32*m**3*n**3 + 4*u*u*S)),
        ('I5  (m + n)^2 - 2mn = S  (so 2mn is a square mod p | S)', Z((m + n)**2 - 2*m*n - S)),
        ('I6  a at m = 0: -4 n^6;  at n = 0: -4 m^6;  at m = n: 32 n^6;  at m = -n: -32 n^6',
         Z(a.subs(m, 0) + 4*n**6) and Z(a.subs(n, 0) + 4*m**6) and Z(a.subs(m, n) - 32*n**6) and Z(a.subs(m, -n) + 32*n**6)),
        ('I7  H(m, n) - H(m, -n) = -4mn(m^2 - n^2);  H(n, n) = H(-n, n) = 4 n^4;  H(m, n) = n^4 at m = 0, m^4 at n = 0',
         Z(Hs - Hm + 4*m*n*u) and Z(Hs.subs(m, n) - 4*n**4) and Z(Hs.subs(m, -n) - 4*n**4) and Z(Hs.subs(m, 0) - n**4) and Z(Hs.subs(n, 0) - m**4)),
        ('I8  H(m, n) = S^2 - 2mn(m^2 - n^2)  (so H = -2mn(m^2 - n^2) mod S)', Z(Hs - (S*S - 2*m*n*u))),
        ('I9  f((x + 1)/(x - 1)) (x - 1)^6 = 8 f(x), f = (x^2 + 1) H(x, 1)  (UNIFORM_RANK.md, remark after Lemma D)',
         Z(sp.cancel(((x + 1)**2/(x - 1)**2 + 1)*H((x + 1)/(x - 1), 1)*(x - 1)**6) - 8*(x*x + 1)*H(x, 1))),
        ("I10 a^2 - b' = 4b  (so p | b gives a^2 = b' mod p, and p | b' gives a^2 = 4b mod p)", Z(a*a - bp - 4*b)),
    ]
    for name, ok in chk:
        assert ok, name
        say(name + ': identity (sympy expand)')
    return True

# ------------------------------------------------------------------ 2. congruences (finite checks)
def congruences():
    # C1: m, n of opposite parity, 4 | mn:  a/(-4) = 1 mod 8 and H(m, +-n) = 1 mod 8; a/(-4) mod 8 depends on m, n mod 8
    for m in range(16):
        for n in range(16):
            if (m + n) % 2 == 1 and (m * n) % 4 == 0:
                for nn in (n, -n):
                    assert (a_of(m, nn) // -4) % 8 == 1 and a_of(m, nn) % 4 == 0, (m, nn)
                    assert H(m, nn) % 8 == 1
    say('C1  m, n of opposite parity with 4 | mn (all residues mod 16): a = -4 a0 with a0 = 1 mod 8, H(m, +-n) = 1 mod 8')
    # C2: primitive Pythagorean (m, n): the even one is divisible by 4 (2rs with r, s of opposite parity)
    for r in range(2, 60):
        for s in range(1, r):
            if gcd(r, s) == 1 and (r + s) % 2 == 1:
                m, n = r*r - s*s, 2*r*s
                assert n % 4 == 0 and m % 2 == 1
                # m +- n = +-(x^2 - 2y^2) with gcd(x, y) = 1
                assert m + n == (r + s)**2 - 2*s*s and m - n == (r - s)**2 - 2*s*s
    say('C2  m = r^2 - s^2, n = 2rs (r > s coprime, opposite parity; r < 60): 4 | n, m + n = (r + s)^2 - 2s^2,',
        'm - n = (r - s)^2 - 2s^2 (identities, also for (m, n) swapped up to sign; so every odd prime of m +- n has (2/p) = 1)')
    return True

def lemma_controls(seed=20261125, N=150):
    """Lemmas L1, L2, L3 of SELMER_SQUARE.md compared with the exhaustive local solver on random curves outside the
    family (a control of the hand proofs, not part of them)."""
    import random
    rnd = random.Random(seed)
    def img(A, B, p):
        for boost in (0, 6, 12, 20):
            RI.DEPTH[2] = 14 + boost; RI.DEPTH_ODD = 8 + boost
            try: poss, cert = RI.local_image(A, B, p)
            finally: RI.DEPTH[2] = 14; RI.DEPTH_ODD = 8
            if poss == cert: return poss
        raise AssertionError('undecided')
    n1 = n2 = n3 = 0
    for _ in range(N):
        p = rnd.choice([3, 5, 7, 11, 13, 17]); j = rnd.randrange(1, 6)
        # L1: p odd, p not dividing A B, v_p(A^2 - 4B) = j
        A1 = rnd.randrange(1, 10**4)
        D1 = rnd.randrange(1, 10**3)
        if A1 % p == 0 or D1 % p == 0: continue
        A, B = 2 * A1, A1 * A1 - p**j * D1
        if B == 0 or B % p == 0: continue
        pred = UNITS if (j % 2 == 0 and RI.legendre(-2 * A, p) == -1) else {(0, 1)}
        assert img(A, B, p) == pred, ('L1', A, B, p)
        n1 += 1
    for _ in range(N):
        p = rnd.choice([3, 5, 7, 11, 13, 17]); j = rnd.randrange(1, 6)
        # L2: p odd, p not dividing A (A^2 - 4B), v_p(B) = j; p >= 5 or (A/p) = -1; p >= 7 if j is odd
        A = rnd.randrange(1, 10**4); B1 = rnd.choice([-1, 1]) * rnd.randrange(1, 10**3)
        if A % p == 0 or B1 % p == 0: continue
        if not (p >= 5 or RI.legendre(A, p) == -1) or (j % 2 == 1 and p < 7): continue
        B = p**j * B1
        if (A * A - 4 * B) == 0: continue
        pred = ALL4 if (j % 2 == 1 or RI.legendre(A, p) == 1) else UNITS
        assert img(A, B, p) == pred, ('L2', A, B, p)
        n2 += 1
    for _ in range(N // 3):
        # L3: A = -4 a0, a0 = 1 mod 8, B = 4 a0^2 + kappa^2 with v_2(kappa) >= 4; and the isogenous (8 a0, -4 kappa^2)
        a0 = 8 * rnd.randrange(-500, 500) + 1
        kap = 2**rnd.randrange(4, 9) * (2 * rnd.randrange(0, 500) + 1)
        A, B = -4 * a0, 4 * a0 * a0 + kap * kap
        assert img(A, B, 2) == {(0, 1), (1, 1)}, ('L3 E', a0, kap)
        assert img(-2 * A, A * A - 4 * B, 2) == {(0, 1), (0, 7), (1, 1), (1, 7)}, ("L3 E'", a0, kap)
        n3 += 1
    say('C3  lemmas L1 (%d curves), L2 (%d), L3 (%d pairs) agree with the exhaustive local solver on random curves'
        ' outside the family (control of the hand proofs)' % (n1, n2, n3))
    return True

# ------------------------------------------------------------------ 3. the rule
UNITS = {(0, 1), (0, -1)}
ALL4 = {(0, 1), (0, -1), (1, 1), (1, -1)}
def predicted_local(m, nn, F):
    """Predicted local images {p: (im_p(E), im_p(E'))} at the bad primes, from Legendre symbols only."""
    a = a_of(m, nn)
    Pb, Pbp = {}, {}
    for key in ('k', 'H', 'Hm'):
        for p, e in F[key].items():
            if p != 2: Pb[p] = Pb.get(p, 0) + (4*e if key == 'k' else e)
    for key in ('m', 'n', 'mmn', 'mpn'):
        for p, e in F[key].items():
            if p != 2: Pbp[p] = Pbp.get(p, 0) + 1
    assert not set(Pb) & set(Pbp)
    loc = {2: ({(0, 1), (1, 1)}, {(0, 1), (0, 7), (1, 1), (1, 7)})}
    info = {}
    for p, j in Pb.items():
        A = RI.legendre(a, p); assert A != 0
        restricted = (j % 2 == 0 and A == -1)
        loc[p] = (UNITS if restricted else ALL4, UNITS if restricted else {(0, 1)})
        info[p] = ('b', j, A, restricted)
    for p in Pbp:
        sp_ = RI.legendre(-2*a, p); assert sp_ != 0
        loc[p] = ({(0, 1)} if sp_ == 1 else UNITS, ALL4 if sp_ == 1 else UNITS)
        info[p] = ("b'", None, sp_, sp_ == -1)
    return loc, info

def selmer_from_local(gens, loc, positive_only):
    """Classes d in <gens> (squarefree, -1 allowed) whose class at every p of loc is in the image."""
    out = []
    for mask in range(1 << len(gens)):
        d = 1
        for i, g in enumerate(gens):
            if mask >> i & 1: d *= g
        if positive_only and d < 0: continue
        if all(RI.qp_class(d, p) in loc[p] for p in loc): out.append(d)
    return sorted(out, key=abs)

def check_hypotheses(m, nn, F):
    """The facts the lemmas use, on this curve."""
    n = abs(nn)
    S = m*m + n*n; k = isqrt(S); assert k*k == S and gcd(m, n) == 1 and m != n and n != 0
    assert (m*n) % 4 == 0 and (m + n) % 2 == 1
    for p in F['k']: assert p % 4 == 1
    for key in ('H', 'Hm'):
        for p in F[key]: assert p % 12 == 1, (key, p)
    for key in ('mmn', 'mpn'):
        for p in F[key]: assert p % 8 in (1, 7), (key, p)
    a = a_of(m, nn); assert (a // -4) % 8 == 1
    for p in F['k']: assert RI.legendre(a, p) == 1                       # a = 32 m^3 n^3 = (4mn)^2 * 2mn, 2mn = (m+n)^2 mod p
    for p in F['m']:
        if p != 2: assert RI.legendre(-2*a, p) == RI.legendre(2, p)
    for p in F['n']:
        if p != 2: assert RI.legendre(-2*a, p) == RI.legendre(2, p)
    for p in F['mpn']: assert RI.legendre(-2*a, p) == 1
    for p in F['mmn']: assert RI.legendre(-2*a, p) == RI.legendre(-1, p)
    for key in ('H', 'Hm'):
        v = H(m, nn) if key == 'H' else H(m, -nn)
        assert isqrt(v)**2 != v                                         # Lemma D via I9: S H(m, +-n) is not a square
    return True

def independent_solver(data, seed=11, N=60, PMAX=400):
    """The local images of the data recomputed with the separate exhaustive search of sigma_generic.py
    (locally_soluble over Q_p, code shared with nothing above), at every bad prime p < PMAX of N random curves."""
    import random
    import sigma_generic as SG
    R = SG.QO(None)
    rnd = random.Random(seed)
    sample = rnd.sample(data, N)
    nloc = 0
    for r in sample:
        m, nn = r['m'], r['n']
        a = a_of(m, nn); b = 4 * (m*m + nn*nn)**2 * H(m, nn) * H(m, -nn)
        ap, bp = -2 * a, a * a - 4 * b
        for p in r['bad']:
            if p >= PMAX: continue
            L = SG.Local(R, 'rat', p)
            for (A, B, key) in ((a, b, 'locS'), (ap, bp, 'locSp')):
                img = set()
                for delta in RI.class_reps(p):
                    # N^2 = delta M^4 + A M^2 e^2 + (B/delta) e^4, times delta^2: delta^3 M^4 + A delta^2 M^2 e^2 + B delta e^4
                    F = [(delta**3, 0), (0, 0), (A * delta * delta, 0), (0, 0), (B * delta, 0)]
                    if SG.locally_soluble(L, F, cap=60): img.add(RI.qp_class(delta, p))
                assert img == set(map(tuple, r[key][str(p)])), ('independent solver disagrees', m, nn, p, key)
                nloc += 1
    say('C4  the local images of %d random curves at every bad prime p < %d (%d local images) recomputed by the'
        ' separate exhaustive search of sigma_generic.py: all agree' % (N, PMAX, nloc))

def against_data(M):
    path = os.path.join(here, 'selmer_square_%d.json' % M)
    data = json.load(open(path))
    say('rule against the exact descent: %s, %d curves' % (os.path.basename(path), len(data)))
    independent_solver(data, N=min(200, len(data)), PMAX=3000)
    rows = []
    bad = {'fit': 0, 'test': 0}; tot = {'fit': 0, 'test': 0}
    from collections import Counter
    dist = Counter(); s2 = []
    for r in data:
        m, nn = r['m'], r['n']
        F = {k: {int(p): e for p, e in v} for k, v in r['F'].items()}
        check_hypotheses(m, nn, F)
        loc, info = predicted_local(m, nn, F)
        half = 'fit' if m <= M // 2 else 'test'
        tot[half] += 1
        locS = {int(p): set(map(tuple, v)) for p, v in r['locS'].items()}
        locSp = {int(p): set(map(tuple, v)) for p, v in r['locSp'].items()}
        ok = set(loc) == set(locS) and all(loc[p][0] == locS[p] and loc[p][1] == locSp[p] for p in loc)
        b_primes = sorted(set([2]) | {p for key in ('k', 'H', 'Hm') for p in F[key]})
        bp_primes = sorted(set([2]) | {p for key in ('m', 'n', 'mmn', 'mpn') for p in F[key]})
        SelE = selmer_from_local([-1] + b_primes, {p: loc[p][0] for p in loc}, True)
        SelEp = selmer_from_local([-1] + bp_primes, {p: loc[p][1] for p in loc}, False)
        ok = ok and set(SelE) == set(r['S']) and set(SelEp) == set(r['Sp'])
        s = (len(SelE).bit_length() - 1) + (len(SelEp).bit_length() - 1) - 2
        ok = ok and s == r['bound']
        # Theorem S1 on this curve, from the exact descent alone: 2, H, H- in Sel(E), independent; -1 in Sel(E')
        sqf = lambda x: RI.sqfree_part(x)
        hE = {sqf(H(m, nn)), sqf(H(m, -nn))}
        assert 2 in r['S'] and hE <= set(r['S']) and -1 in r['Sp'] and len({1, 2} | hE) == 4 and 1 not in hE
        assert r['bound'] >= 2
        if not ok: bad[half] += 1; say('  EXCEPTION', m, nn)
        dist[s] += 1
        if s == 2: s2.append((m, nn))
        # the Redei-type description: generators left after the restrictions, and the conditions
        genE = [2] + [p for p in b_primes if p != 2 and not info[p][3]]
        genEp = [-1, 2] + [p for p in bp_primes if p != 2 and not info[p][3]]
        rows.append((m, nn, s, len(SelE).bit_length() - 1, len(genE), len(SelEp).bit_length() - 1, len(genEp),
                     ' '.join('%d:%s%s' % (p, info[p][0], ('^%d' % info[p][1]) if info[p][1] else '') +
                              ('[(a/p)=%+d]' % info[p][2] if info[p][0] == 'b' else '[(-2a/p)=%+d]' % info[p][2])
                              for p in sorted(info))))
    for h in ('fit', 'test'):
        say('  %s half (m %s %d): %d curves, %d exceptions (local images at every bad prime, Sel(E), Sel(E\'), bound)'
            % (h, '<=' if h == 'fit' else '>', M // 2, tot[h], bad[h]))
    say('  Selmer bound from the rule: ' + ', '.join('%d: %d' % (k, dist[k]) for k in sorted(dist)) + '  (never 1)')
    say('  bound 2 (the least possible): %d curves: %s' % (len(s2), ' '.join('(%d,%d)' % x for x in s2)))
    assert bad['fit'] == 0 and bad['test'] == 0
    assert 1 not in dist and min(dist) >= 2
    return rows

def negative():
    """A forged record, a forged local image, a forged rule must be refused."""
    import selmer_square as SS
    res = []
    def run(name, fn):
        try:
            fn(); say('  NOT REFUSED: ' + name); res.append(False)
        except AssertionError as ex:
            say('  refused: %s (%s)' % (name, str(ex)[:70])); res.append(True)
    m, n = 12, 5
    r = SS.descend(m, n, 2)
    # (1) a forged Selmer element: a class d that some local test excludes, claimed to be in Sel(E)
    def f1():
        forged = sorted(set(r['S']) | {13})
        SS.verify_record(m, n, 2, forged, r['Sp'])
    run('forged element 13 added to Sel(E) of E+ (12,5)', f1)
    def f1b():
        where = SS.obstruction(m, n, 2, 13, 'E')
        say('    (C_13 of E+ (12,5) is certified insoluble at %s)' % where)
        assert where is None, 'C_13 has no point at %s' % where
    run('the class 13 claimed everywhere locally soluble on E+ (12,5)', f1b)
    def f3b():
        where = SS.obstruction(m, n, 2, 3, "E'")
        say("    (C_3 of E' of (12,5) is certified insoluble at %s)" % where)
        assert where is None, 'C_3 has no point at %s' % where
    run("the class 3 claimed everywhere locally soluble on E' of (12,5)", f3b)
    def f2b():
        assert SS.obstruction(m, n, 2, 2, 'E') is None                  # 2 is in Sel(E): no obstruction (control)
        assert SS.obstruction(m, n, 2, 2, 'E') is not None, 'the class 2 has no local obstruction (as Theorem S1 says)'
    run('an obstruction claimed for the class 2 on E+ (12,5)', f2b)
    def f2():
        forged = [d for d in r['S'] if d not in (2, )]
        SS.verify_record(m, n, 2, forged, r['Sp'])
    run('class 2 removed from Sel(E) of E+ (12,5)', f2)
    def f3():
        forged = sorted(set(r['Sp']) | {3})
        SS.verify_record(m, n, 2, r['S'], forged)
    run("forged element 3 added to Sel(E') of E+ (12,5)", f3)
    # (2) a wrong local image: the rule with im_2(E) = {1, 2, 5, 10}
    def f4():
        F = r['F']; loc, info = predicted_local(m, n, F)
        loc[2] = ({(0, 1), (1, 1), (0, 5), (1, 5)}, loc[2][1])
        locS = {p: set(map(tuple, v)) for p, v in r['locS'].items()}
        assert all(loc[p][0] == locS[p] for p in loc), 'forged im_2(E) differs from the exact local image'
    run('a wrong 2-adic image im_2(E) = {1, 2, 5, 10}', f4)
    def f5():
        F = r['F']; loc, info = predicted_local(m, n, F)
        p = 17                                                          # 17 | m + n: split, im_17(E) = {1}
        loc[p] = (UNITS, loc[p][1])
        locS = {q: set(map(tuple, v)) for q, v in r['locS'].items()}
        assert all(loc[q][0] == locS[q] for q in loc), 'forged im_17(E) differs from the exact local image'
    run('a wrong local image im_17(E) = units at 17 | m + n', f5)
    # (3) a wrong local image inside the exact descent: solvable() made to answer True for every disc at 17
    def f6():
        saved = RI.solvable
        RI.solvable = lambda forms, p, depth=None: True if p == 17 else saved(forms, p, depth)
        try: SS.descend(m, n, 2)
        finally: RI.solvable = saved
    run('the local solver forced to accept every class at 17 (product formula must fail)', f6)
    # (4) a hypothesis of the lemmas violated: a non-square-type ratio fed to the rule
    def f7():
        check_hypotheses(5, 2, {'k': {}, 'H': {}, 'Hm': {}, 'm': {5: 1}, 'n': {2: 1}, 'mmn': {3: 1}, 'mpn': {7: 1}})
    run('a ratio of the wrong type (5,2) given to the rule', f7)
    return all(res), res

if __name__ == '__main__':
    if '--negative' in sys.argv:
        ok, res = negative()
        say('negative tests: %d of %d refused' % (sum(res), len(res)))
        if WRITE: open(os.path.join(here, 'selmer_square_rule_negative.txt'), 'w').write('\n'.join(out) + '\n')
        sys.exit(0 if ok else 1)
    M = 5000
    if '--max' in sys.argv: M = int(sys.argv[sys.argv.index('--max') + 1])
    say('selmer_square_rule.py: the lemmas of SELMER_SQUARE.md, checked; then the rule against selmer_square_%d.json' % M)
    identities(); congruences(); lemma_controls()
    rows = against_data(M)
    if '--tsv-to' in sys.argv:
        # for CI: the .tsv written to the path named here (under /tmp), not into the repository
        with open(sys.argv[sys.argv.index('--tsv-to') + 1], 'w') as f:
            f.write('m\tn_signed\tbound\tdim Sel(E)\tgens E\tdim Sel(E\')\tgens E\'\todd bad primes p:(b or b\')^v[symbol]\n')
            for row in rows: f.write('\t'.join(map(str, row)) + '\n')
    if WRITE:
        with open(os.path.join(here, 'selmer_square_rule_%d.txt' % M), 'w') as f: f.write('\n'.join(out) + '\n')
        with open(os.path.join(here, 'selmer_square_rule_%d.tsv' % M), 'w') as f:
            f.write('m\tn_signed\tbound\tdim Sel(E)\tgens E\tdim Sel(E\')\tgens E\'\todd bad primes p:(b or b\')^v[symbol]\n')
            for row in rows: f.write('\t'.join(map(str, row)) + '\n')
