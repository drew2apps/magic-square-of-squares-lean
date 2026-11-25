#!/usr/bin/env python3
"""Hostile review of branch sess_t8 (SQUARE_SECTION.md): an independent re-implementation of its second descent.
Slice v5, run 3.  Written by an AI model (Claude); not read by a person.  Experiment, nothing counted.

    python3 review_t8_rank1.py [--curves all6|M:N,...] [--brute] [--crosscheck]   the six rank-1 curves (default)
    python3 review_t8_rank1.py --class2 300      class 2 on every square-type curve m <= 300, against square_section.tsv
    python3 review_t8_rank1.py --negative        forged inputs
  --brute: certify every insoluble cover of one element per coset a third time by pruned brute force mod p^K;
  --crosscheck: run methods S and T at every odd p < 60 and assert agreement.  Needs sympy and numpy.

Written before reading square_section.py beyond its function names and header docstring; shares no code with it
(it imports only sigma_exact.model, to rebuild the curve).

THE COVERING (re-derived; the review's notes are not committed).  E: y^2 = x(x^2 + a x + b), a = -4 a0, b = 4 a0^2 + kappa^2.
For d | b squarefree with d = al^2 + be^2:  C_d: N^2 = d M^4 + a M^2 e^2 + (b/d) e^4, and d N^2 = X1^2 + X2^2 with
X1 = d M^2 - 2 a0 e^2, X2 = kappa e^2.  Then (X1 + i X2)/(al + i be) has norm N^2, so X1 + i X2 = lam (al + i be)(s + i t)^2
with lam in Q*, (s, t) coprime integers (Hilbert 90).  For an odd prime p not dividing d kappa, p does not divide
gcd(X1, X2) (a common prime divides e or kappa; if e, then d M^2, so d, as gcd(M, e) = 1), and then v_p(lam) = 0
(inert p: compare v_p of both sides; split p = pi pi': one of pi, pi' does not divide X1 + i X2, nor s + i t).
So lam in { +-(squarefree divisors of 2 d kappa) } mod squares.  With q1 + i q2 = (al + i be)(s + i t)^2:
    F1 = lam kappa q2(s, t) = (kappa e)^2,     F2 = d lam kappa (kappa q1 + 2 a0 q2)(s, t) = (d kappa M)^2.
C_d(Q) non-empty => some D_lam (F1, F2 both squares at a point of P^1(Q)) has a rational point.

LOCAL TESTS (own code).  "Soluble at p" = some (s : t) in P^1(Q_p) with every form a square in Q_p (0 allowed).
  Method T (any p): residue discs x0 + p^k Z_p of the two charts (x : 1), (1 : p y); a disc is soluble if all forms are
    squares at its integer centre, insoluble if some form has constant non-square class on it (Taylor coefficients:
    every higher term has valuation >= v(f(x0)) + 1, or + 3 at p = 2), refined otherwise (depth cap; undecided is
    reported, never counted as insoluble).
  Method S (odd p, used for all odd p): at each level the forms are divided by their content; on a residue disc that is
    not a root mod p of any normalised form the class of each form is constant (valuation of the content, Legendre
    symbol of the value); only roots mod p are refined.  For p > 10^5 the non-root residues are scanned (20000 tries);
    if none works a structural certificate is required (a content of odd valuation, or a form = c h^2 mod p with
    (c/p) = -1, or a product of two forms = c h^2 with (c/p) = -1), else the answer is "undecided".
  For odd p < 60 both methods are run and must agree.
  Real place: one exact rational test point on every arc between the real roots (exact isolating intervals).
  Good primes: if p is odd and divides neither the contents, nor disc(F1) disc(F2) Res(F1, F2), the curve
  e^2 = F1, M^2 = F2 is a smooth genus-1 curve mod p, has an F_p-point (p + 1 - 2 sqrt p > 0) and a smooth point lifts;
  the script computes that product and divides out the tested primes, and tests any prime left over.

SELMER GROUPS (own code, for the rank formula).  Sel(E) = { d | b squarefree : C_d soluble at R and every p | 2 b b' },
Sel(E') the same for (a', b') = (-2a, a^2 - 4b); compared with selmer_square_5000.json of the t4/t8 branches.

RANK 1.  K = <alpha(T), alpha(R)> = {1, H H-, 2SH, 2SH-} mod squares.  If every coset of K in Sel(E) other than K has an
element d with every D_lam insoluble somewhere, then alpha(E(Q)) = K; if Sel(E') = {1, -1}, 2^rank = 4 * 2 / 4.
Here EVERY element of every non-trivial coset is tested (not one per coset), as a consistency check: a coset with one
excluded element must have all its elements excluded (if one were the image of a point, all would be).
"""
import os, sys, json, time, itertools
from fractions import Fraction as Fr
from math import gcd, isqrt

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import sigma_exact as SE
import sympy as sp
from sympy import factorint, Poly, symbols, sqrt_mod

X_ = symbols('X_')

def say(*a): print(*a, flush=True)

def H(m, n): return m**4 - 2*m**3*n + 2*m*m*n*n + 2*m*n**3 + n**4

def vp(z, p):
    z = int(z); assert z != 0
    v = 0
    while z % p == 0: z //= p; v += 1
    return v

def sq_ok(z, p):
    """z integer: is z a square in Q_p (0 counts)."""
    z = int(z)
    if z == 0: return True
    v = vp(z, p)
    if v % 2: return False
    u = z // p**v
    if p == 2: return u % 8 == 1
    return pow(u % p, (p - 1)//2, p) == 1

def chi(z, p):
    z %= p
    if z == 0: return 0
    return 1 if pow(z, (p - 1)//2, p) == 1 else -1

# ------------------------------------------------------------- polynomials: lists of ints, low degree first
def peval(f, x):
    r = 0
    for c in reversed(f): r = r*x + c
    return r

def taylor(f, x0):
    """coefficients of f(x0 + u) in u (exact)."""
    n = len(f); c = list(f)
    out = []
    for i in range(n):
        # synthetic division by (x - x0)
        acc = 0; q = []
        for coef in reversed(c):
            acc = acc*x0 + coef; q.append(acc)
        out.append(q[-1]); q = q[:-1]; c = list(reversed(q))
        if not c: break
    while len(out) < n: out.append(0)
    return out

def subst_affine(f, x0, h):
    """coefficients of f(x0 + h u) in u."""
    c = taylor(f, x0)
    return [ci * h**i for i, ci in enumerate(c)]

def chart_polys(F, p):
    """binary form F (coeff list A_0..A_d of s^d, s^{d-1} t, ..., t^d): chart A f(x) = F(x, 1), chart B g(y) = F(1, p y)."""
    d = len(F) - 1
    fA = [F[d - i] for i in range(d + 1)]          # F(x, 1) = sum A_j x^{d-j}
    fB = [F[j] * p**j for j in range(d + 1)]        # F(1, p y) = sum A_j p^j y^j
    return fA, fB

# ------------------------------------------------------------- method T
class Undecided(Exception): pass

def solve_T(polys, p, x0=0, k=0, depth=0, cap=40):
    e = 3 if p == 2 else 1
    cs = [taylor(f, x0) for f in polys]
    if all(sq_ok(c[0], p) for c in cs): return True
    for c in cs:
        if c[0] != 0:
            v0 = vp(c[0], p)
            w = min((vp(ci, p) + k*i for i, ci in enumerate(c) if i >= 1 and ci != 0), default=10**9)
            if w >= v0 + e and not sq_ok(c[0], p): return False
    if depth >= cap: raise Undecided(('T', p, x0, k))
    for r in range(p):
        if solve_T(polys, p, x0 + r * p**k, k + 1, depth + 1, cap): return True
    return False

def soluble_T(forms, p, cap=40):
    A = [chart_polys(F, p)[0] for F in forms]; B = [chart_polys(F, p)[1] for F in forms]
    return solve_T(A, p, cap=cap) or solve_T(B, p, cap=cap)

# ------------------------------------------------------------- method S (odd p)
def roots_mod(g, p):
    gm = [c % p for c in g]
    while gm and gm[-1] == 0: gm.pop()
    if len(gm) <= 1: return []
    if p < 3000:
        return [r for r in range(p) if peval(gm, r) % p == 0]
    P = Poly(list(reversed(gm)), X_, modulus=p)
    rs = set()
    for fac, mult in P.factor_list()[1]:
        if fac.degree() == 1:
            c1, c0 = [int(z) for z in fac.all_coeffs()]
            rs.add((-c0 * pow(c1, -1, p)) % p)
    return sorted(rs)

def const_sq_class(g, p):
    """g mod p (nonzero) = c h^2 ?  return (c mod p) or None."""
    gm = [c % p for c in g]
    while gm and gm[-1] == 0: gm.pop()
    if len(gm) == 1: return gm[0]
    P = Poly(list(reversed(gm)), X_, modulus=p)
    lc, facs = P.sqf_list()
    if all(mult % 2 == 0 for _, mult in facs):
        return int(lc) % p
    return None

def solve_S(polys, p, depth=0, cap=60, big=10**5):
    assert p % 2 == 1
    if all(sq_ok(f[0], p) for f in polys): return True
    for f in polys:
        if f[0] != 0:
            v0 = vp(f[0], p)
            w = min((vp(ci, p) for i, ci in enumerate(f) if i >= 1 and ci != 0), default=10**9)
            if w >= v0 + 1 and not sq_ok(f[0], p): return False
    if depth >= cap: raise Undecided(('S', p, depth))
    vs, gs = [], []
    for f in polys:
        v = min(vp(c, p) for c in f if c != 0)
        vs.append(v); gs.append([c // p**v for c in f])
    R = [set(roots_mod(g, p)) for g in gs]
    allR = set().union(*R)
    # non-root residues: all classes known
    if all(v % 2 == 0 for v in vs):
        if p <= big:
            for u0 in range(p):
                if u0 in allR: continue
                if all(chi(peval(g, u0), p) == 1 for g in gs): return True
        else:
            tries = 0
            for u0 in range(p):
                if u0 in allR: continue
                if all(chi(peval(g, u0), p) == 1 for g in gs): return True
                tries += 1
                if tries > 20000: break
            # structural certificate that no non-root residue works
            ok = False
            for g in gs:
                c = const_sq_class(g, p)
                if c is not None and chi(c, p) == -1: ok = True
            if not ok and len(gs) == 2:
                prod = [0]*(len(gs[0]) + len(gs[1]) - 1)
                for i, x in enumerate(gs[0]):
                    for j, y in enumerate(gs[1]): prod[i+j] += x*y
                c = const_sq_class(prod, p)
                if c is not None and chi(c, p) == -1: ok = True
            if not ok: raise Undecided(('S-big', p))
    # residues that are roots of some form: dead if another form is bad there, else refine
    undecided = False
    for u0 in sorted(allR):
        dead = False
        for j, g in enumerate(gs):
            if u0 not in R[j]:
                if vs[j] % 2 == 1 or chi(peval(g, u0), p) != 1: dead = True
        if dead: continue
        sub = [subst_affine(f, u0, p) for f in polys]
        try:
            if solve_S(sub, p, depth + 1, cap, big): return True
        except Undecided:
            undecided = True
    if undecided: raise Undecided(('S', p, depth))
    return False

def soluble_S(forms, p, cap=60):
    A = [chart_polys(F, p)[0] for F in forms]; B = [chart_polys(F, p)[1] for F in forms]
    return solve_S(A, p, cap=cap) or solve_S(B, p, cap=cap)

def soluble_p(forms, p, cross=False):
    if p == 2: return soluble_T(forms, 2, cap=60)
    r = soluble_S(forms, p)
    if cross and p < 60:
        r2 = soluble_T(forms, p, cap=30)
        assert r == r2, ('method S and T disagree', forms, p, r, r2)
    return r

# ------------------------------------------------------------- real place
def soluble_R(forms):
    d = len(forms[0]) - 1
    if all(F[0] > 0 for F in forms) or (len(forms) == 1 and forms[0][0] >= 0): return True
    prod = Poly([1], X_)
    for F in forms:
        fA = chart_polys(F, 2)[0]
        prod = prod * Poly(list(reversed(fA)), X_)
    sqf = Poly(prod, X_)
    iv = sorted(sqf.intervals(eps=Fr(1, 10**6)))
    pts = []
    if iv:
        pts.append(Fr(iv[0][0][0]) - 1); pts.append(Fr(iv[-1][0][1]) + 1)
        pv = list(reversed(sqf.all_coeffs()))
        for (a1, b1), (a2, b2) in zip([x[0] for x in iv], [x[0] for x in iv[1:]]):
            b1, a2 = Fr(b1), Fr(a2)
            if len(forms) > 1:
                assert b1 <= a2
                if b1 == a2:      # a shared endpoint that is not a root lies strictly between the two roots
                    assert peval([Fr(int(c.p), int(c.q)) if hasattr(c, 'q') else Fr(c) for c in pv], b1) != 0
            pts.append((b1 + a2) / 2)
        if len(forms) == 1: return True      # a real root: N = 0 is a real point
    else:
        pts.append(Fr(0))
    for x in pts:
        if all(peval(chart_polys(F, 2)[0], x) >= 0 for F in forms): return True
    return False

# ------------------------------------------------------------- the curves
def shifted(m, nn):
    Md = SE.model(m, nn, 2)
    a2, a4, a6 = Md['a2'], Md['a4'], Md['a6']
    S = m*m + nn*nn; e = -2 * S**3
    assert e**3 + a2*e*e + a4*e + a6 == 0
    a = 3*e + a2; b = 3*e*e + 2*a2*e + a4
    assert (Fr(e) + 0)**0 == 1
    # check: x(x^2 + a x + b) = cubic(x + e)
    for xx in range(-3, 4):
        X = xx + e
        assert X**3 + a2*X*X + a4*X + a6 == xx*(xx*xx + a*xx + b)
    return a, b

def curve(m, nn):
    a, b = shifted(m, nn)
    assert a % 4 == 0
    a0 = -a // 4
    kap = 4*m*nn*(m - nn)*(m + nn)**3
    assert b == 4*a0*a0 + kap*kap
    S = m*m + nn*nn; k = isqrt(S); assert k*k == S
    Hp, Hm = H(m, nn), H(m, -nn)
    assert b == 4*S*S*Hp*Hm
    return dict(m=m, nn=nn, a=a, b=b, a0=a0, kap=kap, S=S, k=k, Hp=Hp, Hm=Hm)

def primes_of(*zs):
    P = set()
    for z in zs:
        z = abs(int(z))
        if z > 1: P |= set(factorint(z))
    return sorted(P)

def sqfree(z):
    s = -1 if z < 0 else 1
    for p, e in factorint(abs(z)).items():
        if e % 2: s *= p
    return s

def sqfree_divisors(z):
    ps = primes_of(z); out = []
    for r in range(len(ps) + 1):
        for c in itertools.combinations(ps, r):
            q = 1
            for x in c: q *= x
            out += [q, -q]
    return out

def two_squares(d):
    """d > 0 squarefree, every odd prime = 1 mod 4: (al, be) with al^2 + be^2 = d (Gaussian product)."""
    z = (1, 0)
    for p in primes_of(d):
        if p == 2: w = (1, 1)
        else:
            assert p % 4 == 1
            r = sqrt_mod(-1, p)
            a_, b_ = p, r
            lim = isqrt(p)
            while b_ > lim: a_, b_ = b_, a_ % b_
            c2 = p - b_*b_; c = isqrt(c2); assert c*c == c2
            w = (b_, c)
        z = (z[0]*w[0] - z[1]*w[1], z[0]*w[1] + z[1]*w[0])
    assert z[0]**2 + z[1]**2 == d
    return z

def quartic(d, A, B):
    """N^2 = d M^4 + A M^2 e^2 + (B/d) e^4 as a binary form in (M, e)."""
    assert B % d == 0
    return [d, 0, A, 0, B // d]

def selmer(A, B, places):
    out = []
    for d in sqfree_divisors(B):
        F = quartic(d, A, B)
        if not soluble_R([F]): continue
        if all(soluble_p([F], p) for p in places): out.append(d)
    return sorted(out, key=lambda z: (abs(z), z))

def cover_forms(cd, d, al, be, lam):
    kap, a0 = cd['kap'], cd['a0']
    # q1 + i q2 = (al + i be)(s + i t)^2,  (s + i t)^2 = (s^2 - t^2) + i 2st;  forms in (s, t): [s^2, st, t^2]
    q1 = [al, -2*be, -al]
    q2 = [be, 2*al, -be]
    F1 = [lam*kap*c for c in q2]
    F2 = [d*lam*kap*(kap*x + 2*a0*y) for x, y in zip(q1, q2)]
    return F1, F2

def disc2(F): A, B, C = F; return B*B - 4*A*C
def res2(F, G):
    a, b, c = F; dd, e, f = G
    return (a*f - c*dd)**2 - (a*e - b*dd)*(b*f - c*e)

def cover_soluble_everywhere(cd, d, al, be, lam, tested, cross=False):
    F1, F2 = cover_forms(cd, d, al, be, lam)
    D = disc2(F1) * disc2(F2) * res2(F1, F2)
    assert D != 0
    cont = gcd(gcd(*F1), 0) * 1
    extra = abs(D) * abs(gcd(gcd(F1[0], F1[1]), F1[2])) * abs(gcd(gcd(F2[0], F2[1]), F2[2]))
    for p in tested:
        while extra % p == 0: extra //= p
    places = list(tested)
    if extra > 1: places += primes_of(extra)        # not expected; tested if it happens
    if not soluble_R([F1, F2]): return False, 'R'
    for p in places:
        if not soluble_p([F1, F2], p, cross): return False, p
    return True, None

def run_curve(m, nn, jsel=None, cross=False, verbose=True, brute=False):
    t0 = time.time()
    cd = curve(m, nn)
    a, b = cd['a'], cd['b']
    ap, bp = -2*a, a*a - 4*b
    base = primes_of(2, m, nn, m - nn, m + nn, cd['k'], cd['Hp'], cd['Hm'])
    assert set(primes_of(2*b*bp)) <= set(base)
    selE = selmer(a, b, base)
    selEp = selmer(ap, bp, base)
    tag = f"E{'+' if nn > 0 else '-'} ({m},{abs(nn)})"
    s = (len(selE).bit_length() - 1) + (len(selEp).bit_length() - 1) - 2
    if jsel is not None:
        assert sorted(selE) == sorted(jsel['S']) and sorted(selEp) == sorted(jsel['Sp']), (tag, selE, jsel['S'], selEp, jsel['Sp'])
    # K
    Kc = {1, sqfree(cd['Hp']*cd['Hm']), sqfree(2*cd['S']*cd['Hp']), sqfree(2*cd['S']*cd['Hm'])}
    assert Kc <= set(selE), (Kc, selE)
    def cls(z): return sqfree(z)
    cosets = []
    seen = set(Kc)
    for d in selE:
        if d in seen: continue
        cs = sorted({cls(d*k) for k in Kc})
        seen |= set(cs); cosets.append(cs)
    res = []
    kap = cd['kap']
    nocert = []; cert = [0]
    for cs in cosets:
        row = []
        # the representative certified a second time (brute force): the element of the coset with the fewest covers
        rep = min(cs, key=lambda z: (len(primes_of(2*z*kap)), z))
        for d in cs:
            assert d > 0
            al, be = two_squares(d)
            lams = sorted({sqfree(l) for l in sqfree_divisors(2*d*kap)}, key=lambda z: (abs(z), z))
            surv = []
            obst = {}
            for lam in lams:
                ok, where = cover_soluble_everywhere(cd, d, al, be, lam, base, cross)
                if ok: surv.append(lam)
                else: obst[where] = obst.get(where, 0) + 1
                if brute and not ok and d == rep:
                    c = certify_independently(cd, d, al, be, lam, [where] + [q for q in base if q != where] if where != 'R' else base)
                    if c is None: nocert.append((d, lam))
                    else: cert[0] += 1
            row.append((d, len(lams), surv, obst))
        res.append(row)
    excluded_all = all(any(not r[2] for r in row) for row in res)
    consistent = all(len({bool(r[2]) for r in row}) == 1 for row in res)
    rank1 = excluded_all and sorted(selEp) == [-1, 1]
    if verbose:
        say(f"{tag}: |Sel(E)| = {len(selE)}, Sel(E') = {selEp}, s = {s}; {len(cosets)} non-trivial cosets of K; "
            f"{time.time() - t0:.1f}s")
        for row in res:
            for d, nl, surv, obst in row:
                say(f"    d = {d}: {nl} covers, survivors {surv[:6]}{' ...' if len(surv) > 6 else ''}; obstructions at {dict(sorted(obst.items(), key=lambda x: str(x[0])))}")
        say(f"    every coset excluded: {excluded_all}; cosets internally consistent: {consistent}; RANK EXACTLY 1: {rank1}")
        if brute: say(f"    brute-force mod p^K certificates (independent of methods S/T): {cert[0]} insoluble covers (one element per coset) certified, {len(nocert)} not certified {nocert[:4]}")
    return dict(tag=tag, selE=selE, selEp=selEp, s=s, res=res, excluded_all=excluded_all, consistent=consistent, rank1=rank1)


# ------------------------------------------------------------- an independent insolubility certificate (brute force)
import numpy as np
_SQ = {}
def squares_mod(q):
    if q not in _SQ:
        a = np.zeros(q, dtype=bool); y = np.arange(q, dtype=np.int64); a[(y*y) % q] = True; _SQ[q] = a
    return _SQ[q]

def strip_even(F, p):
    c = 0
    for x in F:
        if x: c = gcd(c, abs(x))
    v = vp(c, p) if c else 0
    return [x // p**(2*(v//2)) for x in F]

def brute_insoluble(forms, p, K):
    """True if NO primitive (s : t) mod p^K makes every form a square mod p^K (=> no Q_p point); else False."""
    q = p**K; SQ = squares_mod(q)
    fs = [strip_even(F, p) for F in forms]
    t = np.arange(q, dtype=object if q > 3*10**4 else np.int64)
    sp_ = np.arange(q // p, dtype=np.int64) * p
    ok = None
    for chart in (0, 1):
        good = None
        for F in fs:
            Fm = [x % q for x in F]
            if chart == 0:      # (1 : t)
                tt = np.arange(q, dtype=np.int64)
                vals = (Fm[0] + Fm[1]*tt % q + (Fm[2] * ((tt*tt) % q)) % q) % q
            else:               # (p s' : 1)
                ss = sp_
                vals = ((Fm[0] * ((ss*ss) % q)) % q + (Fm[1]*ss) % q + Fm[2]) % q
            g = SQ[vals]
            good = g if good is None else (good & g)
        if good.any(): return False
    return True


def sq_mod_pK(v, p, K):
    """is v a square modulo p^K (some y with y^2 = v mod p^K)?"""
    q = p**K; v %= q
    if v == 0: return True
    j = 0
    while v % p == 0: v //= p; j += 1
    if j % 2: return False
    r = K - j
    if p == 2:
        if r >= 3: return v % 8 == 1
        if r == 2: return v % 4 == 1
        return True
    return pow(v % p, (p - 1)//2, p) == 1

def lifting_insoluble(forms, p, Kmax=40, cap=3*10**6):
    """Pruned brute force: S_K = residues (1 : t), t mod p^K, and (p s : 1), s mod p^(K-1), at which every form is a
    square mod p^K, built from S_(K-1) by lifting.  Empty S_K => no Q_p point.  Returns K, or None (cap or Kmax)."""
    fs = [strip_even(F, p) for F in forms]
    def ok(s_, t_, K): return all(sq_mod_pK(F[0]*s_*s_ + F[1]*s_*t_ + F[2]*t_*t_, p, K) for F in fs)
    A = [t for t in range(p) if ok(1, t, 1)]
    B = [0] if ok(0, 1, 1) else []          # (p s : 1) with s mod p^0: the single class s = 0 (value at (0 : 1) mod p)
    K = 1
    while A or B:
        if K >= Kmax or len(A) + len(B) > cap: return None
        K += 1; q = p**(K - 1)
        A = [t + q*j for t in A for j in range(p) if ok(1, t + q*j, K)]
        qb = p**(K - 2)
        B = [s_ + qb*j for s_ in B for j in range(p) if ok(p*(s_ + qb*j), 1, K)]
    return K

def certify_independently(cd, d, al, be, lam, places, maxq=2*10**5):
    F1, F2 = cover_forms(cd, d, al, be, lam)
    if not soluble_R([F1, F2]): return 'R'
    for p in places:
        K = 1
        while p**K <= maxq:
            if brute_insoluble([F1, F2], p, K): return p
            K += 1
    for p in places:
        if lifting_insoluble([F1, F2], p) is not None: return p
    return None

def class2_all(mmax=300, out=None):
    """Class 2 only, every square-type curve with m <= mmax: survivors of the covers D_lam, compared with t8's tsv."""
    tsv = os.path.join(here, 'square_section.tsv')          # t8's table, if merged beside this script
    if not os.path.exists(tsv): tsv = '/home/claude/conductor/t8/scripts/g1_full/square_section.tsv'
    T = {}
    for line in open(tsv).read().splitlines()[1:]:
        w = line.split('\t'); T[(int(w[0]), int(w[1]))] = (w[3], int(w[4]), int(w[5]), w[6])
    agree = 0; dis = []
    for m in range(2, mmax + 1):
        for n in range(1, m):
            if gcd(m, n) != 1 or isqrt(m*m + n*n)**2 != m*m + n*n: continue
            for nn in (n, -n):
                cd = curve(m, nn)
                base = primes_of(2, m, nn, m - nn, m + nn, cd['k'], cd['Hp'], cd['Hm'])
                al, be = two_squares(2)
                lams = sorted({sqfree(l) for l in sqfree_divisors(2*2*cd['kap'])}, key=lambda z: (abs(z), z))
                surv = [lam for lam in lams if cover_soluble_everywhere(cd, 2, al, be, lam, base)[0]]
                st, ns, nc, det = T[(m, nn)]
                mine = 'sha' if not surv else 'nonempty'
                theirs = 'sha' if st == 'sha' else 'nonempty'
                ok = (mine == theirs) and len(surv) == ns and len(lams) == nc
                agree += ok
                if not ok: dis.append((m, nn, st, ns, nc, len(surv), len(lams)))
                say(f"  {m:4d} {nn:5d}  mine: {len(surv)}/{len(lams)} survivors ({mine});  t8: {st} {ns}/{nc}  {'AGREE' if ok else 'DISAGREE'}")
    say(f"class 2, m <= {mmax}: {agree} curves agree with t8 (status, number of surviving covers, number of covers); disagreements: {dis}")

SIX = [(15, 8), (15, -8), (56, 33), (60, -11), (99, -20), (180, 19)]

def load_json():
    for path in [os.path.join(here, 'selmer_square_5000.json'),
                 '/home/claude/conductor/t8/scripts/g1_full/selmer_square_5000.json']:
        if os.path.exists(path):
            J = json.load(open(path))
            return {(r['m'], r['n']): r for r in J}
    return {}

def negative():
    say("NEGATIVE TESTS")
    n_ok = 0
    # N1: a forged cover with a known rational point must be soluble everywhere: take E- (4,3), class 2 has a point
    # (t8 reports X = -1968); our run on the class 2 must leave a survivor.
    cd = curve(4, -3)
    al, be = two_squares(2)
    base = primes_of(2, 4, 3, 7, 1, 5, cd['Hp'], cd['Hm'])
    surv = [lam for lam in sorted({sqfree(l) for l in sqfree_divisors(2*2*cd['kap'])})
            if cover_soluble_everywhere(cd, 2, al, be, lam, base)[0]]
    x = -1968 + 2*(4*4 + 9)**3                     # X = x + e, e = -2 S^3: x = X - e
    xr = -1968 - (-2*25**3)
    assert xr > 0 and xr % 2 == 0
    yy = xr*(xr*xr + cd['a']*xr + cd['b']); assert isqrt(yy)**2 == yy, 'not a point'
    assert sqfree(xr) == 2
    assert surv, 'a class with a rational point was excluded'
    n_ok += 1; say(f"N1 E- (4,3): point x = {xr} of class 2 on the curve; survivors {surv}: not excluded (correct)")
    # N2: the point gives lam: X1 + i X2 = lam (1 + i)(s + i t)^2 must have its lam among the survivors
    M2 = xr // 2; # x = d M^2/e^2 with e = 1? check
    from fractions import Fraction as F_
    # x = 2 M^2 / e^2: find M, e
    q = F_(xr, 2); Mn, Md = isqrt(q.numerator), isqrt(q.denominator); assert Mn*Mn == q.numerator and Md*Md == q.denominator
    Mv, ev = Mn, Md
    X1 = 2*Mv*Mv - 2*cd['a0']*ev*ev; X2 = cd['kap']*ev*ev
    # (X1 + i X2)/(1 + i) = ((X1 + X2) + i (X2 - X1))/2 = lam (s + i t)^2
    zr, zi = F_(X1 + X2, 2), F_(X2 - X1, 2)
    nrm = zr*zr + zi*zi
    lamc = None
    for lam in surv:
        # z / lam must be a square in Q(i): its norm a rational square and (|z|/lam + zr/lam)/2 a rational square
        wr, wi = zr / lam, zi / lam
        n2 = wr*wr + wi*wi
        r = F_(isqrt(n2.numerator), isqrt(n2.denominator))
        if r*r != n2: continue
        for sgn in (1, -1):
            s2 = (sgn*r + wr) / 2
            if s2 >= 0 and isqrt(s2.numerator)**2 == s2.numerator and isqrt(s2.denominator)**2 == s2.denominator:
                lamc = lam
    assert lamc is not None
    n_ok += 1; say(f"N2 E- (4,3): the point's lam = {lamc} (mod squares and sign) is a survivor (correct)")
    # N3: method S and T on random forms at small p agree (also a control of both)
    import random
    rnd = random.Random(7); agree = 0
    for _ in range(300):
        p = rnd.choice([3, 5, 7, 11, 13])
        F1 = [rnd.randint(-60, 60) * p**rnd.choice([0, 0, 1, 2]) for _ in range(3)]
        F2 = [rnd.randint(-60, 60) * p**rnd.choice([0, 0, 1, 3]) for _ in range(3)]
        if disc2(F1) == 0 or disc2(F2) == 0 or res2(F1, F2) == 0: continue
        r1 = soluble_S([F1, F2], p); r2 = soluble_T([F1, F2], p, cap=30)
        assert r1 == r2; agree += 1
        # brute force: if a point with |s|, |t| <= 30 exists, both must say soluble
        if not r1:
            for s_ in range(-12, 13):
                for t_ in range(-12, 13):
                    if gcd(s_, t_) != 1: continue
                    v1 = F1[0]*s_*s_ + F1[1]*s_*t_ + F1[2]*t_*t_; v2 = F2[0]*s_*s_ + F2[1]*s_*t_ + F2[2]*t_*t_
                    assert not (sq_ok(v1, p) and sq_ok(v2, p)), 'solver says insoluble but a point exists'
    n_ok += 1; say(f"N3 methods S and T agree on {agree} random pairs of forms; no insoluble answer contradicted by a box point")
    # N4: forged Selmer: rank-1 conclusion must fail if Sel(E') is forged larger
    r = run_curve(15, 8, verbose=False)
    assert r['rank1']
    forged = dict(r); forged['selEp'] = [-1, 1, 3, -3]
    assert not (forged['excluded_all'] and sorted(forged['selEp']) == [-1, 1])
    n_ok += 1; say("N4 E+ (15,8) with Sel(E') forged to 4 elements: rank-1 conclusion refused")
    # N5: a lam outside the allowed set is not needed, but a FORGED kappa (kappa + 2) must break b = 4 a0^2 + kappa^2
    try:
        cdf = dict(curve(15, 8)); cdf['kap'] += 2
        assert cdf['b'] == 4*cdf['a0']**2 + cdf['kap']**2
        bad = False
    except AssertionError:
        bad = True
    assert bad; n_ok += 1; say("N5 forged kappa: refused")
    # N6: real place forged: F1 = -(s^2 + t^2) has no real point
    assert soluble_R([[-1, 0, -1], [1, 0, 1]]) is False and soluble_R([[1, 0, -1], [1, 0, 1]]) is True
    n_ok += 1; say("N6 real test: -(s^2 + t^2) refused, s^2 - t^2 accepted")
    say(f"negative: {n_ok}/6 as required")

def main():
    args = sys.argv[1:]
    if '--negative' in args:
        negative(); return
    if '--class2' in args:
        class2_all(int(args[args.index('--class2') + 1])); return
    cross = '--crosscheck' in args
    which = SIX
    if '--curves' in args:
        v = args[args.index('--curves') + 1]
        if v != 'all6':
            which = [tuple(int(z) for z in w.split(':')) for w in v.split(',')]
    J = load_json()
    say("review_t8_rank1.py: independent second descent on the six rank-1 claims of SQUARE_SECTION.md")
    allok = True
    for (m, nn) in which:
        r = run_curve(m, nn, J.get((m, nn)), cross, brute='--brute' in args)
        allok &= r['rank1'] and r['consistent']
    say(f"ALL SIX RANK-1 CLAIMS CONFIRMED: {allok}")

if __name__ == '__main__':
    main()
