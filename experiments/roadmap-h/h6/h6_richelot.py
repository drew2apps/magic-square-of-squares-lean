"""H6 (part 3): an upper bound for the rank of the Jacobian of a genus-2 curve y^2 = G1 G2 G3 over Q (three rational
factors of degree <= 2), by descent via a Richelot isogeny (EXPERIMENT, written by an AI model; not read by a person;
nothing here changes the census).
C: y^2 = G1 G2 G3, C': Delta y^2 = H1 H2 H3 with H1 = [G2,G3], H2 = [G3,G1], H3 = [G1,G2] ([A,B] = A'B - AB') and
Delta the determinant of the coefficients of G1, G2, G3 (Bost-Mestre; the twist of C' is fixed by checking that C and
C' have the same number of points over F_p and F_{p^2} at good primes: isogenous Jacobians have the same L-polynomial).
phi: J(C) -> J(C') has kernel {0, T1, T2, T3} (T_i = the divisor of the roots of G_i, all rational), phi' its dual.
delta : J(C)(Q)/phi'J(C')(Q) -> (Q*/Q*^2)^2, P1 + P2 - (fibre) |-> (G1(x1)G1(x2), G2(x1)G2(x2)); delta' likewise on
J(C') with H1, H2.  The Selmer groups are cut out by the local images W_p, W'_p at p in S = {2, primes of the
discriminant and leading coefficients} and infinity; outside S both are the unramified classes.
Local images are found by sampling Q_p-points (and, if needed, pairs of conjugate points over a quadratic extension;
PARI's nfislocalpower decides whether such a pair is rational).  They are CERTIFIED by the product formula
|W_p| |W'_p| = |H^1(Q_p, (Z/2)^2)| = 16 (p odd), 64 (p = 2), 4 (R): the two images are exact annihilators under the
pairing given by the Weil pairing (Tate local duality), so the sampled subgroups, which lie inside the true images,
are the true images exactly when the product of their sizes reaches that number.
Then rank J(C)(Q) <= dim Sel(delta) + dim Sel(delta') - 2 (J(C)(Q)[2] contains T1, T2).
Theory quoted, not proved here: the Richelot construction, the description of the Kummer maps by (G1, G2), local
duality, unramified images at odd primes of good reduction.  The arithmetic is exact (Python integers and fractions),
except for nfislocalpower where it is used.
usage: h6_richelot.py m n [--write]   (reads h6_split_<m>_<n>.json; writes h6_richelot_<m>_<n>.txt)
       h6_richelot.py --controls [--write]   (h6_richelot_controls.txt)    h6_richelot.py --negative [--write]"""
import sys, os, json, random, itertools
from fractions import Fraction as Fr
import sympy, cypari2
pari = cypari2.Pari(); pari.allocatemem(10**9)
HERE = os.path.dirname(os.path.abspath(__file__))
X = sympy.Symbol('x')

def padd(a, b):
    n = max(len(a), len(b)); return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]
def pmul(a, b):
    r = [0]*(len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i+j] += x*y
    return r
def pneg(a): return [-c for c in a]
def pder(a): return [i*a[i] for i in range(1, len(a))] or [0]
def ev(a, x):
    r = 0
    for c in reversed(a): r = r*x + c
    return r
def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0: a.pop()
    return a
def wr(A, B): return trim(padd(pmul(pder(A), B), pneg(pmul(A, pder(B)))))

# ---- square classes ----
def vp(a, p):
    a = Fr(a); assert a != 0
    v = 0; nu, de = a.numerator, a.denominator
    while nu % p == 0: nu //= p; v += 1
    while de % p == 0: de //= p; v -= 1
    return v, nu * de                                           # unit part (up to squares: nu/de ~ nu*de)
def cls(a, p):
    """coordinates of a in Q_p*/Q_p*^2 (p = 0 is the real place)"""
    a = Fr(a)
    if p == 0: return (1 if a < 0 else 0,)
    v, u = vp(a, p)
    if p == 2:
        u %= 8; return (v % 2, 1 if u % 4 == 3 else 0, 1 if u in (3, 5) else 0)
    return (v % 2, 0 if pow(u % p, (p - 1)//2, p) == 1 else 1)
def issq_p(a, p):
    if a == 0: return False
    return all(c == 0 for c in cls(a, p))
def hdim(p): return 1 if p == 0 else 3 if p == 2 else 2

# ---- F_2 linear algebra ----
def rank2(vs):
    rows = [int(''.join(map(str, v)), 2) for v in vs if any(v)]
    basis = []
    for r in rows:
        for b in basis: r = min(r, r ^ b)
        if r: basis.append(r)
    return len(basis), basis
def span_size(vs): return 2**rank2(vs)[0]

def nullspace_dim(M, ncols):
    """dimension of {e in F_2^ncols : M e = 0}, M a list of rows (lists of 0/1 of length ncols)"""
    rows = [int(''.join(map(str, r)), 2) for r in M]
    basis = []
    for r in rows:
        for b in basis: r = min(r, r ^ b)
        if r: basis.append(r)
    return ncols - len(basis)

# ---- point counts (for the twist of C') ----
def count(f, p):
    N = 0
    for x in range(p):
        a = ev(f, x) % p
        N += 1 if a == 0 else (2 if pow(a, (p-1)//2, p) == 1 else 0)
    d = len(trim([c % p for c in f])) - 1
    lc = trim([c % p for c in f])[-1]
    if d == 6: N += 2 if pow(lc, (p-1)//2, p) == 1 else 0
    elif d == 5: N += 1
    else: return None
    return N

class Rich:
    def __init__(s, G):
        s.G = [trim(g) for g in G]
        s.f = trim(pmul(s.G[0], pmul(s.G[1], s.G[2])))
        assert len(s.f) - 1 in (5, 6)
        g3 = [[(g + [0, 0, 0])[k] for k in range(3)] for g in s.G]
        s.Delta = sympy.Matrix(g3).det()
        assert s.Delta != 0, 'Delta = 0: the Jacobian splits'
        s.H = [wr(s.G[1], s.G[2]), wr(s.G[2], s.G[0]), wr(s.G[0], s.G[1])]
        base = trim(pmul(s.H[0], pmul(s.H[1], s.H[2])))
        # choose the twist (Delta or -Delta) by point counts
        s.fd = None
        for sg in (1, -1):
            fd = [sg * int(s.Delta) * c for c in base]
            primes = [p for p in sympy.primerange(3, 400) if s.good(p) and s.good(p, fd)][:12]
            assert len(primes) >= 8, 'too few good primes for C and its dual'
            if all(count(s.f, p) == count(fd, p) for p in primes): s.fd = fd; s.Hd = [[sg*int(s.Delta)*c for c in s.H[0]], s.H[1], s.H[2]]; break
        assert s.fd is not None, 'no twist of the Richelot dual matches the point counts'
        s.checked = primes
    def good(s, p, f=None):
        P = sympy.Poly(sum(c * X**i for i, c in enumerate(s.f if f is None else f)), X)
        return P.LC() % p != 0 and sympy.discriminant(P) % p != 0
    def badS(s):
        P = sympy.Poly(sum(c * X**i for i, c in enumerate(s.f)), X)
        n = abs(int(sympy.discriminant(P) * P.LC()))
        return sorted(set([2]) | set(int(q) for q in sympy.factorint(n)))

NFC = {}
def near_roots(G, p):
    """the roots mod p of the factors G_i (each made p-primitive first), lifted p-adically to precision 12 by PARI's
    polrootspadic where that works: only points to sample near (a finder; certification is by the product formula)"""
    out = set()
    for g in G:
        g = [int(c) for c in trim(g)]
        if len(g) < 2: continue
        from math import gcd
        c = 0
        for v in g: c = gcd(c, v)
        while c % p == 0: g = [v // p for v in g]; c //= p
        P = pari('Pol(%s, x)' % str(list(reversed(g))))
        if pari.poldegree(pari.Mod(1, p) * P) < 1: continue
        for r in pari.polrootsmod(P, p): out.add(int(pari.lift(r)))
        try:
            for r in pari.polrootspadic(P, p, 12):
                if pari.valuation(r, p) >= 0: out.add(int(pari.lift(pari.Mod(pari.truncate(r), p**12))))
        except cypari2.PariError: pass
    return sorted(out)
def sample_image(G, f, p, rnd, target_extra=None, tries=4000, quad=False):
    """subgroup of (Q_p*/Q_p*^2)^2 spanned by delta of divisors P - P0 (Q_p-points) and of rational conjugate pairs"""
    vals = []
    def addx(x):
        y2 = ev(f, x)
        if y2 == 0: return
        if p == 0:
            if y2 < 0: return
        elif not issq_p(y2, p): return
        a, b = ev(G[0], x), ev(G[1], x)
        if a == 0 or b == 0: return
        vals.append(cls(a, p) + cls(b, p))
    roots = near_roots(G, p) if p else []
    for _ in range(tries):
        if p == 0:
            x = Fr(rnd.randrange(-10**6, 10**6), rnd.randrange(1, 10**4))
        elif roots and rnd.random() < 0.5:      # near a root mod p of a factor (where the classes change)
            x = Fr(rnd.choice(roots)) + Fr(rnd.randrange(-p**4, p**4)) * Fr(p)**rnd.randrange(1, 9)
        else:
            e = rnd.randrange(-4, 5); x = Fr(rnd.randrange(-p**8, p**8)) * Fr(p)**e
        addx(x)
    # points at infinity (degree 6, leading coefficient a square): delta = leading coefficients
    if len(f) == 7 and (issq_p(f[-1], p) if p else f[-1] > 0):
        a, b = trim(G[0])[-1], trim(G[1])[-1]
        if len(trim(G[0])) == 3 and len(trim(G[1])) == 3: vals.append(cls(a, p) + cls(b, p))
    pts = list(vals)
    gens = [tuple((u + w) % 2 for u, w in zip(v, pts[0])) for v in pts[1:]] if pts else []
    if quad and p != 0:
        # conjugate pairs x1, x2 = a +- b sqrt(D), D running over small representatives of the non-trivial classes of
        # Q_p*/Q_p*^2 (so nfinit is only ever called on y^2 - D with D small)
        if p == 2: Ds = [-1, 2, -2, 3, -3, 6, -6]
        else:
            u = next(k for k in range(2, p) if pow(k, (p - 1)//2, p) == p - 1)
            Ds = [u, p, u*p]
        for _ in range(tries // 4):
            if p < 100 and rnd.random() < 0.5:
                # small p: a generic pair x^2 - sx x + nx (nfinit on its discriminant, which is small enough here)
                e = rnd.randrange(-2, 3)
                sx = Fr(rnd.randrange(-p**6, p**6)) * Fr(p)**e; nx = Fr(rnd.randrange(-p**6, p**6)) * Fr(p)**(2*e)
                D = sx*sx - 4*nx
                if D == 0 or issq_p(D, p): continue
                Dn = D.numerator * D.denominator
                nf = pari.nfinit('y^2 - (%d)' % Dn)
                x1 = pari('Mod((%s) + y/(%d), y^2 - (%d))/2' % (str(sx), D.denominator, Dn))
            else:
                Dn = rnd.choice(Ds)
                if Dn not in NFC: NFC[Dn] = pari.nfinit('y^2 - (%d)' % Dn)
                nf = NFC[Dn]
                e = rnd.randrange(-2, 3)
                a = Fr(rnd.randrange(-p**6, p**6)) * Fr(p)**e; b = Fr(rnd.randrange(1, p**6)) * Fr(p)**rnd.randrange(-2, 3)
                if roots and rnd.random() < 0.5:
                    a = Fr(rnd.choice(roots)) + Fr(rnd.randrange(-p**4, p**4)) * Fr(p)**rnd.randrange(1, 7)
                    b = Fr(rnd.randrange(1, p**4)) * Fr(p)**rnd.randrange(0, 7)
                sx, nx = 2*a, a*a - b*b*Dn
                x1 = pari('Mod((%s) + (%s)*y, y^2 - (%d))' % (str(a), str(b), Dn))
            fx = sum(pari(str(c)) * x1**i for i, c in enumerate(f))
            if fx == 0: continue
            pr = pari.idealprimedec(nf, p)[0]
            if not int(pari.nfislocalpower(nf, pr, fx, 2)): continue
            h = [nx, -sx, 1]
            r1 = sympy.resultant(sum(c*X**i for i, c in enumerate(h)), sum(c*X**i for i, c in enumerate(G[0])), X)
            r2 = sympy.resultant(sum(c*X**i for i, c in enumerate(h)), sum(c*X**i for i, c in enumerate(G[1])), X)
            if r1 == 0 or r2 == 0: continue
            gens.append(cls(Fr(r1), p) + cls(Fr(r2), p))
    return gens

def local_images(R, p, rnd):
    tot = {0: 4, 2: 64}.get(p, 16)
    gW, gW2 = [], []
    for rnd_quad in (False, True, True, True, True, True):
        gW += sample_image(R.G, R.f, p, rnd, quad=rnd_quad)
        gW2 += sample_image(R.Hd, R.fd, p, rnd, quad=rnd_quad)
        a, b = span_size(gW), span_size(gW2)
        assert a * b <= tot, ('local images too large: the descent maps are wrong', p, a, b, tot)
        if a * b == tot: return rank2(gW)[1], rank2(gW2)[1], (a, b), rnd_quad
    return None, None, (a, b), True

def local_coords(a, p):
    return cls(Fr(a), p)

def selmer(S, W_by_p):
    """basis of the Selmer group: pairs (a1, a2) in Q(S,2)^2 whose local classes lie in W_p at every p in S and at
    infinity.  Returned as lists of exponent vectors (e1, e2) over the basis (-1, primes of S)."""
    B = [-1] + S; nb = len(B)
    rows = []
    for p in [0] + S:
        k = hdim(p); dimv = 2*k
        Wv = [[(b >> (dimv - 1 - i)) & 1 for i in range(dimv)] for b in W_by_p[p]]
        ann = [c for c in itertools.product((0, 1), repeat=dimv)
               if any(c) and all(sum(x*y for x, y in zip(c, w)) % 2 == 0 for w in Wv)]
        ann = [[(b >> (dimv - 1 - i)) & 1 for i in range(dimv)] for b in rank2(ann)[1]]
        loc = [local_coords(b, p) for b in B]
        for c in ann:
            row = []
            for comp in (0, 1):
                for j in range(nb): row.append(sum(c[comp*k + i] * loc[j][i] for i in range(k)) % 2)
            rows.append(row)
    # nullspace over F_2
    ncols = 2*nb; piv = {}; R = [r[:] for r in rows]
    red = []
    for r in R:
        for (pc, pr) in red:
            if r[pc]: r = [(x + y) % 2 for x, y in zip(r, pr)]
        if any(r):
            pc = r.index(1)
            red = [(c0, [(x + y) % 2 for x, y in zip(rr, r)] if rr[pc] else rr) for c0, rr in red]
            red.append((pc, r))
    pivots = {pc: r for pc, r in red}
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for fc in free:
        v = [0]*ncols; v[fc] = 1
        for pc, r in pivots.items():
            if r[fc]: v[pc] = 1
        basis.append(v)
    for v in basis:
        for r in rows: assert sum(x*y for x, y in zip(r, v)) % 2 == 0
    return B, basis

def in_W(vec, W, p):
    k = hdim(p); dimv = 2*k
    n, b = rank2([[(w >> (dimv - 1 - i)) & 1 for i in range(dimv)] for w in W])
    return rank2([[(w >> (dimv - 1 - i)) & 1 for i in range(dimv)] for w in W] + [list(vec)])[0] == n

def prod_at_roots(A, B):
    """prod of B over the roots of A (A of degree 1 or 2; a root at infinity of a degree-1 A is not included):
    Res(A, B) / lc(A)^deg B, exact"""
    A, B = trim(A), trim(B)
    r = sympy.resultant(sum(c*X**i for i, c in enumerate(A)), sum(c*X**i for i, c in enumerate(B)), X)
    return Fr(int(r)) / Fr(A[-1])**(len(B) - 1)

def torsion_images(G):
    """delta(T_i), T_i = the divisor of the roots of G_i (degree 2 each): coordinates (G1(T_i), G2(T_i)) in Q*/Q*^2,
    with the vanishing coordinate replaced through G1 G2 G3 = y^2.  Degree-1 factors (a root at infinity) are skipped."""
    out = []
    for i in range(3):
        if len(trim(G[i])) != 3: continue
        j, k = [t for t in range(3) if t != i]
        vals = {j: prod_at_roots(G[i], G[j]), k: prod_at_roots(G[i], G[k])}
        vals[i] = vals[j] * vals[k]
        out.append((vals[0], vals[1]))
    return out

def factor_groups(f):
    """all ways of writing f as c G1 G2 G3 with rational G_i of degree <= 2 (degree 5: one G_i of degree 1)"""
    P = sympy.Poly(sum(c * X**i for i, c in enumerate(f)), X)
    c, fl = sympy.factor_list(P.as_expr())
    facs = []
    for g, e in fl:
        facs += [sympy.Poly(g, X)] * e
    if any(g.degree() > 2 for g in facs): return []
    out = set()
    for perm in itertools.permutations(range(len(facs))):
        groups = [[], [], []]; ok = True; gi = 0
        # greedy fill by degree in permutation order
        for idx in perm:
            while gi < 3 and sum(facs[j].degree() for j in groups[gi]) + facs[idx].degree() > 2: gi += 1
            if gi == 3: ok = False; break
            groups[gi].append(idx)
        if not ok or any(not g for g in groups): continue
        key = tuple(sorted(tuple(sorted(g)) for g in groups))
        out.add(key)
    res = []
    for key in sorted(out):
        Gs = []
        for g in key:
            q = sympy.Poly(1, X)
            for j in g: q = q * facs[j]
            Gs.append([int(v) for v in reversed(q.all_coeffs())])
        Gs[2] = [int(c) * v for v in Gs[2]] if c == int(c) else None
        if Gs[2] is None: continue
        if sorted(len(trim(g)) - 1 for g in Gs) in ([2, 2, 2], [1, 2, 2]) and len(trim(pmul(Gs[0], pmul(Gs[1], Gs[2])))) == len(trim(f)):
            res.append(Gs)
    return res

def bound(f, seed=1, verbose=False):
    best = None; out = []
    for G in factor_groups(f):
        try: R = Rich(G)
        except AssertionError as e:
            out.append('grouping %s: %s' % (G, e.args[0])); continue
        S = R.badS(); rnd = random.Random(seed)
        W, W2, ok = {}, {}, True; quadused = False
        for p in [0] + S:
            a, b, sz, q = local_images(R, p, rnd)
            if a is None: ok = False; out.append('grouping %s: local images at p = %d not certified (sizes %s)' % (G, p, sz)); break
            W[p], W2[p] = a, b; quadused |= q
        if not ok: continue
        B, sel1 = selmer(S, W); _, sel2 = selmer(S, W2)
        s1, s2 = len(sel1), len(sel2)
        # 2^r |X(Q)[2]| = |X/phi'Y| |Y/phi X| / [ker phi'(Q) : ker phi'(Q) cap phi X(Q)], and the index is the size of
        # delta'(ker phi'(Q)) (delta' is injective on Y/phi X); symmetrically with X and Y exchanged
        def tdim(Gs):
            vs = []
            for a, b in torsion_images(Gs):
                v = []
                for p in [0] + S: v += list(cls(a, p) + cls(b, p))
                for p in sympy.primefactors(abs(a.numerator*a.denominator*b.numerator*b.denominator)):
                    if p not in S: v += list(cls(a, p) + cls(b, p))
                vs.append(v)
            # dimension of the span in Q*/Q*^2 x Q*/Q*^2: by the square classes of the pairs
            sq = [(sympy.factorint(abs(a.numerator*a.denominator)), sympy.factorint(abs(b.numerator*b.denominator)), a < 0, b < 0) for a, b in torsion_images(Gs)]
            primes = sorted(set(q for x in sq for d in x[:2] for q in d))
            rows = [[int(x[2])] + [x[0].get(q, 0) % 2 for q in primes] + [int(x[3])] + [x[1].get(q, 0) % 2 for q in primes] for x in sq]
            return rank2(rows)[0]
        t1, t2 = tdim(R.G), tdim(R.Hd)
        r = s1 + s2 - 2 - max(t1, t2)
        # global points (a finder: PARI's hyperellratpoints; each checked exactly): their images must be in Sel
        npts = 0
        for curve, GG, WW in ((R.f, R.G, W), (R.fd, R.Hd, W2)):
            pts = [(Fr(int(pari.numerator(P[0])), int(pari.denominator(P[0]))), Fr(int(pari.numerator(P[1])), int(pari.denominator(P[1]))))
                   for P in pari.hyperellratpoints(pari('Pol(%s, x)' % str(list(reversed(curve)))), 300)]
            pts = [(x, y) for x, y in pts if y*y == ev(curve, x) and ev(GG[0], x) != 0 and ev(GG[1], x) != 0]
            for (x, y) in pts[1:]:
                x0 = pts[0][0]
                a, b = ev(GG[0], x) * ev(GG[0], x0), ev(GG[1], x) * ev(GG[1], x0)
                for p in [0] + S:
                    assert in_W(cls(a, p) + cls(b, p), WW[p], p), ('a global point is outside the local image', p)
                npts += 1
        out.append('grouping %s: S = %s; dim Sel(phi-hat side) %d, dim Sel(phi side) %d; dim delta(T) %d, dim delta-dual(T-dual) %d; rank J(Q) <= %d%s; %d global point differences checked in Sel' % (
            [len(trim(g)) - 1 for g in G], S, s1, s2, t1, t2, r, '' if not quadused else ' (quadratic points used)', npts))
        best = r if best is None else min(best, r)
    return best, out

def main(m, n):
    pieces = json.load(open(os.path.join(HERE, 'h6_split_%d_%d.json' % (m, n))))
    lines = ['H6: Richelot descent on the genus-2 pieces of (%d,%d) (h6_split). Written by an AI model; not read by a person. Experiment, not counted.' % (m, n),
             'For each piece: the groupings of its factors, the Selmer dimensions and the rank bound; the multiplicity of the piece in Jac(cover).']
    seen = {}
    for p in pieces:
        if p['genus'] != 2: continue
        f = trim(p['form'])
        b, out = bound(f)
        lines.append('%-14s %-14s best bound: %s' % (p['S'], p['path'], b))
        lines += ['    ' + o for o in out]
        sys.stderr.write(lines[-1 - len(out)] + '\n')
    return lines

def controls():
    """y^2 = (x^2 + x + 1)(x^2 - x + 1)(x^2 + c) is even, so its Jacobian is isogenous to E1 x E2 with
    E1: y^2 = (w^2 + w + 1)(w + c), E2: y^2 = w (w^2 + w + 1)(w + c); but this grouping has Delta = 2c - 2 != 0, so the
    Richelot descent runs.  Its bound must be at least PARI's rank of E1 plus that of E2 (small conductors, bounds equal).
    The images of the small global points must lie in the Selmer groups (asserted inside bound())."""
    lines = []
    for c in (2, 3, -2, 5, -3, 7, 10, -6, 13, -11):
        G = [[1, 1, 1], [1, -1, 1], [c, 0, 1]]
        f = trim(pmul(G[0], pmul(G[1], G[2])))
        e1 = pari.ellinit(pari('ellfromeqn(y^2 - (x^2+x+1)*(x+%d))' % c))
        e2 = pari.ellinit(pari('ellfromeqn(y^2 - x*(x^2+x+1)*(x+%d))' % c))
        r1, r2 = pari.ellrank(e1, 2), pari.ellrank(e2, 2)
        lo, hi = int(r1[0]) + int(r2[0]), int(r1[1]) + int(r2[1])
        bd, out = bound(f)
        lines.append('c = %d: PARI rank E1 + rank E2 in [%d,%d]; Richelot bound %s: %s' % (c, lo, hi, bd,
                     'consistent' if bd is not None and bd >= hi else 'INCONSISTENT' if bd is not None else 'not certified'))
        lines += ['    ' + o for o in out]
    return lines

def negative():
    """forged inputs: (1) a wrong descent map (G1 replaced by G1 + 1) must break the local product formula or the
    point check; (2) a forged local image (one generator dropped) must leave the product below the target, so it is
    not certified; (3) a wrong twist of the dual must be refused by the point counts."""
    lines = []
    G = [[1, 0, 1], [2, 1, 1], [-3, 1, 2]]
    f = trim(pmul(G[0], pmul(G[1], G[2])))
    R = Rich(G)
    # (1)
    R1 = Rich(G); R1.G = [[2, 0, 1], G[1], G[2]]
    bad = None
    for p in [0] + R.badS():
        try: local_images(R1, p, random.Random(1))
        except AssertionError as e: bad = (p, e.args[0][0]); break
    lines.append('wrong descent map: %s' % ('refused at p = %d (%s)' % bad if bad else 'NOT refused by the product formula'))
    # (2)
    W, W2, sz, q = local_images(R, 5, random.Random(1))
    tot = 16
    lines.append('dropped generator at p = 5 (a bad prime): product %d instead of %d: %s' % (2**(len(W) - 1) * 2**len(W2), tot,
                 'not certified' if 2**(len(W) - 1) * 2**len(W2) < tot else 'CERTIFIED'))
    # (3)
    base = trim(pmul(R.H[0], pmul(R.H[1], R.H[2])))
    fw = [-int(R.Delta) * c for c in base]
    ps = [p for p in sympy.primerange(3, 400) if R.good(p) and R.good(p, fw)][:12]
    lines.append('wrong twist of the dual (-Delta): %s' % ('refused by the point counts' if any(count(f, p) != count(fw, p) for p in ps) else 'ACCEPTED'))
    return lines

if __name__ == '__main__':
    argv = sys.argv
    if '--controls' in argv: L = controls(); fn = 'h6_richelot_controls.txt'
    elif '--negative' in argv: L = negative(); fn = 'h6_richelot_negative.txt'
    else:
        m, n = int(argv[1]), int(argv[2]); L = main(m, n); fn = 'h6_richelot_%d_%d.txt' % (m, n)
    s = '\n'.join(L) + '\n'; fn = os.path.join(HERE, fn)
    if '--write' in argv: open(fn, 'w').write(s)
    else:
        sys.stdout.write(s)
        if os.path.exists(fn): print('# matches' if open(fn).read() == s else '# DIFFERS from ' + os.path.basename(fn))
