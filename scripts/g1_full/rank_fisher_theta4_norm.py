#!/usr/bin/env python3
"""FINDER, not certifying: points of the two cone conics of the level-3 curve C3' (rank_fisher5.py) of E+ of
(945,248), E' side, over L = Q(sqrt 5820647) (class group [4, 2, 2]) and L' = Q(sqrt -5820647) ([604, 2, 2]), for
rank_fisher_theta4_gen.py; and the reason for the level-3 twist xi = 17794 used there.  The lattice finder of
rank_fisher_theta4_find.py / conic_lattice_find.py found no point on these conics (60000 short vectors, ~9e6 pairs).

    python3 rank_fisher_theta4_norm.py                 re-check the committed hints exactly (rank_fisher5.analyse)
                                                       and two refusals; --write: rank_fisher_theta4_norm.txt
    python3 rank_fisher_theta4_norm.py --find          find the points again (about 1-2 minutes) and compare;
                                                       --find --write: write them to rank_fisher_theta4_gen_hints.json
    python3 rank_fisher_theta4_norm.py --obstruction   the dual obstruction of the level-3 twists xi = 1271 s,
                                                       s = 1, 14, 17, 238 (finder points, sampled local classes);
                                                       --write: rank_fisher_theta4_norm_obstruction.txt

Written by an AI model (Claude); not read by a person.  Experiment, not counted.

WHAT RESTS ON THIS SCRIPT: NOTHING.  Its points are inputs (hints); rank_fisher_theta4_gen.py re-checks them
exactly through rank_fisher5.analyse (the point on its cone, tangency of l_1, l_3 along the ruling, f = N(l_1) and
f_dag = N(l_3) rational, q and c by exact linear algebra), and the level-4 curve C4 built from them is checked
locally soluble at every place of S and at R there.  The default mode of this script runs the same exact check.

THE FINDER (PARI is used with every assumption it likes: bnfinit flag 1 (class group and units under GRH, not
certified), rnfisnorm (GRH), floating point (mpmath) for lattice weights and enumeration radii only):
 (1) the conic h0 Z^2 + b20 e^2 + b11 e g + b02 g^2 = 0 over L (rank_fisher5 step (5); the V conic likewise)
     becomes X^2 - A Y^2 - B Z^2 = 0, A = -4 b20 h0, B = b11^2 - 4 b20 b02 (X = 2 b20 e + b11 g, Y = Z, Z = g);
 (2) A, B are reduced modulo squares with the class group: (A) = J I^2, I = K0 (alpha) with K0 a small ideal in
     the class of I (bnfisprincipal with forced generator), A -> A / alpha^2 (unit-balanced); likewise B;
 (3) a Legendre-type descent over L: for P | J_A, u = sqrt(B) v mod P cuts a rank-4 lattice of (u, v) in O_L^2 on
     which u^2 - B v^2 is divisible by J_A; among its short vectors (LLL, qfminim) the one whose cofactor has the
     smallest square-free part gives A' = (u^2 - B v^2) A / s^2, an equivalent conic (and a map back of points);
     the same for B; a random walk over both sides until N(J_A) N(J_B) <= 1e11 (5-10 steps here);
 (4) the small conic by rnfisnorm over L(sqrt A) (absolute discriminant about 1e22, a few seconds);
 (5) the point is large (hundreds to thousands of digits): it is made small with the quaternion algebra
     H = (A, B)_L.  xi = X0 + Y0 i + Z0 j is a zero divisor; the left ideal H xi meets O_L<1, i, j, k> in a rank-4
     lattice; for a short v there, the left ideal H conj(v) gives the point (z x - A w y, z y - w x, z^2 - A w^2)
     (v = x + y i + z j + w k after conjugation), small; done on the final conic, after the unwinding on the
     class-reduced conic, then mapped to the original conic and divided by a generator of the gcd ideal.
 Hints: U point 35-45 digits per coordinate (Z, g), V point 36-56; e, c about 90.

THE LEVEL-3 TWIST.  rank_fisher_curves.json's E'-side C3 of (945,248) has lam = -349, xi = 1271.  On this side
S_1 (dim 4, rank_indep: <-1, 14, 17, 97826>) is larger than S_2 = <-1, -5820647> (dim 2), and a level-3 lift lifts
to level 4 iff its Cassels-Tate pairing with every eta of S_1 vanishes (rank_fisher.py docstring (6); Poitou-Tate:
the Selmer set of C3' for its pushout form f0 is non-empty iff Ob(eta) = prod_v (f0(P_v), eta)_v = 1 for all eta in
S_1).  --obstruction computes Ob from sampled local points (a sampler: the classes of f0 at each place are a coset
of the local image, and the symbol is checked constant on the sampled coset):
    xi = 1271:            Ob(14) = +1, Ob(17) = -1, Ob(97826) = -1   (no level-4 twist: the json's C3 does not lift)
    xi = 17794 = 1271 * 14: Ob = 1 on S_1                             (used by rank_fisher_theta4_gen.py)
    xi = 21607 = 1271 * 17: Ob(14) = Ob(17) = -1;  xi = 302498 = 1271 * 238: Ob(14) = -1, Ob(17) = +1
 Ob is 1 on S_2 in every case (known classes and y0: Thm 2.1), and it moves by Theta_1(s, .) under xi -> xi s,
 Theta_1 alternating of rank 2 on S_1 / S_2, as the reading of Fisher Sect. 2 in rank_fisher.py predicts.  The
 twist xi = 17794 changes only the third equation of C3 (xi w^2 = l_X l_Y); D_lam, U, V, l_X, l_Y, p_X, p_Y are
 the json's.  It is stored in the hints file under the key '945,248 E+|E'|xi'.
"""
import os, sys, json, time, random, math
from fractions import Fraction
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
sys.set_int_max_str_digits(0)
import sympy as sp
import mpmath
import cypari2
pari = cypari2.Pari()
try:
    pari.allocatemem(2 * 10 ** 9)
except Exception:
    pass
import rank_fisher5 as R5

WRITE = '--write' in sys.argv
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(s)

HINTS_GEN = os.path.join(here, 'rank_fisher_theta4_gen_hints.json')
NAME, SIDE, XI = '945,248', "E'", 17794
t_ = sp.Symbol('t')

# ---------------------------------------------------------------- the quadratic field
class Fld:
    def __init__(self, D0):
        self.D0 = D0
        self.pol = pari('y^2 - (%d)' % D0)
        self.bnf = pari.bnfinit(self.pol, 1)
        self.nf = self.bnf
        zk = pari('nfinit(y^2 - (%d)).zk' % D0)
        self.w = pari.Mod(zk[1], self.pol)          # O_L = Z + Z w
        self.real = D0 > 0
        self._crep = {}
    def el(self, z):
        return pari.Mod(pari.lift(z), self.pol) if str(pari.type(z)) == 't_POLMOD' else pari.Mod(z, self.pol)
    def coords(self, z):
        v = pari.nfalgtobasis(self.nf, pari.lift(z))
        return [int(v[0]), int(v[1])]
    def from_coords(self, c):
        return pari.Mod(c[0], self.pol) + c[1] * self.w
    def emb(self, z):
        """embeddings (mpmath; real field: two values, the cancelling one from the exact norm)."""
        zl = pari.lift(z)
        a = Fraction(str(pari.polcoef(zl, 0, 'y'))); b = Fraction(str(pari.polcoef(zl, 1, 'y')))
        dig = max(len(str(abs(a.numerator))), len(str(abs(b.numerator))), 30) + 20
        with mpmath.workdps(dig):
            A_ = mpmath.mpf(a.numerator) / a.denominator; B_ = mpmath.mpf(b.numerator) / b.denominator
            s = mpmath.sqrt(abs(self.D0))
            if self.real:
                nq = a * a - b * b * self.D0; n = mpmath.mpf(nq.numerator) / nq.denominator
                if a * b >= 0:
                    p = A_ + B_ * s; out = [p, n / p if p != 0 else A_ - B_ * s]
                else:
                    m_ = A_ - B_ * s; out = [n / m_ if m_ != 0 else A_ + B_ * s, m_]
                return [mpmath.mpf(+x) for x in out]
            return [mpmath.mpc(A_, B_ * s)]

def factor_int(n, known=()):
    n = abs(int(n)); out = {}
    for p in known:
        while n % p == 0: n //= p; out[p] = out.get(p, 0) + 1
    if n > 1:
        fa = pari.factor(n)
        for i in range(len(fa[0])):
            p = int(fa[0][i]); assert pari.isprime(p)
            out[p] = out.get(p, 0) + int(fa[1][i])
    return out

def ideal_fact(K, z, known=()):
    n = abs(int(pari.norm(z)))
    res = []
    for p in sorted(factor_int(n, known)):
        for P in pari.idealprimedec(K.nf, p):
            e = int(pari.idealval(K.nf, pari.lift(z), P))
            if e: res.append((P, e))
    return res

def class_rep(K, cl):
    """an ideal of small norm in the class cl: the smallest degree-1 prime ideal of norm < 3000 in it, else idealred."""
    key = tuple(int(c) for c in cl)
    if key in K._crep: return K._crep[key]
    best = None
    for q in range(2, 3000):
        if not pari.isprime(q): continue
        for Pq in pari.idealprimedec(K.nf, q):
            if int(Pq[3]) == 1 and pari.bnfisprincipal(K.bnf, Pq, 0) == cl: best = Pq; break
        if best is not None: break
    if best is None:
        I = pari.idealhnf(K.nf, 1)
        for g_, c_ in zip(K.bnf.bnf_get_gen(), key): I = pari.idealmul(K.nf, I, pari.idealpow(K.nf, g_, c_))
        best = pari.idealred(K.nf, I)
        assert pari.bnfisprincipal(K.bnf, best, 0) == cl
    K._crep[key] = best
    return best

def sqred(K, z, known=()):
    """z = zz s^2 with zz integral, (zz) = J K0^2 (J square-free, K0 small), zz unit-balanced (real field)."""
    zl = pari.lift(z)
    den = pari.denominator(pari.content(zl))
    zz = K.el(zl * den * den); s = K.el(1) / den
    I = pari.idealhnf(K.nf, 1); J = []
    for P, e in ideal_fact(K, zz, known):
        if e >= 2: I = pari.idealmul(K.nf, I, pari.idealpow(K.nf, P, e // 2))
        if e % 2: J.append(P)
    cl = pari.bnfisprincipal(K.bnf, I, 0)
    if any(int(c) for c in cl): I = pari.idealdiv(K.nf, I, class_rep(K, cl))
    ip = pari.bnfisprincipal(K.bnf, I, 3)
    assert not any(int(c) for c in ip[0])
    al = K.el(pari.nfbasistoalg(K.nf, ip[1]))
    zz = zz / al ** 2; s = s * al
    if K.real:
        eps = K.el(pari.nfbasistoalg(K.nf, K.bnf.bnf_get_fu()[0]))
        f1, f2 = [abs(x) for x in K.emb(eps)]
        lf = float(mpmath.log(f1 / f2))
        for _ in range(200):
            e1, e2 = [abs(x) for x in K.emb(zz)]
            k_ = round(float(mpmath.log(e1 / e2)) / 2 / lf)
            if k_ == 0: break
            zz = zz / eps ** (2 * k_); s = s * eps ** k_
    assert zz * s ** 2 == z
    Jn = 1
    for P in J: Jn *= int(pari.idealnorm(K.nf, P))
    return zz, s, J, Jn

# ---------------------------------------------------------------- the descent
def meet(K, Lat, phi, P):
    H = pari.idealhnf(K.nf, P)
    Hinv = pari.matsolve(H, pari.matid(2))
    R_ = Hinv * pari.matrix(2, 4, phi[0] + phi[1]) * Lat
    d_ = pari.denominator(R_)
    K_ = pari.matkerint(pari.matconcat([R_ * d_, d_ * pari.matid(2)]))
    nc = int(pari.matsize(K_)[1])
    return Lat * pari.mathnf(pari.matrix(4, nc, [K_[i, j] for i in range(4) for j in range(nc)]))

def mat_of(K, el):
    c0 = K.coords(K.el(el)); c1 = K.coords(K.el(el) * K.w)
    return [[c0[0], c1[0]], [c0[1], c1[1]]]

def lattice_for(K, J, B):
    """{(u, v) in O_L^2 : u = r v mod P (r^2 = B mod P), or u in P where P | B, for the odd P of J}."""
    Lat = pari.matid(4)
    for P in J:
        if int(P[0]) == 2: continue
        vB = int(pari.idealval(K.nf, pari.lift(B), P))
        if vB % 2 == 1: phi = [[1, 0, 0, 0], [0, 1, 0, 0]]
        elif vB == 0:
            modpr = pari.nfmodprinit(K.nf, P)
            b = pari.nfmodpr(K.nf, pari.lift(B), modpr)
            assert pari.issquare(b), 'B is not a square modulo P: the conic is not soluble at P'
            r = K.el(pari.nfbasistoalg(K.nf, pari.nfmodprlift(K.nf, pari.sqrt(b), modpr)))
            Mr = mat_of(K, -r)
            phi = [[1, 0, Mr[0][0], Mr[0][1]], [0, 1, Mr[1][0], Mr[1][1]]]
        else: continue
        Lat = meet(K, Lat, phi, P)
    return Lat

def emb_rows(K, wts, nblk):
    """rows of a real matrix M with |M x|^2 = sum_sigma sum_blk wts[blk][sigma] |sigma(x_blk)|^2 (x_blk in O_L)."""
    w0 = K.emb(K.el(1)); ww = K.emb(K.w)
    rows = []
    for s_ in range(len(w0)):
        for blk in range(nblk):
            sw = mpmath.sqrt(wts[blk][s_])
            if K.real:
                r = [mpmath.mpf(0)] * (2 * nblk); r[2 * blk] = sw * w0[s_]; r[2 * blk + 1] = sw * ww[s_]; rows.append(r)
            else:
                rr = [mpmath.mpf(0)] * (2 * nblk); ri = [mpmath.mpf(0)] * (2 * nblk)
                for i, b_ in enumerate((w0[s_], ww[s_])):
                    rr[2 * blk + i] = sw * mpmath.sqrt(2) * mpmath.re(b_); ri[2 * blk + i] = sw * mpmath.sqrt(2) * mpmath.im(b_)
                rows.append(rr); rows.append(ri)
    return rows

def short_vectors(K, Lat, rows, nvec, Rmul=4):
    """short vectors of the integer lattice Lat (columns) for |rows x|^2: LLL on scaled integers, qfminim."""
    n = int(pari.matsize(Lat)[0]); k = int(pari.matsize(Lat)[1])
    with mpmath.workdps(150):
        Lm = [[int(Lat[i, j]) for j in range(k)] for i in range(n)]
        ML = [[sum(rows[r][a] * Lm[a][j] for a in range(n)) for j in range(k)] for r in range(len(rows))]
        G0 = [[sum(ML[r][i] * ML[r][j] for r in range(len(rows))) for j in range(k)] for i in range(k)]
        sc = mpmath.mpf(10) ** 40 / abs(mpmath.det(mpmath.matrix(G0))) ** (mpmath.mpf(1) / (2 * k))
        Mi = pari.matrix(len(rows), k, [int(mpmath.nint(ML[r][j] * sc)) for r in range(len(rows)) for j in range(k)])
        U = pari.qflll(Mi)
        Uu = [[int(U[i, j]) for j in range(k)] for i in range(k)]
        MU = [[sum(ML[r][a] * Uu[a][j] for a in range(k)) for j in range(k)] for r in range(len(rows))]
        G = [[sum(MU[r][i] * MU[r][j] for r in range(len(rows))) for j in range(k)] for i in range(k)]
        m0 = G[0][0]
        G2 = pari.matrix(k, k, [int(mpmath.nint(G[i][j] / m0 * 10 ** 15)) for i in range(k) for j in range(k)])
    vs = pari.qfminim(G2, int(G2[0, 0]) * Rmul, nvec, 2)[2]
    LU = Lat * U
    out = []
    for j in range(int(pari.matsize(vs)[1])):
        x = LU * pari.matrix(k, 1, [vs[i, j] for i in range(k)])
        out.append([int(x[i, 0]) for i in range(n)])
    return out

def descend(K, A, B, rounds=40, nvec=1500, target=10 ** 11, known=(), seed=1, trials=6):
    """random-walk descent: each round tries both sides and takes the step with the smallest N(J_A) N(J_B) (a
    non-decreasing step with probability 1/2); returns the best state and the steps leading to it."""
    rng = random.Random(seed)
    A, sA, JA, nA = sqred(K, A, known); B, sB, JB, nB = sqred(K, B, known)
    steps = [dict(side='init', sA=sA, sB=sB)]
    best = (nA * nB, 1, A, B, nA, nB)
    for rd in range(rounds):
        if nA * nB <= target: break
        opts = []
        for side in ('A', 'B'):
            X, JX, nX, Y = (A, JA, nA, B) if side == 'A' else (B, JB, nB, A)
            Lat = lattice_for(K, JX, Y)
            eY = [abs(x) for x in K.emb(Y)]
            cands = []
            for trial in range(trials):
                mu = [mpmath.mpf(1), mpmath.exp(rng.uniform(-4, 4))] if K.real else [mpmath.mpf(1)]
                tB = mpmath.exp(rng.uniform(-2, 2)) if trial else mpmath.mpf(1)
                rows = emb_rows(K, [[mu[s_] for s_ in range(len(mu))], [mu[s_] * tB * eY[s_] for s_ in range(len(mu))]], 2)
                for vec in short_vectors(K, Lat, rows, nvec):
                    u = K.from_coords(vec[0:2]); v = K.from_coords(vec[2:4])
                    n = u * u - Y * v * v
                    if n == 0: continue
                    Nn = abs(int(pari.norm(n)))
                    sc = Nn // nX if Nn % nX == 0 else Nn
                    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47):
                        while sc % (p * p) == 0: sc //= p * p
                    cands.append((sc, rng.random(), u, v, n))
            cands.sort(key=lambda c_: (c_[0], c_[1]))
            for sc, _, u, v, n in cands[:3]:
                Xn, s, Jn, nn = sqred(K, n * X, known)
                opts.append((nn * (nB if side == 'A' else nA), side, u, v, n, Xn, s, Jn, nn))
        opts.sort(key=lambda o: o[0])
        o = opts[0]
        if o[0] >= nA * nB and rng.random() < 0.5 and len(opts) > 1: o = opts[rng.randrange(min(4, len(opts)))]
        prod, side, u, v, n, Xn, s, Jn, nn = o
        steps.append(dict(side=side, u=u, v=v, n=n, s=s, A=A, B=B))
        if side == 'A': A, JA, nA = Xn, Jn, nn
        else: B, JB, nB = Xn, Jn, nn
        if nA * nB < best[0]: best = (nA * nB, len(steps), A, B, nA, nB)
    _, ln, A, B, nA, nB = best
    return A, B, steps[:ln], nA, nB

def reduce_point(K, P):
    """a projective point divided by a generator of (a small ideal times) its gcd ideal, unit-balanced."""
    X, Y, Z = [K.el(c) for c in P]
    def clear(X, Y, Z):
        den = pari.lcm(pari.lcm(pari.denominator(pari.content(pari.lift(X))), pari.denominator(pari.content(pari.lift(Y)))),
                       pari.denominator(pari.content(pari.lift(Z))))
        return X * den, Y * den, Z * den
    X, Y, Z = clear(X, Y, Z)
    I = pari.idealadd(K.nf, pari.idealadd(K.nf, pari.lift(X), pari.lift(Y)), pari.lift(Z))
    cl = pari.bnfisprincipal(K.bnf, I, 0)
    if any(int(c) for c in cl): I = pari.idealdiv(K.nf, I, pari.idealred(K.nf, I))
    g = K.el(pari.nfbasistoalg(K.nf, pari.bnfisprincipal(K.bnf, I, 3)[1]))
    X, Y, Z = clear(X / g, Y / g, Z / g)
    if K.real:
        eps = K.el(pari.nfbasistoalg(K.nf, K.bnf.bnf_get_fu()[0]))
        f1, f2 = [abs(x) for x in K.emb(eps)]
        lf = float(mpmath.log(f1 / f2))
        for _ in range(100):
            m = [max(abs(x) for x in ys) for ys in zip(*[K.emb(c) for c in (X, Y, Z) if c != 0])]
            k_ = round(float(mpmath.log(m[0] / m[1])) / 2 / lf)
            if k_ == 0: break
            X, Y, Z = X / eps ** k_, Y / eps ** k_, Z / eps ** k_
    return X, Y, Z

def digits(P):
    return max(len(str(pari.lift(c))) for c in P)

def unwind(K, steps, sol):
    """the point of the conic before the descent steps (the 'init' step excluded), reduced after each step."""
    X, Y, Z = reduce_point(K, sol)
    for st in reversed(steps):
        if st['side'] == 'A':
            u, v, n, s, B = st['u'], st['v'], st['n'], st['s'], st['B']
            X, Y, Z = (u * X + B * v * Z) * s, n * Y, (u * Z + v * X) * s
        elif st['side'] == 'B':
            u, v, n, s, A = st['u'], st['v'], st['n'], st['s'], st['A']
            X, Y, Z = (u * X + A * v * Y) * s, (u * Y + v * X) * s, n * Z
        else: continue
        X, Y, Z = reduce_point(K, (X, Y, Z))
    return X, Y, Z

def qmul(A, B, p, q):
    """product in the quaternion algebra (A, B): i^2 = A, j^2 = B, k = i j = -j i."""
    a1, b1, c1, d1 = p; a2, b2, c2, d2 = q
    return (a1 * a2 + A * b1 * b2 + B * c1 * c2 - A * B * d1 * d2, a1 * b2 + b1 * a2 - B * c1 * d2 + B * d1 * c2,
            a1 * c2 + c1 * a2 + A * b1 * d2 - A * d1 * b2, a1 * d2 + d1 * a2 + b1 * c2 - c1 * b2)

def small_points(K, A, B, P0, nvec=200):
    """small points of X^2 - A Y^2 - B Z^2 = 0 from one point P0 (see the docstring, (5))."""
    X0, Y0, Z0 = [K.el(c) for c in P0]
    assert X0 ** 2 - A * Y0 ** 2 - B * Z0 ** 2 == 0
    xib = (X0, -Y0, -Z0, K.el(0))
    cols = []
    for m in range(4):
        for b_ in (K.el(1), K.w):
            eta = [K.el(0)] * 4; eta[m] = b_
            col = []
            for c in qmul(A, B, eta, xib):
                cc = pari.nfalgtobasis(K.nf, pari.lift(c)); col += [cc[0], cc[1]]
            cols.append(col)
    Mq = pari.matrix(8, 8, [cols[j][i] for i in range(8) for j in range(8)])
    ker = pari.matkerint(Mq * pari.denominator(Mq))          # (H xi) n O_L<1, i, j, k>, rank 4
    assert int(pari.matsize(ker)[1]) == 4
    wts = [[abs(x) for x in K.emb(c)] for c in (K.el(1), A, B, A * B)]
    pts = []
    for vec in short_vectors(K, ker, emb_rows(K, wts, 4), nvec, Rmul=6):
        x, y, z, w = [K.from_coords(vec[2 * m:2 * m + 2]) for m in range(4)]
        x, y, z, w = x, -y, -z, -w                             # conj(v)
        P = (z * x - A * w * y, z * y - w * x, z * z - A * w * w)
        if all(c == 0 for c in P): continue
        assert P[0] ** 2 - A * P[1] ** 2 - B * P[2] ** 2 == 0
        pts.append(reduce_point(K, P))
    pts.sort(key=digits)
    return pts

# ---------------------------------------------------------------- the conics of C3'
def curve():
    js = json.load(open(os.path.join(here, 'rank_fisher_curves.json')))
    cv = next(c for c in js if c['name'].split(' ')[0] == NAME)
    return cv, cv['sides'][SIDE]

def sqfree_split(n):
    fa = pari.factor(abs(n)); D0 = -1 if n < 0 else 1; m = 1
    for p_, ex in zip(fa[0], fa[1]):
        D0 *= int(p_) ** (int(ex) % 2); m *= int(p_) ** (int(ex) // 2)
    assert m * m * D0 == n
    return D0, m

def conic(PAIR, xi):
    """A1, B1, b20, b11, h0 of the U (or V) conic of C3' with level-3 twist xi, as PARI elements of Q(sqrt D0)."""
    cv, S = curve()
    R5.HINTS.clear()
    R = R5.analyse(cv['name'], SIDE, dict(S['level3_curve'], xi=str(xi)), dd=int(S['dd']))
    if PAIR == 'U': h0, Bq, v1, v2, dd = R['h0'], R['Bq'], R5.e, R5.g, R['dU']
    else: h0, Bq, v1, v2, dd = R['coneV'][3], R['coneV'][4], R5.c, R5.d, R['dV']
    D0, m = sqfree_split(dd)
    def topari(expr):
        P_ = sp.Poly(sp.expand(expr), t_)
        return pari('Mod((%s) + (%s)*%d*y, y^2 - (%d))' % (P_.coeff_monomial(1), P_.coeff_monomial(t_), m, D0))
    P = sp.Poly(Bq, v1, v2)
    b20, b11, b02 = [topari(P.coeff_monomial(mm)) for mm in (v1 ** 2, v1 * v2, v2 ** 2)]
    h0p = topari(h0)
    return dict(cv=cv, D0=D0, m=m, b20=b20, b11=b11, b02=b02, h0=h0p, A1=-4 * b20 * h0p, B1=b11 ** 2 - 4 * b20 * b02)

def find_hint(PAIR, xi, target=10 ** 11, log=print):
    from cysignals.alarm import alarm, AlarmInterrupt, cancel_alarm
    t0 = time.time()
    Cn = conic(PAIR, xi)
    K = Fld(Cn['D0'])
    A1, B1 = Cn['A1'], Cn['B1']
    known = set()
    for z in (A1, B1):
        zl = pari.lift(z); den = pari.denominator(pari.content(zl))
        known |= set(factor_int(pari.norm(pari.Mod(zl * den * den, K.pol))))
    known = tuple(sorted(known))
    sol = None
    for seed in range(1, 12):
        A, B, steps, nA, nB = descend(K, A1, B1, target=target, known=known, seed=seed)
        log('  %s, descent seed %d: N(J_A), N(J_B) of %d, %d digits after %d steps' % (PAIR, seed, len(str(nA)), len(str(nB)), len(steps) - 1))
        if nA * nB > target: continue
        swap = nB < nA
        al, be = (B, A) if swap else (A, B)
        alarm(60)
        try:
            r = pari.rnfisnorm(pari.rnfisnorminit(K.pol, pari('x^2') - pari.lift(al), 2), pari.lift(be), 0)
        except AlarmInterrupt:
            log('  %s: rnfisnorm over 60 s, next seed' % PAIR); continue
        finally:
            cancel_alarm()
        if r[1] != 1: continue
        a = pari.lift(r[0])
        a0 = K.el(pari.polcoef(a, 0, 'x')); a1 = K.el(pari.polcoef(a, 1, 'x'))
        P0 = (a0, K.el(1), a1) if swap else (a0, a1, K.el(1))
        if digits(P0) > 20000: continue
        sol = P0; break
    assert sol is not None, 'no solution of the norm equation'
    assert sol[0] ** 2 - A * sol[1] ** 2 - B * sol[2] ** 2 == 0
    X, Y, Z = unwind(K, steps[1:], small_points(K, A, B, sol)[0])
    sA, sB = steps[0]['sA'], steps[0]['sB']
    X, Y, Z = small_points(K, A1 / sA ** 2, B1 / sB ** 2, (X, Y, Z))[0]
    X, Y, Z = reduce_point(K, (X, Y / sA, Z / sB))
    assert X ** 2 - A1 * Y ** 2 - B1 * Z ** 2 == 0
    ee = (X - Cn['b11'] * Z) / (2 * Cn['b20'])
    assert Cn['h0'] * Y ** 2 + Cn['b20'] * ee ** 2 + Cn['b11'] * ee * Z + Cn['b02'] * Z ** 2 == 0
    out = lambda z: [str(pari.polcoef(pari.lift(z), 0, 'y')), str(pari.polcoef(pari.lift(z), 1, 'y'))]
    key = '%s|%s' % (Cn['cv']['name'], SIDE)
    if PAIR == 'U': return key, dict(D0=Cn['D0'], m=Cn['m'], Z=out(Y), e=out(ee), g=out(Z)), time.time() - t0
    return key + '|V', dict(D0=Cn['D0'], m=Cn['m'], Z=out(Y), c=out(ee), d=out(Z)), time.time() - t0

# ---------------------------------------------------------------- modes
def committed():
    H = json.load(open(HINTS_GEN))
    cv, S = curve()
    key = '%s|%s' % (cv['name'], SIDE)
    return cv, S, key, {k: H[k] for k in (key, key + '|V')}, int(H[key + '|xi'])

def recheck():
    cv, S, key, hints, xi = committed()
    C0 = S['level3_curve']
    say('rank_fisher_theta4_norm.py: the committed conic points of %s, %s side (rank_fisher_theta4_gen_hints.json), '
        're-checked exactly' % (cv['name'], SIDE))
    say('C3 (rank_fisher_curves.json): lam = %s, xi = %s; used with xi = %d (the hint \'%s|xi\'); xi / %s = %d mod squares'
        % (C0['lam'], C0['xi'], xi, key, C0['xi'], R5.core_small(xi * int(C0['xi']))))
    for k_ in (key, key + '|V'):
        h = hints[k_]
        K = Fld(int(h['D0']))
        say('  %s: over Q(sqrt %d) (class group %s under GRH, PARI bnfinit: for information), m = %s; digits of the '
            'coordinates %s' % (k_, int(h['D0']), K.bnf.bnf_get_cyc(), h['m'],
                                [len(s_) for kk in ('Z', 'e', 'g', 'c', 'd') if kk in h for s_ in h[kk]]))
    R5.HINTS.clear(); R5.HINTS.update(hints)
    R = R5.analyse(cv['name'], SIDE, dict(C0, xi=str(xi)), dd=int(S['dd']))
    assert R['pushout'] is not None and R['pushout3'] is not None and R['comp'] is not None
    say('  rank_fisher5.analyse, steps (1)-(8), exact: P1 on the cone over Gamma_1, l_1 tangent along the ruling, '
        'f = N(l_1) rational; P3 on the cone over Gamma_3, l_3 tangent, f_dag rational; q through the 8 points and '
        'f f_dag = c q^2 mod I(C3\') with c of %d / %d digits' % (len(str(abs(R['comp']['c'].p))), len(str(R['comp']['c'].q))))
    n = ok = 0
    for what, C_, hs in (('the same points on the json\'s C3 (xi = %s): another conic' % C0['xi'], C0, hints),
                         ('P1 moved off its cone (Z -> Z + 1)', dict(C0, xi=str(xi)),
                          dict(hints, **{key: dict(hints[key], Z=[str(Fraction(hints[key]['Z'][0]) + 1), hints[key]['Z'][1]])}))):
        n += 1
        R5.HINTS.clear(); R5.HINTS.update(hs)
        try:
            R5.analyse(cv['name'], SIDE, C_, dd=int(S['dd'])); say('  NOT REFUSED: ' + what)
        except AssertionError as ex:
            ok += 1; say('  refused: %s (%s)' % (what, str(ex)[:60]))
    say('refusals: %d of %d' % (ok, n))
    assert ok == n

def find():
    cv, S, key, hints, xi = committed()
    out = {}
    for PAIR in ('U', 'V'):
        k_, h, dt = find_hint(PAIR, xi)
        print('%s: found in %.1f s, digits %s; equal to the committed hint: %s'
              % (k_, dt, [len(s_) for kk in ('Z', 'e', 'g', 'c', 'd') if kk in h for s_ in h[kk]], h == hints[k_]), flush=True)
        out[k_] = h
    if WRITE:
        H = json.load(open(HINTS_GEN))
        H.update(out); H[key + '|xi'] = str(xi)
        json.dump(H, open(HINTS_GEN, 'w'), indent=1, sort_keys=True)
        print('written to', HINTS_GEN)

def obstruction():
    import rank_fisher_theta4 as T
    import rank_indep as RI, rank_second as R2
    cv, S, key, hints, xi = committed()
    C0 = S['level3_curve']
    m_, n_ = [int(x) for x in NAME.split(',')]
    a, b, G = R2.sigma_curve(m_, n_, 'E+')
    RD = RI.first_descent(a, b)
    places = [2, 3, 5, 7, 11, 13, 17, 31, 41, 61, 97, 157, 349, 977, 1193, 245411, 1493911429, 100067394997]
    say('rank_fisher_theta4_norm.py --obstruction: %s, %s side, lam = %s: the level-3 twists xi = %s s, s in S_1 / S_2 '
        '(a sampler: finder points, sampled local classes of f0)' % (cv['name'], SIDE, C0['lam'], C0['xi']))
    say('  S_1 of this side (rank_indep): %s; S_2 = <-1, -5820647> (rank_fisher.txt); places: R and %s (B0 and 11)'
        % (RD['Sp'], places))
    for s in (1, 14, 17, 238):
        xi_ = R5.core_small(int(C0['xi']) * s)
        k_, h, dt = find_hint('U', xi_, log=lambda s_: None)
        R5.HINTS.clear(); R5.HINTS[k_] = h
        C = dict(C0, xi=str(xi_))
        R = R5.analyse(cv['name'], SIDE, C, dd=int(S['dd']))
        s0, f0 = T.poly_int(sp.expand(R['pushout']['f']), T.VS)
        D = dict(QA=sp.expand(R['QA']), QB=sp.expand(R['QB']), f0=f0)
        M = T.Model(D['QA'], D['QB'])
        cls = {}
        for p in places:
            Kp = {}
            for sd in range(1, 9):
                Kp.update(T.f0_classes(M, D, C, p, n=150, seed=sd))
                if Kp and sd >= 2: break
            assert Kp, 'no sampled point at %d' % p
            cls[p] = sorted(Kp)
        Rs = set(T.real_signs(M, D, C=C)); assert len(Rs) == 1
        ob = {}
        for eta in (-1, 14, 17, 97826, -5820647):
            prod = -1 if (Rs == {-1} and eta < 0) else 1
            for p in places:
                hs = {T.hilbert(k, T.int_class(eta, p), p) for k in cls[p]}
                assert len(hs) == 1, 'the symbol varies on the sampled classes'
                prod *= hs.pop()
            ob[eta] = prod
        say('  xi = %d (s = %d): sign of f0 on C3\'(R) %s; Ob(eta) = prod_v (f0(P_v), eta)_v: %s -> %s'
            % (xi_, s, Rs.pop(), ' '.join('%d:%+d' % (e_, v_) for e_, v_ in ob.items()),
               'level-4 Selmer set non-empty' if all(v_ == 1 for v_ in ob.values()) else 'empty'))

if __name__ == '__main__':
    if '--find' in sys.argv: find()
    elif '--obstruction' in sys.argv:
        obstruction()
        if WRITE: open(os.path.join(here, 'rank_fisher_theta4_norm_obstruction.txt'), 'w').write('\n'.join(OUT) + '\n')
    else:
        recheck()
        if WRITE: open(os.path.join(here, 'rank_fisher_theta4_norm.txt'), 'w').write('\n'.join(OUT) + '\n')
