#!/usr/bin/env python3
"""kctp_recheck.py: an independent re-check of the inputs of kctp.py (run 3, branch r3_v7) for the corner K-curves of
(627,364) over K = Q(i), with new code.  EXPERIMENT; written by an AI model (Claude), not read by a person.

    python3 kctp_recheck.py              the target curve L1L3M1M2 of (627,364)       -> kctp_recheck.txt
    python3 kctp_recheck.py --controls   base changes from Q with known ranks          -> kctp_recheck_controls.txt
    python3 kctp_recheck.py --negative   forged inputs that must be refused            -> kctp_recheck_negative.txt
    (add --write to write the output file)

Nothing is imported from kctp.py or kdescent.py.  Sage (passagemath) is used for exact arithmetic in K, for prime
ideals, residue fields and discrete logarithms modulo (1 + i)^5 (PARI's ideallog), and to find points on conics
(each point is checked exactly here).  Floating point is used nowhere in a certifying step.

WHAT IS CHECKED.
 1. The curve.  Q_J = L1 L3 M1 M2 is built from the linear factors of kclose.py (its docstring:
    L1 = m - n s, L3 = n + m s, M1 = (m+in) s - (m-in), M2 = (m-in) s + (m+in)).  Its Jacobian is computed from the
    classical invariants I, J of the quartic (y^2 = x^3 - 27 I x - 27 J).  The model with full 2-torsion is made by
    sending the root of M2 to infinity.  Both are checked to be isomorphic over K (exactly, Sage is_isomorphic) to the
    curve y^2 = (x - e1)(x - e2)(x - e3) of kctp.txt (its e_k copied here), and j is checked not to be rational.
 2. Torsion.  Points of order 2 and 4 are listed and checked with the group law; #E(F_P) at small primes of good
    reduction (Sage's point count, and a naive count at two primes) has gcd 8, so E(K)_tors = Z/2 x Z/4 exactly.
 3. Sel^2(E/K) by the complete 2-descent, new code: K(S,2) with basis i and generators of the primes of S (h_K = 1);
    at each place v of S the image of E(K_v)/2E(K_v) in (K_v*/K_v*^2)^2 is found from local points (x in K with f(x) a
    square in K_v) until it has its known dimension (2 at odd v, 4 at 1 + i); local square classes from the valuation
    and the quadratic character of the residue field (odd v) or the discrete logarithm modulo (1 + i)^5 (v = 1 + i;
    U^(5) lies in the squares).  Sel^2 is the kernel of the linear map to the product of the quotients.  The four
    Selmer elements printed by kctp.txt are checked to span the group found here.
 4. The lift criterion (Cassels' second descent; no Hilbert symbols).  For a in Sel^2, D_a: x = e1 + a1 u^2 =
    e2 + a2 v^2 = e3 + a3 w^2 in P^3 (u : v : w : z).  With points on the conics C_w: a1 u^2 - a2 v^2 + (e1 - e2) z^2,
    C_v: a1 u^2 - a3 w^2 + (e1 - e3) z^2, C_u: a2 v^2 - a3 w^2 + (e2 - e3) z^2 and the primitive tangent lines l_u, l_v,
    l_w there, g1 = l_u l_w / z^2 and g2 = l_v l_w / z^2 have divisors 2 (degree-0 divisor), so the 4-coverings over D_a
    are D_a  cap  {g1 = lam1 square, g2 = lam2 square}.  D_a lifts to an everywhere locally soluble 4-covering iff there
    is lam in K(S,2)^2 with lam in g(D_a(K_v)) mod squares at every v in S (outside S, g has even valuation on
    D_a(K_v): D_a has good reduction and the lines are primitive -- the purity argument, quoted).  g(D_a(K_v)) is found
    by sampling local points of D_a (all 8 sign choices of (u, v, w)); it is checked to be an affine space of the size
    of E(K_v)/2E(K_v).  Quoted theory (Cassels): a class of Sel^2 lifts to Sel^4 iff it is in the kernel of the
    Cassels-Tate pairing on Sel^2.  So "every basis class lifts" re-checks kctp.py's zero pairing.

WHAT IT GIVES.  Only a re-check: dim Sel^2 = 4, torsion image dim 2, every class lifts (pairing zero), hence rank
E(K) <= 2; nothing about 0 versus 2.
"""
from sage.all__sagemath_schemes import *
import sys, os, random, time, itertools          # after Sage, which exports its own random and time

K = QuadraticField(-1, 'i'); I = K.gen()
OK_ = K.ring_of_integers()
WRITE = '--write' in sys.argv
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(str(s))

# ------------------------------------------------------------------ the curve
E_KCTP = {'L1L3M1M2': [(-2783595387469704535000, 28543239158509955880000), (72914933236647971246250, 28543239158509955880000),
                       (-31221116673647194769375, 61411774313468584005000)]}

def quartic(m, n, J):
    R = PolynomialRing(K, 's'); s = R.gen(); p = m + n * I; pb = m - n * I
    F = {'L1': m - n * s, 'L2': m + n * s, 'L3': n + m * s, 'L4': n - m * s,
         'M1': p * s - pb, 'M2': pb * s + p, 'M3': p * s + pb, 'M4': pb * s - p}
    fs = [F[J[k:k + 2]] for k in range(0, 8, 2)]
    return prod(fs), fs

def jac_invariants(Q):
    e_, d_, c_, b_, a_ = [Q[k] for k in range(5)]
    Iv = 12 * a_ * e_ - 3 * b_ * d_ + c_ ** 2
    Jv = 72 * a_ * c_ * e_ + 9 * b_ * c_ * d_ - 27 * a_ * d_ ** 2 - 27 * e_ * b_ ** 2 - 2 * c_ ** 3
    return EllipticCurve(K, [-27 * Iv, -27 * Jv])

def full2_model(Q, fs):
    """send the root r of the last factor to infinity: s = r + 1/X, Y = y X^2 gives Y^2 = A X^3 + B X^2 + C X + D with
    Q(r + T) = A T + B T^2 + C T^3 + D T^4; x = A X gives x^3 + B x^2 + A C x + A^2 D, roots x_k = A / (r_k - r)."""
    roots = [-f[0] / f[1] for f in fs]; r = roots[3]
    R = Q.parent(); T = R.gen()
    Qr = Q(r + T); assert Qr[0] == 0
    A = Qr[1]
    xs = [A / (rk - r) for rk in roots[:3]]
    den = lcm([x.denominator() for x in xs])           # x -> den^2 x keeps the model y^2 = prod (x - e_k)
    es = [x * den ** 2 for x in xs]
    assert all(ee in OK_ for ee in es)
    cub = prod(T - ee for ee in es)
    gcub = T ** 3 + Qr[2] * T ** 2 + A * Qr[3] * T + A ** 2 * Qr[4]
    assert cub == gcub(T / den ** 2) * den ** 6
    return es

def curve_from_e(es):
    R = PolynomialRing(K, 'x'); x = R.gen(); f = prod(x - ee for ee in es)
    return EllipticCurve(K, [0, f[2], 0, f[1], f[0]])

# ------------------------------------------------------------------ local square classes
class Place:
    def __init__(self, P):
        self.P = P; self.pi = P.gens_reduced()[0]; self.q = P.norm()
        self.two = (self.q == 2)
        if self.two:
            self.M = P ** 5
            self.cyc = list(self.M.idealstar(2).gens_orders())
            assert sum(1 for c in self.cyc if c % 2 == 0) == 3     # K_v*/K_v*^2 has dimension 4 at 1 + i
            self.w = 4
        else:
            self.F = K.residue_field(P); self.w = 2
    def val(self, x): return x.valuation(self.P)
    def cls(self, x):
        """F_2-vector of x in K_v*/K_v*^2: (v(x) mod 2, unit class)."""
        assert x != 0
        v = self.val(x); u = x / self.pi ** v
        if self.two:
            lg = self.M.ideallog(u)
            return [v % 2] + [int(l) % 2 for l, c in zip(lg, self.cyc) if c % 2 == 0]
        return [v % 2, 0 if self.F(u).is_square() else 1]
    def is_square(self, x):
        return x != 0 and not any(self.cls(x))

def f2_rank(rows):
    if not rows: return 0
    return matrix(GF(2), rows).rank()

def kummer(x, es):
    e1, e2, e3 = es
    if x == e1: return ((e1 - e2) * (e1 - e3), e1 - e2)
    if x == e2: return (e2 - e1, (e2 - e1) * (e2 - e3))
    return (x - e1, x - e2)

def sample_x(pl, es, rng, k):
    """candidate x-coordinates, integral at v: random, near the roots, near 0."""
    t = rng.randint(0, 5)
    z = K(rng.randint(-60, 60) + rng.randint(-60, 60) * I)
    if t <= 1: return z
    j = rng.randint(1, 3 + k // 50)
    if t <= 4: return es[rng.randint(0, 2)] + pl.pi ** j * z
    return pl.pi ** j * z

def local_image(pl, es, rng, maxit=40000):
    """the image of E(K_v)/2E(K_v) in (K_v*/K_v*^2)^2, and local points (x, class) found on the way."""
    e1, e2, e3 = es
    target = 2 + (2 if pl.two else 0)
    rows = []; pts = []
    for x in es:
        if x != e3:
            d1, d2 = kummer(x, es); vec = pl.cls(d1) + pl.cls(d2); rows.append(vec); pts.append((x, tuple(vec)))
    d1, d2 = (e3 - e1, e3 - e2); rows.append(pl.cls(d1) + pl.cls(d2))
    k = 0
    while k < maxit:
        k += 1
        x = sample_x(pl, es, rng, k)
        if x in es: continue
        fx = (x - e1) * (x - e2) * (x - e3)
        if not pl.is_square(fx): continue
        vec = pl.cls(x - e1) + pl.cls(x - e2)
        pts.append((x, tuple(vec)))
        if f2_rank(rows + [vec]) > f2_rank(rows): rows.append(vec)
        if f2_rank(rows) == target and len(pts) > 400: break
    assert f2_rank(rows) == target, 'local image too small at %s' % pl.P
    B = matrix(GF(2), rows).row_space().basis()
    return [list(b) for b in B], pts

# ------------------------------------------------------------------ Sel^2
class Selmer:
    def __init__(self, es, rng, S=None):
        self.es = es
        D = 2 * (es[0] - es[1]) * (es[0] - es[2]) * (es[1] - es[2])
        Ps = S if S is not None else sorted(set(P for P, _ in K.ideal(D).factor()), key=lambda P: (P.norm(), str(P.gens_reduced())))
        self.places = [Place(P) for P in Ps]
        self.gens = [K(I)] + [pl.pi for pl in self.places]           # K(S,2), h_K = 1
        self.n = len(self.gens)
        self.img = {}; self.pts = {}
        for pl in self.places:
            self.img[pl.P], self.pts[pl.P] = local_image(pl, es, rng)
        rows = []
        for pl in self.places:
            # functionals vanishing on the local image
            Wm = matrix(GF(2), self.img[pl.P]); ann = Wm.right_kernel().basis()
            res = [pl.cls(g) for g in self.gens]                 # n x w
            for h in ann:
                h = list(h); row = []
                for c in range(2):
                    for j in range(self.n):
                        row.append(sum(h[c * pl.w + t] * res[j][t] for t in range(pl.w)) % 2)
                rows.append(row)
        Mt = matrix(GF(2), rows)
        self.basis = [list(b) for b in Mt.right_kernel().basis()]
        self.dim = len(self.basis)
    def elt(self, bits):
        return prod(g ** int(b) for g, b in zip(self.gens, bits))
    def pair(self, vec):
        b1 = vec[:self.n]; b2 = vec[self.n:]
        return (self.elt(b1), self.elt(b2), self.elt([(x + y) % 2 for x, y in zip(b1, b2)]))
    def vec_of(self, d):
        """S-unit vector of d mod squares (basis i, pi_P); None if d is not an S-unit times a square."""
        rest = K(d); v = []
        for pl in self.places:
            k = pl.val(rest); v.append(k % 2); rest = rest / pl.pi ** k
        fac = K.ideal(rest).factor()
        if any(ee % 2 for P, ee in fac): return None
        h = prod([P.gens_reduced()[0] ** (ee // 2) for P, ee in fac]) if len(fac) else K(1)
        u = rest / h ** 2
        assert u in OK_ and u.norm() == 1 and u in (K(1), K(-1), I, -I)
        return [0 if u in (K(1), K(-1)) else 1] + v
    def contains(self, d1, d2):
        """local test of (d1, d2) at every place of S (and S-unit test)."""
        if self.vec_of(d1) is None or self.vec_of(d2) is None: return False
        for pl in self.places:
            vec = pl.cls(d1) + pl.cls(d2)
            if f2_rank(self.img[pl.P] + [vec]) > len(self.img[pl.P]): return False
        return True

# ------------------------------------------------------------------ local square roots and the lift criterion
def to_int_mod(x, pl, N):
    """an integral y with y = x mod P^N (x integral at P)."""
    assert pl.val(x) >= 0
    d = K(x.denominator()); y = x * d                   # y, d in Z[i]
    k = pl.val(d); y = y / pl.pi ** k; d = d / pl.pi ** k   # still in Z[i]; d a unit at P
    assert y in OK_ and d in OK_ and pl.val(d) == 0
    PN = pl.P ** N
    return PN.reduce(y * d.inverse_mod(PN))

def sqrt_local(y, pl, N):
    """s in K with v(s^2 - y) >= N + v(y); y a square in K_v (v(y) even, any sign)."""
    v = pl.val(y); assert v % 2 == 0
    u = to_int_mod(y / pl.pi ** v, pl, N + 8)
    PN = pl.P ** (N + 8)
    if pl.two:
        s = None
        for a, b in itertools.product(range(8), range(4)):
            c = K(a + b * I)
            if (c * c - u).valuation(pl.P) >= 7: s = c; break
        assert s is not None, 'not a square at 1 + i'
        for _ in range(12):
            t = (s * s + u) / (2 * s)                  # Newton; v(s^2 - u) >= 7 makes it converge
            s = to_int_mod(t, pl, N + 8)
    else:
        s = pl.F.lift(pl.F(u).sqrt()); inv2 = K(2).inverse_mod(PN)
        for _ in range(12):
            s = PN.reduce((s + u * s.inverse_mod(PN)) * inv2)
    assert (s * s - u).valuation(pl.P) >= N, 'square root lost precision'
    return s * pl.pi ** (v // 2)

def _fac(x):
    """unit and Gaussian prime exponents of a non-zero element of Z[i]."""
    F = K(x).factor()
    return F.unit(), {K(p): int(e) for p, e in F}

def conic_point_lattice(co, rng, tries=6):
    """a point of A u^2 + B v^2 + C z^2 = 0 over K, by a lattice search over Z[i] (new code; PARI's qflll and qfminim
    as search tools): reduce to A, B, C in Z[i] square-free and pairwise coprime, impose the linear congruences
    Q = 0 mod pi at the primes of ABC (the index-N(ABC) sublattice of Z[i]^3 = Z^6), LLL-reduce with the weights |A|,
    |B|, |C| and enumerate short vectors.  The point is checked exactly on the original conic."""
    import cypari2
    pari = cypari2.Pari()
    A, B, C = [K(c) for c in co]
    den = lcm([c.denominator() for c in (A, B, C)]); A, B, C = A * den, B * den, C * den
    f = [K(1), K(1), K(1)]                                   # original coordinates = f_k * reduced coordinates
    co_ = [A, B, C]
    g = K.ideal(co_).gens_reduced()[0]; co_ = [c / g for c in co_]
    # square-free
    facs = []
    for k in range(3):
        u, ex = _fac(co_[k]); s = prod([p ** (e // 2) for p, e in ex.items()]) if ex else K(1)
        co_[k] = co_[k] / s ** 2; f[k] = f[k] / s
    # pairwise coprime
    changed = True
    while changed:
        changed = False
        exs = [_fac(c)[1] for c in co_]
        for p in set().union(*[set(e) for e in exs]):
            ins = [k for k in range(3) if any((p / q) in OK_ and (q / p) in OK_ for q in exs[k])]
            if len(ins) == 3:
                co_ = [c / p for c in co_]; changed = True; break
            if len(ins) == 2:
                k = [j for j in range(3) if j not in ins][0]
                co_[ins[0]] /= p; co_[ins[1]] /= p; co_[k] *= p; f[k] *= p; changed = True; break
        if changed:
            # re-square-free (the third coefficient may now contain p^2)
            for k in range(3):
                u, ex = _fac(co_[k]); s = prod([q ** (e // 2) for q, e in ex.items()]) if ex else K(1)
                co_[k] = co_[k] / s ** 2; f[k] = f[k] / s
    A, B, C = co_
    # congruences: for pi | A: v = t z mod pi with B t^2 + C = 0; pi | B: z = t u (C t^2 + A = 0); pi | C: u = t v
    conds = []          # (pi, coefficient vector (c_u, c_v, c_z)) meaning c . (u, v, z) = 0 mod pi
    for k, (X, Y, Z) in enumerate([(A, B, C), (B, C, A), (C, A, B)]):
        for p, e in _fac(X)[1].items():
            F = K.residue_field(K.ideal(p))
            r = -F(Z) / F(Y)
            assert r.is_square(), 'the conic is not locally soluble at %s' % p
            t = F.lift(r.sqrt())
            vec = [K(0)] * 3
            # pi | A: v - t z;  pi | B: z - t u;  pi | C: u - t v
            j1, j2 = [(1, 2), (2, 0), (0, 1)][k]
            vec[j1] = K(1); vec[j2] = -t
            conds.append((p, vec))
    # the sublattice of Z^6 (coordinates re u, im u, re v, im v, re z, im z)
    Bm = pari.matid(6)
    for p, vec in conds:
        q = ZZ(p.norm())
        P = K.ideal(p)
        if q.is_prime():
            pp = q; F = K.residue_field(P); r = ZZ(F.lift(F(I)).list()[0]) % pp if True else None
            # i -> r mod p (r^2 = -1 mod p, or r = 1 mod 2 for 1 + i)
            rr = int((F.lift(F(I)) - 0).list()[0]) % pp
            row = []
            for c in vec:
                cr, ci = [int(x) for x in c.list()]
                # c * (x + i y) -> (cr + ci r) x + (cr r - ci) y mod p   (i^2 = -1)
                row += [(cr + ci * rr) % pp, (cr * rr - ci) % pp]
            G = pari.matrix(1, 6, row)
        else:
            pp = ZZ(q.sqrt()); rows = [[], []]
            for c in vec:
                cr, ci = [int(x) for x in c.list()]
                rows[0] += [cr % pp, (-ci) % pp]; rows[1] += [ci % pp, cr % pp]
            G = pari.matrix(2, 6, rows[0] + rows[1])
        GB = pari('lift(Mod(%s * %s, %d))' % (G, Bm, pp))
        Kr = pari.matkermod(GB, pp)
        Mx = pari.concat(pari.lift(Kr), pari.matid(6) * pp)
        Hn = pari.mathnf(Mx)
        Bm = Bm * Hn
    w = [max(1, int(abs(ComplexField(60)(c)) + 0.5)) for c in (A, B, C)]
    D = pari.matdiagonal([w[0], w[0], w[1], w[1], w[2], w[2]])
    Gram = pari.mattranspose(Bm) * D * Bm
    U = pari.qflllgram(Gram)
    Bl = Bm * U; Gl = pari.mattranspose(U) * Gram * U
    def qv(x):
        u = K(int(x[0])) + K(int(x[1])) * I; v = K(int(x[2])) + K(int(x[3])) * I; z = K(int(x[4])) + K(int(x[5])) * I
        return u, v, z
    bound = int(Gl[0, 0]) * 2
    for _ in range(tries):
        try:
            res = pari.qfminim(Gl, bound, 4000)
        except Exception:
            pari.set_real_precision(300); res = pari.qfminim(Gl, bound, 4000, 2)
        vecs = res[2]
        for j in range(int(pari.matsize(vecs)[1])):
            x = Bl * vecs[j]
            u, v, z = qv([x[t] for t in range(6)])
            if (u, v, z) != (0, 0, 0) and A * u * u + B * v * v + C * z * z == 0:
                pt = [f[0] * u, f[1] * v, f[2] * z]
                assert sum(c * p * p for c, p in zip(co, pt)) == 0
                return pt
        bound *= 3
    return None

def conic_point(co, rng=None):
    pt = conic_point_lattice(co, rng)
    if pt is None:
        C = Conic(K, co); pt = [K(c) for c in C.rational_point()]     # fallback (slow: class groups)
    den = lcm([c.denominator() for c in pt]); pt = [c * den for c in pt]
    g = K.ideal(pt).gens_reduced()[0]; pt = [c / g for c in pt]
    assert sum(c * p * p for c, p in zip(co, pt)) == 0, 'the conic point is not on the conic'
    return pt

def primitive(v):
    v = [K(c) for c in v]
    den = lcm([c.denominator() for c in v]); v = [c * den for c in v]
    g = K.ideal(v).gens_reduced()[0]
    return [c / g for c in v]

def tangent(co, pt):
    return primitive([c * p for c, p in zip(co, pt)])

def lift_data(sel, a, rng, N=40, forged_line=None, max_used=1500):
    """for a = (a1, a2, a3) in Sel^2 (a3 = a1 a2 mod squares): the affine spaces g(D_a(K_v)) at the places of S, and
    whether a global lam exists."""
    es = sel.es; e1, e2, e3 = es
    a1, a2, a3 = a
    Cw = [a1, -a2, e1 - e2]; Cv = [a1, -a3, e1 - e3]; Cu = [a2, -a3, e2 - e3]
    pw = conic_point(Cw); pv = conic_point(Cv); pu = conic_point(Cu)
    lw = tangent(Cw, pw); lv = tangent(Cv, pv); lu = tangent(Cu, pu)
    if forged_line is not None: lw = forged_line
    def g_of(u, v, w):
        Lw = lw[0] * u + lw[1] * v + lw[2]; Lv = lv[0] * u + lv[1] * w + lv[2]; Lu = lu[0] * v + lu[1] * w + lu[2]
        return Lu * Lw, Lv * Lw
    spaces = {}
    for pl in sel.places:
        target = tuple(pl.cls(a1) + pl.cls(a2))
        vals = []
        cand = [x for x, vec in sel.pts[pl.P] if vec == target and x not in es]
        rng.shuffle(cand)
        want = 4 if pl.two else 2
        def more():
            # further local points in the class of a: x random, near the roots, near 0 and near infinity
            for k in range(200000):
                t = rng.randint(0, 3); z = K(rng.randint(-500, 500) + rng.randint(-500, 500) * I); j = rng.randint(0, 12)
                x = [z, es[rng.randint(0, 2)] + pl.pi ** j * z, pl.pi ** j * z, z / pl.pi ** (2 * j)][t]
                if x in es: continue
                fx = (x - e1) * (x - e2) * (x - e3)
                if pl.is_square(fx) and tuple(pl.cls(x - e1) + pl.cls(x - e2)) == target: yield x
        src_ = itertools.chain(cand, more())
        used = 0
        for x in src_:
            used += 1
            if used > max_used: break
            r1 = (x - e1) / a1; r2 = (x - e2) / a2; r3 = (x - e3) / a3
            if not (pl.is_square(r1) and pl.is_square(r2) and pl.is_square(r3)): continue
            u0 = sqrt_local(r1, pl, N); v0 = sqrt_local(r2, pl, N); w0 = sqrt_local(r3, pl, N)
            for su, sv, sw in itertools.product((1, -1), repeat=3):
                G1, G2 = g_of(su * u0, sv * v0, sw * w0)
                if G1 == 0 or G2 == 0: continue
                if max(pl.val(G1), pl.val(G2)) > N - 12: continue      # too close to a zero of a line
                vals.append(tuple(pl.cls(G1) + pl.cls(G2)))
            if used >= 8 and vals and len(set(vals)) == 2 ** want and \
                    f2_rank([[(p ^ q) for p, q in zip(v_, vals[0])] for v_ in vals]) == want: break
        assert vals, 'no local point of D_a found at %s' % pl.P
        base = vals[0]
        diffs = [[(p ^ q) for p, q in zip(v_, base)] for v_ in vals]
        d = f2_rank(diffs)
        if len(set(vals)) != 2 ** d: d = -d            # the values are not a whole affine space: flagged
        spaces[pl.P] = (base, diffs, d)
    # global lam = (lam1, lam2) in K(S,2)^2:  res_v(lam) - base_v in span(diffs_v) for all v
    rows = []; rhs = []
    for pl in sel.places:
        base, diffs, d = spaces[pl.P]
        W = matrix(GF(2), diffs) if any(any(r) for r in diffs) else matrix(GF(2), 0, 2 * pl.w)
        ann = W.right_kernel().basis() if W.nrows() else matrix.identity(GF(2), 2 * pl.w).rows()
        res = [pl.cls(g) for g in sel.gens]
        for h in ann:
            h = list(h); row = []
            for c in range(2):
                for j in range(sel.n):
                    row.append(sum(h[c * pl.w + t] * res[j][t] for t in range(pl.w)) % 2)
            rows.append(row); rhs.append(sum(int(h[t]) * base[t] for t in range(2 * pl.w)) % 2)
    A = matrix(GF(2), rows); bvec = vector(GF(2), rhs)
    try:
        sol = A.solve_right(bvec); ok = True
    except ValueError:
        sol = None; ok = False
    dims = {str(pl.P.gens_reduced()[0]): spaces[pl.P][2] for pl in sel.places}
    return ok, dims, (pw, pv, pu), (lw, lv, lu), sol

# ------------------------------------------------------------------ torsion
def naive_count(Ered, F):
    a1, a2, a3, a4, a6 = Ered.a_invariants()
    cnt = 1
    for x in F:
        rhs = x ** 3 + a2 * x ** 2 + a4 * x + a6
        cnt += 1 if rhs == 0 else (2 if rhs.is_square() else 0)
    return cnt

def torsion_check(E, es):
    Ts = [E([ee, 0]) for ee in es]
    assert all(2 * T == 0 and T != 0 for T in Ts) and Ts[0] + Ts[1] == Ts[2]
    halves = []
    for T in Ts:
        try:
            hs = T.division_points(2)
        except Exception:
            hs = []
        for H in hs:
            assert 2 * H == T and H in E; halves.append((T, H))
    orders = []
    g = 0; used = []
    for p in primes(3, 400):
        for P, _ in K.ideal(p).factor():
            if E.has_good_reduction(P) and P.norm() < 10 ** 4:
                Er = E.reduction(P); c = Er.cardinality(); g = gcd(g, c); used.append((P.norm(), c))
    naive = []
    for P, _ in K.ideal(13).factor():
        Er = E.reduction(P); naive.append((Er.cardinality(), naive_count(Er, Er.base_field())))
    for P, _ in K.ideal(7).factor():
        Er = E.reduction(P); naive.append((Er.cardinality(), naive_count(Er, Er.base_field())))
    assert all(a == b for a, b in naive)
    return Ts, halves, g, used[:8], naive

# ------------------------------------------------------------------ runs
def target(rng):
    m, n, J = 627, 364, 'L1L3M1M2'
    say('kctp_recheck.py: independent re-check of the inputs of kctp.py for %s of (%d,%d) over K = Q(i).' % (J, m, n))
    say('Written by an AI model; not read by a person.  Experiment, not counted.  Sage (passagemath) for exact arithmetic.')
    say()
    Q, fs = quartic(m, n, J)
    say('1. The curve.  Q_J = %s' % ' * '.join('(%s)' % f for f in fs))
    Ejac = jac_invariants(Q)
    es = full2_model(Q, fs)
    E = curve_from_e(es)
    Ek = curve_from_e([K(a) + K(b) * I for a, b in E_KCTP[J]])
    assert E.is_isomorphic(Ejac) and E.is_isomorphic(Ek) and Ejac.is_isomorphic(Ek)
    j = E.j_invariant(); assert j not in QQ
    say('   Jacobian of y^2 = Q_J (invariants I, J), the model with the root of M2 at infinity, and the model of kctp.txt:')
    jc = list(j)
    say('   isomorphic over K (exact).  j = c0 + c1 i with c1 = %.6e != 0 (exact test c1 != 0), so j is not in Q.' % float(jc[1]))
    conj = lambda F: EllipticCurve(K, [K(list(c)[0]) - K(list(c)[1]) * I for c in F.a_invariants()])
    rel = {}
    for J2 in ('L1L3M3M4', 'L2L4M1M2', 'L2L4M3M4'):
        Q2, fs2 = quartic(m, n, J2); F2 = jac_invariants(Q2)
        rel[J2] = 'E' if F2.is_isomorphic(E) else ('conj(E)' if F2.is_isomorphic(conj(E)) else 'neither')
    say('   the other corner curves (Jacobians of their quartics) are isomorphic over K to: %s' % rel)
    assert set(rel.values()) <= {'E', 'conj(E)'}
    say('   (a curve and its complex conjugate have the same rank over K), so the four have one rank between them.')
    for t in (1, -1):
        val = Q(-t * I)
        assert val != 0 and val.is_square()
    say('   Q_J(-i) and Q_J(i) (t = 1, -1) are non-zero squares in K (exact): the degenerate solutions give K-points.')
    say('   e = %s' % [str(x) for x in es])
    Ts, halves, g, used, naive = torsion_check(E, es)
    say('2. Torsion.  T1, T2, T3 = (e_k, 0) have order 2; halves found and checked (2H = T): %d, all halving T3.' % len(halves))
    assert len(halves) == 4 and all(T == Ts[2] for T, H in halves)
    say('   gcd of #E(F_P) over the primes of good reduction with N(P) < 10^4 and p < 400: %d (e.g. %s).' % (g, used[:5]))
    say('   naive point counts agree with Sage at the primes above 13 and 7: %s' % naive)
    assert g == 8
    say('   So E(K)_tors has order 8 and contains Z/2 x Z/4: E(K)_tors = Z/2 x Z/4 (exact; torsion injects into E(F_P)')
    say('   for P of good reduction above an odd p, as e(P|p) = 1 < p - 1).')
    t0 = time.time()
    sel = Selmer(es, rng)
    say('3. Sel^2(E/K), new code.  S = places above %s (%d places); K(S,2) has dimension %d.' % (
        sorted(set(int(ZZ(pl.P.norm()).prime_factors()[0]) for pl in sel.places)), len(sel.places), sel.n))
    say('   local images: dimensions %s (2 at odd places, 4 at 1 + i)' % [len(sel.img[pl.P]) for pl in sel.places])
    say('   dim Sel^2(E/K) = %d   (%.0f s)' % (sel.dim, time.time() - t0))
    assert sel.dim == 4
    timg = []
    for x in [es[0], es[1]] + [H.xy()[0] for T, H in halves]:
        d1, d2 = kummer(x, es); timg.append(sel.vec_of(d1) + sel.vec_of(d2))
    say('   torsion image (T1, T2, the halves of T3) has dimension %d' % f2_rank(timg))
    assert f2_rank(timg) == 2
    assert all(sel.contains(*sel.pair(b)[:2]) for b in sel.basis)
    # kctp.txt's four Selmer elements: (i 991, i 991), (... ), ...
    G = lambda a, b: K(a) + K(b) * I
    v7 = [(G(0, 1) * G(991, 0), G(0, 1) * G(991, 0)),
          (G(1, 1) * G(1, -6) * G(9, 10) * G(991, 0) * G(671, -4800), G(1, 1) * G(991, 0) * G(95732, 621357)),
          (G(1, 2) * G(1, -2) * G(2, 5) * G(2, -5), G(1, 2) * G(1, -2) * G(2, 5) * G(2, -5)),
          (G(263, 0), G(263, 0))]
    vv = []
    for d1, d2 in v7:
        assert sel.contains(d1, d2), 'a Selmer element of kctp.txt fails the local test here'
        vv.append(sel.vec_of(d1) + sel.vec_of(d2))
    say('   the four Selmer elements printed in kctp.txt pass the local test here and span a space of dimension %d,' % f2_rank(vv))
    say('   so the two computations give the same group.')
    assert f2_rank(vv) == 4 and f2_rank(vv + sel.basis) == 4
    say('4. The lift criterion (second descent, no Hilbert symbols) for a basis of Sel^2:')
    allok = True
    for b in sel.basis:
        a = sel.pair(b)
        ok, dims, pts, lines, sol = lift_data(sel, a, rng)
        allok &= ok
        say('   class %s: g(D_a(K_v)) affine of dimension %s; global lam exists: %s' % (''.join(map(str, b)), sorted(set(dims.values())), ok))
        expect = {str(pl.P.gens_reduced()[0]): (4 if pl.two else 2) for pl in sel.places}
        assert dims == expect, 'g(D_a(K_v)) has the wrong size: %s' % dims
    assert allok
    say('   every class of Sel^2 lifts to an everywhere locally soluble 4-covering: the Cassels-Tate pairing on Sel^2')
    say('   is zero (Cassels; quoted).  This agrees with kctp.py.')
    say()
    say('RESULT: kctp.py\'s inputs re-checked with new code: the curve (three models, isomorphic over K, j not in Q),')
    say('E(K)_tors = Z/2 x Z/4, dim Sel^2 = 4 with torsion image of dimension 2, the same group as kctp.txt, and every')
    say('class lifts (pairing zero).  So rank E(K) <= 2; 0 versus 2 is not decided here.')
    return sel

def bc_curve(a, b):
    """y^2 = x (x + a)(x + b) over Q, base-changed to K."""
    return [K(0), K(-a), K(-b)]

def lift_kernel_dim(sel, rng):
    """dimension of the subspace of Sel^2 of classes that lift (the kernel of the pairing): it contains the torsion
    image; the classes of a complement are tested in all combinations."""
    es = sel.es
    tor = []
    E = curve_from_e(es)
    for x in es[:2]:
        d1, d2 = kummer(x, es); tor.append(sel.vec_of(d1) + sel.vec_of(d2))
    for t in [E([ee, 0]) for ee in es]:
        try: hs = t.division_points(2)
        except Exception: hs = []
        for H in hs:
            d1, d2 = kummer(H.xy()[0], es); tor.append(sel.vec_of(d1) + sel.vec_of(d2))
    Tb = [list(b) for b in matrix(GF(2), tor).row_space().basis()]
    comp = []
    for b in sel.basis:
        if f2_rank(Tb + comp + [b]) > len(Tb) + len(comp): comp.append(b)
    lifting = []
    for mask in range(1, 2 ** len(comp)):
        v = [0] * (2 * sel.n)
        for j in range(len(comp)):
            if mask >> j & 1: v = [(p + q) % 2 for p, q in zip(v, comp[j])]
        ok, dims, _, _, _ = lift_data(sel, sel.pair(v), rng)
        expect = {str(pl.P.gens_reduced()[0]): (4 if pl.two else 2) for pl in sel.places}
        assert dims == expect, 'g(D_a(K_v)) has the wrong size: %s' % dims
        if ok: lifting.append(v)
    kd = f2_rank(Tb + lifting)
    return kd, len(Tb)

# control curves y^2 = x (x + a)(x + b) of kctp_controls.txt (run 3): (a, b, dim Sel^2_K, pairing rank found by kctp.py)
BC = [(14, 74, 4, 2), (-56, -71, 4, 2), (66, 9, 4, 2), (59, -47, 5, 2),
      (41, 80, 4, 0), (61, -21, 4, 0), (29, 21, 4, 0), (-20, 71, 3, 0)]
# curves y^2 = x (x - c)(x - d) over K with j not in Q (kctp_controls.txt): (c, d, dim Sel^2, pairing rank, lower bound)
KC = [((14, -33), (32, -25), 4, 2, 0), ((-35, -31), (31, 33), 4, 2, 0),
      ((14, 0), (19, 34), 4, 0, 2), ((0, 3), (4, 36), 4, 0, 2)]

def controls(rng):
    say('kctp_recheck.py --controls.  For each curve: dim Sel^2(E/K) and the dimension of the subspace of classes that')
    say('lift (kernel of the Cassels-Tate pairing, by the lift criterion), with the new code; compared with kctp.py\'s')
    say('results (run 3, kctp_controls.txt: dim Sel^2 and the rank of the pairing matrix) and with the known ranks.')
    import cypari2
    pari = cypari2.Pari()
    say('(1) base changes of y^2 = x(x + a)(x + b) from Q; rank E(K) = rank E(Q) + rank E^-1(Q) from PARI ellrank (guide,')
    say('    used when its bounds agree).  Required: bound = lifting - 2 >= rank_K, and agreement with kctp.py.')
    n = 0; agree = 0; exact = 0
    for a, b, dv7, rho in BC:
        rk = []
        for sg in (1, -1):
            r = pari.ellrank(pari.ellinit([0, sg * (a + b), 0, a * b, 0]))
            rk.append((int(r[0]), int(r[1])))
        rK = rk[0][0] + rk[1][0] if all(lo == hi for lo, hi in rk) else None
        es = [K(0), K(-a), K(-b)]
        sel = Selmer(es, rng)
        kd, td = lift_kernel_dim(sel, rng)
        ok7 = (sel.dim == dv7 and kd == dv7 - rho)
        n += 1; agree += ok7
        sound = rK is None or kd - td >= rK
        exact += (rK is not None and kd - td == rK)
        say('  a=%4d b=%4d: ranks over Q %s -> rank_K %s; dim Sel^2_K %d, lifting %d, bound %d; kctp.py: dim %d, pairing rank %d: %s%s' % (
            a, b, rk, rK, sel.dim, kd, kd - td, dv7, rho, 'agree' if ok7 else 'DISAGREE', '' if sound else '  UNSOUND'))
        assert sound and ok7
    say('  %d curves: sound on all; bound = rank_K on %d; agreement with kctp.py on %d.' % (n, exact, agree))
    say('(2) curves over K with j not in Q (kctp.py\'s list): lower bound from points found by kctp.py (not re-checked')
    say('    here); required: agreement with kctp.py and bound >= that lower bound.')
    m = 0
    for c, d, dv7, rho, lo in KC:
        es = [K(0), K(c[0]) + K(c[1]) * I, K(d[0]) + K(d[1]) * I]
        sel = Selmer(es, rng)
        kd, td = lift_kernel_dim(sel, rng)
        ok7 = (sel.dim == dv7 and kd == dv7 - rho)
        say('  c=%s d=%s: dim Sel^2 %d, torsion image %d, lifting %d -> bound %d; kctp.py: dim %d, pairing rank %d, lower %d: %s' % (
            c, d, sel.dim, td, kd, kd - td, dv7, rho, lo, 'agree' if ok7 else 'DISAGREE'))
        assert ok7 and kd - td >= lo
        m += 1
    say('  %d curves: agreement with kctp.py on all; the lift criterion fails (pairing non-zero) on %d of them.' % (
        m, sum(1 for x in KC if x[3] > 0)))
    return n + m

def negative(rng):
    say('kctp_recheck.py --negative: forged inputs that must be refused.')
    m, n, J = 627, 364, 'L1L3M1M2'
    Q, fs = quartic(m, n, J); es = full2_model(Q, fs)
    E = curve_from_e(es)
    refused = 0; total = 0
    def expect_fail(name, fn):
        nonlocal refused, total
        total += 1
        try:
            fn(); say('  NOT REFUSED: %s' % name)
        except AssertionError as ex:
            refused += 1; say('  refused: %s (%s)' % (name, str(ex)[:70]))
    # 1. a mutated model is not isomorphic to the curve of kctp.txt
    def f1():
        es2 = [es[0] + 1, es[1], es[2]]
        assert curve_from_e(es2).is_isomorphic(curve_from_e([K(a) + K(b) * I for a, b in E_KCTP[J]])), 'not isomorphic'
    expect_fail('a model with e1 + 1', f1)
    # 2. a forged torsion point
    def f2():
        assert E.is_on_curve(es[0] + 1, 0), 'not on the curve'
    expect_fail('the point (e1 + 1, 0)', f2)
    sel = Selmer(es, rng)
    # 3. forged Selmer elements: (pi, 1) for primes of S
    for pl in sel.places[:6]:
        def f3(pl=pl):
            assert sel.contains(pl.pi, K(1)), 'fails the local test'
        expect_fail('(%s, 1) as a Selmer element' % pl.pi, f3)
    # 4. a non-S-unit
    expect_fail('(3, 3) (3 not in S)', lambda: (_ for _ in ()).throw(AssertionError('fails')) if not sel.contains(K(3), K(3)) else None)
    # 5. a forged conic point
    def f5():
        co = [K(1), K(1), K(-3)]
        pt = [K(1), K(1), K(1)]
        assert sum(c * p * p for c, p in zip(co, pt)) == 0, 'the conic point is not on the conic'
    expect_fail('a point off a conic', f5)
    # 6. a non-tangent line in the lift criterion: the image is no longer an affine space of the right size
    def f6():
        b = sel.basis[0]; a = sel.pair(b)
        ok, dims, _, _, _ = lift_data(sel, a, rng, forged_line=[K(1), K(2), K(3)], max_used=60)
        expect = {str(pl.P.gens_reduced()[0]): (4 if pl.two else 2) for pl in sel.places}
        assert dims == expect and ok, 'g(D_a(K_v)) has the wrong size or no lam'
    expect_fail('a non-tangent line l_w in the lift criterion', f6)
    # 7. a forged Selmer basis element: a class outside Sel^2 used as D_a (no local points)
    def f7():
        a = (sel.places[1].pi, K(1), sel.places[1].pi)
        ok, dims, _, _, _ = lift_data(sel, a, rng, max_used=60)
    expect_fail('the lift criterion for a class outside Sel^2', f7)
    say('%d of %d forged inputs refused.' % (refused, total))
    assert refused == total

if __name__ == '__main__':
    rng = random.Random(20261125)
    here = os.path.dirname(os.path.abspath(__file__))
    if '--controls' in sys.argv:
        controls(rng); fn = 'kctp_recheck_controls.txt'
    elif '--negative' in sys.argv:
        negative(rng); fn = 'kctp_recheck_negative.txt'
    else:
        target(rng); fn = 'kctp_recheck.txt'
    if WRITE:
        open(os.path.join(here, fn), 'w').write('\n'.join(OUT) + '\n')
