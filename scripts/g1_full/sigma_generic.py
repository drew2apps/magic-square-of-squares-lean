#!/usr/bin/env python3
"""Uniform rank bounds for the curves E+ and E- of F0 F8 GL (UNIFORM_RANK.md): rank >= 1 at every ratio, rank >= 2 at
every ratio with m^2 + n^2 twice a square.  Exact arithmetic only: sympy for the identities in general m, n, Python
integers, Fractions and pairs (a, b) = a + b w (w^2 = 2 or -1) for the descents.  PARI only in --pari, as a guide.
Written by an AI model; not read by a person.

  python3 sigma_generic.py [--pari] [--write]   # all parts (about 40 s; --pari adds PARI as a guide); --write writes sigma_generic.txt
  python3 sigma_generic.py --negative           # the negative tests and the check of the local search only

Part A  identities in general m, n: the model of E+ (E- at (m, n) is E+ at (m, -n)), its 2-torsion, the points P, R
        (2R = P), Q, the values of the 2-isogeny map, the identity for -theta, the duplication identity.  A9 to A11
        (the Q(i)-classes on the conic) are kept for UNIFORM_RANK.md, section 6; the proof does not use them.
Part B  Lemma B (Euler, d = 14): m^8 + 14m^4n^4 + n^8 is a square only if mn(m^2 - n^2) = 0; 2-isogeny descent over Q.
Part C  Lemma D: w^2 = 2(x^2 + 1) H(x, 1) has only the rational points (+-1, +-4); 2-isogeny descent over Q(sqrt2)
        on an elliptic quotient, torsion by reduction, and the pull-back of the ten points.
Part D  the conclusion checked curve by curve from the model of sigma_exact.py: the 127 ratios m <= 1000 of
        sigma_twice.py, 200 random twice-square ratios with m <= 10^6, the square type, all ratios with m < 60.
Part E  negative tests: forged inputs that the checks must refuse; the local search against a small-point search.
"""
import sys, random
from fractions import Fraction as Fr
from math import gcd, isqrt

LOG = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); LOG.append(s)

# ---------------------------------------------------------------- Part A: identities (sympy)
def part_A():
    from sympy import symbols, expand, I, diff, Poly, conjugate, simplify, cancel, together
    m, n, z, X, h, p, q = symbols('m n z X h p q')
    S = m**2 + n**2
    g = (2*S*z - 8*m*n)*(m**2*n**2*z**2 + (m**2 - n**2)**2)
    def wmodel(s, mm=m, nn=n):
        gg = g.subs({m: mm, n: nn}, simultaneous=True)
        z0 = -s
        qz = expand((z - z0)*gg)
        qh = Poly(expand(qz.subs(z, z0 + h)), h)
        c = [qh.coeff_monomial(h**k) for k in range(5)]
        assert c[0] == 0
        return z0, c, (c[2], c[1]*c[3], c[1]**2*c[4])
    z0, c, (a2, a4, a6) = wmodel(2)
    A = c[1]
    cub = X**3 + a2*X**2 + a4*X + a6
    # A1: A = -4(m+n)^2 S^2 on E+
    assert expand(A + 4*(m + n)**2*S**2) == 0
    say('A1  E+: y^2 = X^3 + a2 X^2 + a4 X + a6 by the recipe of sigma_exact.model(m, n, 2), in general m, n;', 'A = -4(m+n)^2 S^2')
    # A2: E-(m, n) has literally the same Weierstrass model as E+(m, -n)
    _, cm, mod_m = wmodel(-2)
    _, cp, mod_p = wmodel(2, m, -n)
    assert all(expand(u - v) == 0 for u, v in zip(mod_m, mod_p))
    say('A2  the model of E- at (m, n) equals the model of E+ at (m, -n)')
    # A3: the roots of the cubic
    e = -2*S**3
    th = 4*I*m*n*(m + n)**2*(m + I*n)**2
    thb = -4*I*m*n*(m + n)**2*(m - I*n)**2
    assert expand(cub - (X - e)*(X - th)*(X - thb)) == 0
    say('A3  X^3 + a2 X^2 + a4 X + a6 = (X - e)(X - theta)(X - conj theta), e = -2S^3, theta = 4i mn (m+n)^2 (m+in)^2')
    # A4: imaginary part of theta is 4mn(m+n)^2 (m^2 - n^2), nonzero for mn(m^2-n^2) != 0, so E(Q)[2] = {O, T}
    assert expand((th - thb)/(2*I) - 4*m*n*(m + n)**2*(m**2 - n**2)) == 0
    u_ = 16*m**2*n**2*(m + n)**2; v_ = 16*m**2*n**2*(m + n)**4*S**2
    assert expand(cub - (X - e)*(X**2 + u_*X + v_)) == 0
    assert expand(u_**2 - 4*v_ + 64*m**2*n**2*(m - n)**2*(m + n)**6) == 0
    say('A4  the cubic is (X + 2S^3)(X^2 + 16m^2n^2(m+n)^2 X + 16m^2n^2(m+n)^4 S^2), the quadratic factor has',
        'discriminant -64 m^2n^2(m-n)^2(m+n)^6 < 0: E(Q)[2] = {O, T}, T = (-2S^3, 0), when mn(m^2-n^2) != 0')
    # A5: P (t = 1, z = 2) and R with 2R = P
    xP = A/4
    assert expand(xP + (m + n)**2*S**2) == 0
    yP = S**3*(m - n)*(m + n)**2
    assert expand(cub.subs(X, xP) - yP**2) == 0
    Hf = lambda mm, nn: mm**4 - 2*mm**3*nn + 2*mm**2*nn**2 + 2*mm*nn**3 + nn**4
    H = Hf(m, n)
    xR = -4*m*n*S*(m**2 - n**2)
    yR = -8*m*n*(m + n)*S*H
    assert expand(cub.subs(X, xR) - yR**2) == 0
    b2, b4, b6, b8 = 4*a2, 2*a4, 4*a6, 4*a2*a6 - a4**2
    num = xR**4 - b4*xR**2 - 2*b6*xR - b8
    den = 4*xR**3 + b2*xR**2 + 2*b4*xR + b6
    assert expand(num - xP*den) == 0
    # the tangent at R meets the curve again at -P: the slope lam gives y(2R) = -(lam (x(2R) - xR) + yR)
    lam = (3*xR**2 + 2*a2*xR + a4)/(2*yR)
    y2R = cancel(-(lam*(xP - xR) + yR))
    assert expand(cancel(y2R - yP)) == 0
    say('A5  P = (-(m+n)^2 S^2, S^3 (m-n)(m+n)^2) and R = (-4mnS(m^2-n^2), -8mn(m+n)S H(m,n)) are on E+ and 2R = P,',
        'for every ratio; H(m,n) = m^4 - 2m^3n + 2m^2n^2 + 2mn^3 + n^4')
    # A6: the 2-isogeny values (shifted model x -> x - e)
    assert expand(xR - e - 2*S*H) == 0
    fprime_e = expand(diff(cub, X).subs(X, e))
    assert expand(fprime_e - 4*S**2*(m**8 + 14*m**4*n**4 + n**8)) == 0
    assert expand(Hf(m, n)*Hf(m, -n) - (m**8 + 14*m**4*n**4 + n**8)) == 0
    assert expand(fprime_e - (e - th)*(e - thb)) == 0
    say('A6  x(R) - e = 2S H(m,n);  f\'(e) = (e - theta)(e - conj theta) = 4S^2 (m^8 + 14m^4n^4 + n^8)',
        '= 4S^2 H(m,n) H(m,-n)')
    # A7: Q: the constant term a6 = A^2 * 2 m^2 n^2 S, so X = 0 is a rational point iff 2S is a square
    assert expand(a6 - A**2*2*m**2*n**2*S) == 0
    say('A7  a6 = A^2 (2S) m^2 n^2: Q = (0, A mn sqrt(2S)) is rational iff 2S is a square')
    # A8: the Q(i)-classes, general m, n
    w = (m**2 - 2*m*n - n**2) + I*S
    dR = xR - th
    assert expand(dR + 4*m*n*(m + n)*(m + I*n)*w) == 0
    dQ = -th
    dT = e - th
    quartT = m**4 - n**4 + 4*I*m**2*n**2
    assert expand(dT + 2*(m + I*n)**2*quartT) == 0
    say('A8  x(R) - theta = -4mn(m+n)(m+in) w,  w = (m^2-2mn-n^2) + i S;  -theta = -4i mn (m+n)^2 (m+in)^2;',
        'e - theta = -2 (m+in)^2 (m^4 - n^4 + 4i m^2 n^2)')
    # A9: on the conic
    sub = {m: p**2 + 2*p*q - q**2, n: -p**2 + 2*p*q + q**2}
    Sp = expand(S.subs(sub, simultaneous=True))
    assert expand(Sp - 2*(p**2 + q**2)**2) == 0
    assert expand((m + n).subs(sub, simultaneous=True) - 4*p*q) == 0
    assert expand((m + I*n).subs(sub, simultaneous=True) - (1 - I)*(p + I*q)**2) == 0
    W = p**4 + (2 - 2*I)*p**3*q + (-2 + 4*I)*p**2*q**2 + (-2 + 2*I)*p*q**3 + q**4
    Wp = p**4 + (-2 - 2*I)*p**3*q + (-2 - 4*I)*p**2*q**2 + (2 + 2*I)*p*q**3 + q**4
    assert expand(w.subs(sub, simultaneous=True) - (2 + 2*I)*W) == 0
    assert expand(quartT.subs(sub, simultaneous=True) - 4*I*W*Wp) == 0
    # W'(p, q) = conj(W)(p, -q)
    Wc = W.subs(I, -I)
    assert expand(Wp - Wc.subs(q, -q)) == 0
    # the classes with explicit square cofactors (-1 = i^2, -8i = -4 (1+i)^2, 2+2i times 1-i is 4)
    mm, nn = sub[m], sub[n]
    dRp = expand(dR.subs(sub, simultaneous=True))
    assert expand(dRp - (8*I*(p + I*q))**2*(mm*nn*p*q*W)) == 0
    # -theta = -4i mn (m+n)^2 (m+in)^2 = (2i(m+n)(m+in))^2 * (i mn)
    assert expand(dQ - (2*I*(m + n)*(m + I*n))**2*(I*m*n)) == 0
    # e - theta = -2(m+in)^2 * 4i W W' = (2(1+i) i (m+in))^2 * W W'   [(2(1+i)i)^2 = -8i]
    dTp = expand(dT.subs(sub, simultaneous=True))
    assert expand((2*(1 + I)*I)**2 + 8*I) == 0
    assert expand(dTp - (2*(1 + I)*I*(1 - I)*(p + I*q)**2)**2*W*Wp) == 0
    say('A9  on the conic: m+n = 4pq, m+in = (1-i)(p+iq)^2, w = (2+2i) W, m^4-n^4+4i m^2n^2 = 4i W W\',',
        'W\'(p,q) = conj(W)(p,-q); hence, up to squares in Q(i):  delta(R) = mn pq W,  delta(Q) = i mn,  delta(T) = W W\'')
    # A10: norms
    NW = expand(W*Wc)
    NWp = expand(Wp*Wp.subs(I, -I))
    assert expand(NW - Hf(mm, nn)/4) == 0 and expand(NWp - Hf(mm, -nn)/4) == 0
    say('A10 N(W) = H(m,n)/4, N(W\') = H(m,-n)/4 on the conic')
    # A11: resultants of the atoms (common prime factors of two forms at coprime (p,q) divide the resultant)
    from sympy import resultant, Symbol
    t = Symbol('t')
    atoms = {'m': mm, 'n': nn, 'W': W, "W'": Wp}
    names = list(atoms)
    for i1 in range(len(names)):
        for j1 in range(i1 + 1, len(names)):
            F = expand(atoms[names[i1]].subs({p: t, q: 1})); G = expand(atoms[names[j1]].subs({p: t, q: 1}))
            r = expand(resultant(F, G, t))
            from sympy import re as Re, im as Im
            Nr = int(Re(r))**2 + int(Im(r))**2
            while Nr % 2 == 0: Nr //= 2
            assert Nr == 1, (names[i1], names[j1], r)
    for nm, F in atoms.items():
        assert expand(F.subs({p: 0, q: 1}))**2 == 1 and expand(F.subs({p: 1, q: 0}))**2 == 1
    say('A11 the resultants of m, n, W, W\' in pairs are powers of 2 up to units, and each form is a unit at (0,1), (1,0)')
    return True



# ---------------------------------------------------------------- quadratic orders Z[w], w^2 = D, and their completions
class QO:
    """elements a + b w of Z[w] (w^2 = D, D = 2 or -1) or of Z (D = None); exact integers."""
    def __init__(s, D): s.D = D
    def mul(s, x, y):
        if s.D is None: return (x[0]*y[0], 0)
        return (x[0]*y[0] + s.D*x[1]*y[1], x[0]*y[1] + x[1]*y[0])
    def add(s, x, y): return (x[0] + y[0], x[1] + y[1])
    def sub(s, x, y): return (x[0] - y[0], x[1] - y[1])
    def neg(s, x): return (-x[0], -x[1])
    def norm(s, x): return x[0] if s.D is None else x[0]*x[0] - s.D*x[1]*x[1]
    def pw(s, x, k):
        r = (1, 0)
        for _ in range(k): r = s.mul(r, x)
        return r
    def is_zero(s, x): return x[0] == 0 and x[1] == 0

def v_int(a, p):
    if a == 0: return 10**9
    k = 0
    while a % p == 0: a //= p; k += 1
    return k

class Local:
    """the completion of Q (D None), Q(sqrt2) (D = 2) or Q(i) (D = -1) at a prime; kind:
       ('rat', p)       Q_p;   ('ram', None) the prime over 2 (sqrt2, or 1+i);   ('inert', p) an inert odd prime p."""
    def __init__(s, R, kind, p):
        s.R, s.kind, s.p = R, kind, p
        if kind == 'rat':
            s.pi = (p, 0); s.res = [(r, 0) for r in range(p)]; s.prec = 3 if p == 2 else 1
            s.q = p
        elif kind == 'ram':
            s.pi = (0, 1) if R.D == 2 else (1, 1); s.res = [(0, 0), (1, 0)]; s.prec = 5; s.q = 2
        else:
            s.pi = (p, 0); s.res = [(a, b) for a in range(p) for b in range(p)]; s.prec = 1; s.q = p*p
        if s.prec > 1:   # unit squares modulo pi^prec, by enumeration
            s.sq = set()
            for a in range(16):
                for b in range(16 if R.D is not None else 1):
                    x = (a, b)
                    if s.val(x) != 0: continue
                    s.sq.add(s.reduce(R.mul(x, x)))
    def val(s, x):
        R = s.R
        if R.is_zero(x): return 10**9
        if s.kind == 'rat': return v_int(x[0], s.p)
        if s.kind == 'ram': return v_int(abs(R.norm(x)), 2)
        return min(v_int(x[0], s.p), v_int(x[1], s.p))
    def divpi(s, x):
        R = s.R
        if s.kind in ('rat', 'inert'):
            assert x[0] % s.p == 0 and x[1] % s.p == 0
            return (x[0] // s.p, x[1] // s.p)
        if R.D == 2:   # (a + b w)/w = b + (a/2) w
            assert x[0] % 2 == 0
            return (x[1], x[0] // 2)
        # (a + b i)/(1 + i) = ((a + b) + (b - a) i)/2
        assert (x[0] + x[1]) % 2 == 0
        return ((x[0] + x[1]) // 2, (x[1] - x[0]) // 2)
    def reduce(s, x):
        """canonical residue of a unit modulo pi^prec (only used for prec > 1, i.e. at 2)"""
        if s.kind == 'rat': return (x[0] % 8, 0)
        # modulo 8 = (unit) pi^6 for both rings: keep a, b mod 8, then test differences by valuation instead
        return (x[0] % 8, x[1] % 8)
    def is_square(s, x):
        """x != 0 a square in the completion"""
        R = s.R
        v = s.val(x)
        if v % 2: return False
        u = x
        for _ in range(v): u = s.divpi(u)
        if s.prec == 1:
            # odd residue characteristic: u^((q-1)/2) = 1 modulo p
            e = (s.q - 1)//2
            r = (1, 0); b = (u[0] % s.p, u[1] % s.p)
            while e:
                if e & 1: r = R.mul(r, b); r = (r[0] % s.p, r[1] % s.p)
                b = R.mul(b, b); b = (b[0] % s.p, b[1] % s.p); e >>= 1
            return r == (1, 0)
        # residue characteristic 2: u is a square iff u - x^2 has valuation >= prec for some unit x
        for a in range(16):
            for b in range(16 if R.D is not None else 1):
                y = (a, b)
                if s.val(y) != 0: continue
                if s.val(R.sub(u, R.mul(y, y))) >= s.prec: return True
        return False

def binom(n, k):
    from math import comb
    return comb(n, k)

def locally_soluble(L, F, cap=40):
    """Is there (M, e) primitive over the completion with F(M, e) a square (zero allowed)?  F = [f0, ..., f4],
    F(M, e) = f0 M^4 + f1 M^3 e + f2 M^2 e^2 + f3 M e^3 + f4 e^4, coefficients in the order.  Exact search over
    residue classes: a class t0 + pi^k O is decided when every higher Taylor coefficient has valuation at least
    v(c0) + prec (then the whole class has the square class of c0)."""
    R = L.R
    # chart 1: (M, e) = (t, 1), t in O;   chart 2: (M, e) = (1, t), t in pi O
    charts = [([F[4], F[3], F[2], F[1], F[0]], 0), ([F[0], F[1], F[2], F[3], F[4]], 1)]   # poly in t, low degree first
    for poly, k0 in charts:
        stack = [((0, 0), k0, R.pw(L.pi, k0))]
        while stack:
            t0, k, pik = stack.pop()
            assert k <= cap, 'local search did not terminate'
            # Taylor coefficients at t0, scaled by pi^(k j)
            c = []
            for j in range(5):
                acc = (0, 0)
                for i in range(j, 5):
                    acc = R.add(acc, R.mul((poly[i][0]*binom(i, j), poly[i][1]*binom(i, j)), R.pw(t0, i - j)))
                c.append(R.mul(acc, R.pw(pik, j)))
            if R.is_zero(c[0]): return True
            v0 = L.val(c[0])
            # Hensel: v(f(t0)) > 2 v(f'(t0)) gives a root of f within the chart, so a solution with value 0
            fd = (0, 0)
            for i in range(1, 5):
                fd = R.add(fd, R.mul((poly[i][0]*i, poly[i][1]*i), R.pw(t0, i - 1)))
            if not R.is_zero(fd) and v0 > 2*L.val(fd): return True
            if all(L.val(cj) >= v0 + L.prec for cj in c[1:]):
                if L.is_square(c[0]): return True
                continue
            for r in L.res:
                stack.append((R.add(t0, R.mul(r, pik)), k + 1, R.mul(pik, L.pi)))
    return False

def real_soluble(d, a, bd):
    """N^2 = d M^4 + a M^2 e^2 + bd e^4 over R (real numbers d, a, bd)"""
    if d > 0 or bd > 0: return True
    # both <= 0: need z > 0 with d z^2 + a z + bd > 0
    if d == 0 and bd == 0: return a > 0
    return a > 0 and a*a > 4*d*bd


# ---------------------------------------------------------------- exact arithmetic in Q(sqrt2): pairs of Fractions
R2 = QO(2)
def fmul(x, y): return (x[0]*y[0] + 2*x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def fadd(x, y): return (x[0] + y[0], x[1] + y[1])
def fsub(x, y): return (x[0] - y[0], x[1] - y[1])
def finv(x):
    N = x[0]*x[0] - 2*x[1]*x[1]; assert N != 0
    return (Fr(x[0]) / N, Fr(-x[1]) / N)
def fdiv(x, y): return fmul(x, finv(y))
def fsign(x, emb):
    """sign of a + b sqrt2 under sqrt2 -> emb*sqrt2 (emb = +-1), exact"""
    a, b = x[0], emb*x[1]
    if a >= 0 and b >= 0: return 0 if (a == 0 and b == 0) else 1
    if a <= 0 and b <= 0: return -1
    # opposite signs: compare a^2 with 2 b^2
    d = a*a - 2*b*b
    return (1 if a > 0 else -1) if d > 0 else ((1 if b > 0 else -1) if d < 0 else 0)
def fint(x):
    assert Fr(x[0]).denominator == 1 and Fr(x[1]).denominator == 1, x
    return (int(x[0]), int(x[1]))
def fpow(x, k):
    r = (1, 0)
    for _ in range(k): r = fmul(r, x)
    return r

EPS = (1, 1); SQ2 = (0, 1)
def sclass(z):
    """the class of an S-unit z of Q(sqrt2) (S = {sqrt2, 3}) in F*/F*^2 as a vector (sign, eps, sqrt2, 3)"""
    z = (Fr(z[0]), Fr(z[1]))
    N = z[0]*z[0] - 2*z[1]*z[1]
    vs = v_int(abs(N.numerator), 2) - v_int(N.denominator, 2)
    v3 = (v_int(abs(N.numerator), 3) - v_int(N.denominator, 3))
    assert v3 % 2 == 0; v3 //= 2
    u = z
    for _ in range(abs(vs)): u = fdiv(u, SQ2) if vs > 0 else fmul(u, SQ2)
    for _ in range(abs(v3)): u = fdiv(u, (3, 0)) if v3 > 0 else fmul(u, (3, 0))
    Nu = u[0]*u[0] - 2*u[1]*u[1]
    assert abs(Nu.numerator) == 1 and Nu.denominator == 1, ('not an S-unit', z)
    e = 1 if Nu < 0 else 0
    if e: u = fdiv(u, EPS)
    # now u = +- eps^(2j): its sign at the first embedding decides
    sg = 0 if fsign(u, 1) > 0 else 1
    return (sg, e, vs % 2, v3 % 2)
def from_class(c):
    z = (1, 0)
    if c[0]: z = fmul(z, (-1, 0))
    if c[1]: z = fmul(z, EPS)
    if c[2]: z = fmul(z, SQ2)
    if c[3]: z = fmul(z, (3, 0))
    return fint(z)

def descent_sqrt2(A, B, name):
    """2-isogeny descent over Q(sqrt2) for y^2 = x(x^2 + A x + B), A, B in Z[sqrt2] with all primes of B in
    {sqrt2, 3}.  Returns the classes not excluded by local insolubility (real places, sqrt2, 3)."""
    Lr = Local(R2, 'ram', None); L3 = Local(R2, 'inert', 3)
    import itertools
    keep = []
    for c in itertools.product((0, 1), repeat=4):
        dl = from_class(c)
        Bd = fdiv(B, dl)
        if Fr(Bd[0]).denominator != 1 or Fr(Bd[1]).denominator != 1: continue   # dl does not divide B
        Bd = fint(Bd)
        Fq = [dl, (0, 0), A, (0, 0), Bd]
        why = None
        for emb in (1, -1):
            sd, sa, sb = fsign(dl, emb), fsign(A, emb), fsign(Bd, emb)
            # real solubility: exact test
            if not (sd > 0 or sb > 0):
                # need z > 0 with d z^2 + a z + b > 0: a > 0 and a^2 > 4 d b
                a2 = fmul(A, A); db4 = fmul((4, 0), fmul(dl, Bd))
                if not (sa > 0 and fsign(fsub(a2, db4), emb) > 0): why = 'real%+d' % emb; break
        if why is None and not locally_soluble(Lr, Fq): why = 'sqrt2'
        if why is None and not locally_soluble(L3, Fq): why = '3'
        if why is None: keep.append(c)
        say('   %s: delta = %-12s %s' % (name, '%d%+d*sqrt2' % dl, 'kept' if why is None else 'insoluble at ' + why))
    return keep


def count_points_sqrt2(A, B, p):
    """#E(k) for y^2 = x(x^2 + A x + B) over the residue fields of Z[sqrt2] above the odd prime p (good reduction
    assumed and checked).  p = +-1 mod 8: two primes, sqrt2 -> r, -r with r^2 = 2 mod p; else one prime, F_{p^2}."""
    out = []
    if p % 8 in (1, 7):
        r = next(r for r in range(p) if (r*r - 2) % p == 0)
        for rr in (r, p - r):
            a = (A[0] + A[1]*rr) % p; b = (B[0] + B[1]*rr) % p
            disc = (b*b*(a*a - 4*b)) % p
            assert disc != 0
            cnt = 1
            for x in range(p):
                f = (x*(x*x + a*x + b)) % p
                cnt += 1 if f == 0 else (2 if pow(f, (p - 1)//2, p) == 1 else 0)
            out.append(cnt)
    else:
        q = p*p
        def mul(x, y): return ((x[0]*y[0] + 2*x[1]*y[1]) % p, (x[0]*y[1] + x[1]*y[0]) % p)
        def pw(x, e):
            r = (1, 0)
            while e:
                if e & 1: r = mul(r, x)
                x = mul(x, x); e >>= 1
            return r
        a = (A[0] % p, A[1] % p); b = (B[0] % p, B[1] % p)
        a2 = mul(a, a); disc = mul(mul(b, b), ((a2[0] - 4*b[0]) % p, (a2[1] - 4*b[1]) % p))
        assert disc != (0, 0)
        cnt = 1
        for x0 in range(p):
            for x1 in range(p):
                x = (x0, x1)
                f = mul(x, (mul(x, x)[0] + mul(a, x)[0] + b[0], mul(x, x)[1] + mul(a, x)[1] + b[1]))
                f = (f[0] % p, f[1] % p)
                cnt += 1 if f == (0, 0) else (2 if pw(f, (q - 1)//2) == (1, 0) else 0)
        out.append(cnt)
    return out


def fsqrt(z):
    """square root in Q(sqrt2) of z (pair of Fractions), or None"""
    a, b = Fr(z[0]), Fr(z[1])
    if a == 0 and b == 0: return (Fr(0), Fr(0))
    N = a*a - 2*b*b
    if N < 0: return None
    def qsqrt(r):
        if r < 0: return None
        nu, de = isqrt(r.numerator), isqrt(r.denominator)
        return Fr(nu, de) if nu*nu == r.numerator and de*de == r.denominator else None
    n0 = qsqrt(N)
    if n0 is None: return None
    for c2 in ((a + n0)/2, (a - n0)/2):
        c = qsqrt(c2)
        if c is None: continue
        if c == 0:
            d2 = a / 2 if b == 0 else None
            d = qsqrt(d2) if d2 is not None else None
            if d is not None and (c*c + 2*d*d, 2*c*d) == (a, b): return (c, d)
            continue
        d = b / (2*c)
        if (c*c + 2*d*d, 2*c*d) == (a, b): return (c, d)
    return None

class CurveF:
    """y^2 = x^3 + a2 x^2 + a4 x over Q(sqrt2), points as pairs of field elements, None = O"""
    def __init__(s, a2, a4): s.a2, s.a4 = a2, a4
    def on(s, P):
        if P is None: return True
        x, y = P
        return fsub(fmul(y, y), fmul(x, fadd(fadd(fmul(x, x), fmul(s.a2, x)), s.a4))) == (0, 0)
    def add(s, P, Q):
        if P is None: return Q
        if Q is None: return P
        if P[0] == Q[0]:
            if fadd(P[1], Q[1]) == (0, 0): return None
            num = fadd(fadd(fmul((3, 0), fmul(P[0], P[0])), fmul((2, 0), fmul(s.a2, P[0]))), s.a4)
            l = fdiv(num, fmul((2, 0), P[1]))
        else:
            l = fdiv(fsub(Q[1], P[1]), fsub(Q[0], P[0]))
        x = fsub(fsub(fsub(fmul(l, l), s.a2), P[0]), Q[0])
        y = fsub(fmul(l, fsub(P[0], x)), P[1])
        x = (Fr(x[0]), Fr(x[1])); y = (Fr(y[0]), Fr(y[1]))
        return (x, y)

def part_C():
    """Lemma D: the rational points of D2: w^2 = 2(x^2 + 1) H(x, 1) are (+-1, +-4)."""
    from sympy import symbols, sqrt as ssqrt, expand, Poly, simplify, nsimplify
    xx, ss, uu = symbols('x s u')
    r2 = ssqrt(2)
    f = 2*(xx**2 + 1)*(xx**4 - 2*xx**3 + 2*xx**2 + 2*xx + 1)
    al, be = 1 + r2, 1 - r2
    eps = 1 + r2
    G = expand(expand(f.subs(xx, (al - be*ss)/(1 - ss))*(1 - ss)**6))
    Gs = expand(simplify(G))
    gt = expand(16*r2*eps**3*(uu + 1)*(uu**2 + 10*uu + 1))
    assert simplify(expand(Gs.subs(ss**2, uu*eps**2)) - 0) is not None
    # G(s) = 16 sqrt2 eps^3 (u + 1)(u^2 + 10u + 1) with u = (s/eps)^2: check as a polynomial identity in s
    chk = expand(Gs - gt.subs(uu, (ss/eps)**2))
    assert simplify(chk) == 0
    # the chain to the model of the descent: with c = 16 sqrt2 eps^3 = lam^2 d, lam = 4 eps, d = 2 + sqrt2, and
    # v = u + 1, X = c v, x = X/lam^2, (c Y/lam^3)^2 = x (x^2 + 8 d x - 8 d^2) whenever Y^2 = c (u+1)(u^2+10u+1)
    cc = 16*r2*eps**3; lam = 4*eps; dd = 2 + r2
    assert simplify(expand(cc - lam**2*dd)) == 0
    vv = uu + 1; Xc = cc*vv/lam**2
    lhs = (cc/lam**3)**2 * cc*(uu + 1)*(uu**2 + 10*uu + 1)
    rhs = Xc*(Xc**2 + 8*dd*Xc - 8*dd**2)
    assert simplify(expand(lhs - rhs)) == 0
    say('C1  (1-s)^6 f((a - b s)/(1 - s)) = 16 sqrt2 eps^3 (u+1)(u^2+10u+1), u = (s/eps)^2, a, b = 1 +- sqrt2,',
        'eps = 1 + sqrt2: the involution x -> (x+1)/(x-1) of D2 is s -> -s')
    # the descent on E1: y^2 = x(x^2 + 8d x - 8d^2), d = 2 + sqrt2 (= sqrt2 eps)
    d = (2, 1)
    A = fint(fmul((8, 0), d)); B = fint(fmul((-8, 0), fmul(d, d)))
    Ap = (-2*A[0], -2*A[1]); Bp = fint(fsub(fmul(A, A), fmul((4, 0), B)))
    k1 = descent_sqrt2(A, B, 'E1 ')
    k2 = descent_sqrt2(Ap, Bp, "E1'")
    assert sorted(k1) == sorted([(0, 0, 0, 0), sclass(B)]) and sorted(k2) == sorted([(0, 0, 0, 0), sclass(Bp)])
    # E1[2](F) = {O, T}: x^2 + A x + B has discriminant A^2 - 4B = Bp, of class (3): not a square; same for E1'
    assert fsqrt(Bp) is None and fsqrt(fsub(fmul(Ap, Ap), fmul((4, 0), Bp))) is None
    say('C2  2-isogeny descent over Q(sqrt2): im(alpha) = {1, B}, im(alpha\') = {1, B\'} (the images of the 2-torsion',
        'points); 2^rank = |im alpha| |im alpha\'| / 4 = 1: rank E1(Q(sqrt2)) = 0')
    g = 0
    for p in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
        for c in count_points_sqrt2(A, B, p): g = gcd(g, c)
    assert g == 10
    say('C3  #E1(k) over the residue fields of the primes above 5, 7, ..., 31 (good reduction): gcd 10;',
        'so #E1(Q(sqrt2)) divides 10')
    # the map D2 -> E1 and the images of the four rational points of D2
    EPSf = (Fr(1), Fr(1)); alf = (Fr(1), Fr(1)); bef = (Fr(1), Fr(-1))
    c = fmul((16, 0), fmul(SQ2, fpow(EPS, 3)))
    lam = fmul((4, 0), EPS)
    def to_E1(xD, w):
        s = fdiv(fsub(xD, alf), fsub(xD, bef))
        u = fdiv(fmul(s, s), fmul(EPSf, EPSf))
        Y = fmul(w, fpow(fsub((1, 0), s), 3))
        v = fadd(u, (1, 0))
        X, Yp = fmul(c, v), fmul(c, Y)
        return (fdiv(X, fmul(lam, lam)), fdiv(Yp, fpow(lam, 3)))
    E1 = CurveF((Fr(A[0]), Fr(A[1])), (Fr(B[0]), Fr(B[1])))
    pts = [to_E1((Fr(x0), Fr(0)), (Fr(w0), Fr(0))) for x0 in (1, -1) for w0 in (4, -4)]
    assert all(E1.on(P) for P in pts)
    # the group they generate
    grp = [None]
    frontier = list(pts) + [((Fr(0), Fr(0)), (Fr(0), Fr(0)))]
    while frontier:
        P = frontier.pop()
        if P in grp: continue
        grp.append(P)
        for Q in list(grp): frontier.append(E1.add(P, Q))
    assert len(grp) == 10
    say('C4  the images of (+-1, +-4) in D2(Q) and T generate a group of order 10, so E1(Q(sqrt2)) is that group')
    # back to D2: a point (x, y) of E1 comes from D2 iff u = lam^2 x / c - 1 is a square in F (s = +-eps sqrt u);
    # s = infinity would need the leading coefficient 16 sqrt2 eps^-3 of G to be a square: it is not (odd at sqrt2)
    assert fsqrt(fmul((16, 0), fmul(SQ2, finv(fpow(EPS, 3))))) is None
    xs = set()
    for P in grp:
        if P is None: continue
        u = fsub(fdiv(fmul(fmul(lam, lam), P[0]), c), (1, 0))
        r = fsqrt(u)
        if r is None: continue
        for sg in (1, -1):
            s = fmul((sg, 0), fmul(EPSf, r))
            if fsub((1, 0), s) == (0, 0): continue   # s = 1 is x = infinity
            xD = fdiv(fsub(alf, fmul(bef, s)), fsub((1, 0), s))
            xs.add((Fr(xD[0]), Fr(xD[1])))
    xs_rat = sorted(x[0] for x in xs if x[1] == 0)
    # x = infinity (s = 1): u = 1/eps^2, a square; check whether it is the x of a point of E1(F)
    u1 = fdiv((1, 0), fmul(EPS, EPS)); X1 = fdiv(fmul(c, fadd(u1, (1, 0))), fmul(lam, lam))
    inf_pt = any(P is not None and P[0] == X1 for P in grp)
    # x = infinity: D2 has points there over Q(sqrt2) (leading coefficient 2 = sqrt2^2) but none over Q
    assert inf_pt
    say('C5  the finite x-coordinates of D2(Q(sqrt2)) are', sorted([(str(a), str(b)) for a, b in xs]),
        '(and x = infinity, which is not a point over Q: the leading coefficient 2 is not a rational square)',
        '; the rational ones:', [str(x) for x in xs_rat])
    def fval(x0): return 2*(x0*x0 + 1)*(x0**4 - 2*x0**3 + 2*x0*x0 + 2*x0 + 1)
    def qsq(r): r = Fr(r); return r >= 0 and isqrt(r.numerator)**2 == r.numerator and isqrt(r.denominator)**2 == r.denominator
    xs_Q = [x0 for x0 in xs_rat if qsq(fval(x0))]
    say('    of these, f(x) is a rational square only at', [str(x) for x in xs_Q], '(f(0) = 2)')
    assert xs_Q == [-1, 1]
    hits = []
    for b0 in range(1, 801):
        for a0 in range(-800, 801):
            if gcd(a0, b0) == 1 and qsq(fval(Fr(a0, b0))): hits.append((a0, b0))
    assert hits == [(-1, 1), (1, 1)]
    say('C6  (control, not used) a search of x = a/b, |a|, b <= 800, finds only x = -1, 1')
    say('Lemma D: 2(m^2+n^2) H(m,n) is a square (m, n coprime, n != 0) only for m = +-n.')
    return True


def part_B():
    """Lemma B: x^4 + 14x^2y^2 + y^4 = z^2 with gcd(x, y) = 1 only if xy(x^2 - y^2) = 0.  Hence m^8 + 14m^4n^4 + n^8
    is a square only if mn(m^2 - n^2) = 0.  Over Q: E: Y^2 = X(X^2 + 14X + 1) has rank 0 and E(Q) = {O, (0,0), (1, +-4)};
    a solution with y != 0 gives the point (x^2/y^2, xz/y^3), whose abscissa must be 0 or 1."""
    Z = QO(None); L2 = Local(Z, 'rat', 2); L3 = Local(Z, 'rat', 3)
    a, b = 14, 1
    ap, bp = -2*a, a*a - 4*b          # E': Y^2 = X(X^2 - 28X + 192)
    assert bp == 192
    def classes(bb):
        out = []
        for sg in (1, -1):
            for dd in range(1, abs(bb) + 1):
                if bb % dd == 0 and all(dd % (k*k) for k in range(2, isqrt(dd) + 1)):
                    out.append(sg*dd)
        return out
    def test(aa, bb, name):
        kept = []
        for dl in classes(bb):
            Fq = [(dl, 0), (0, 0), (aa, 0), (0, 0), (bb // dl, 0)]
            why = None
            if not real_soluble(dl, aa, bb // dl): why = 'real'
            elif not locally_soluble(L2, Fq): why = '2'
            elif not locally_soluble(L3, Fq): why = '3'
            if why is None: kept.append(dl)
            say('   %s: delta = %4d %s' % (name, dl, 'kept' if why is None else 'insoluble at ' + why))
        return kept
    k1 = test(a, b, 'E ')
    k2 = test(ap, bp, "E'")
    assert k1 == [1]
    # the kept classes of E' are images of points: T' (192 = 3 mod squares) and (8, 16) (8 = 2 mod squares)
    assert sorted(k2) == [1, 2, 3, 6]
    assert 16**2 == 8*(64 - 28*8 + 192)
    say('B1  im(alpha) = {1}, im(alpha\') = {1, 2, 3, 6} (from T\' = (0,0) and (8, 16)): 2^rank = 1*4/4, rank 0')
    # torsion by Nagell-Lutz: integral points with y = 0 or y^2 | disc = b^2 (a^2 - 4b) = 192
    disc = b*b*(a*a - 4*b)
    tors = []
    for y in range(0, isqrt(disc) + 1):
        if y and disc % (y*y): continue
        c0 = -y*y
        for xv in range(-abs(c0) - 20, abs(c0) + 21):
            if xv**3 + a*xv*xv + b*xv + c0 == 0: tors.append((xv, y))
    assert sorted(tors) == [(0, 0), (1, 4)]
    say('B2  Nagell-Lutz: the torsion points are among O, (0,0), (1,+-4), so E(Q) = {O, (0,0), (1,+-4)}')
    say('Lemma B: m^8 + 14m^4n^4 + n^8 is a square only if mn(m^2 - n^2) = 0 (x = m^4/n^4 must be 0 or 1).')
    return True

def gauss_is_square(a, b):
    """is a + b i (integers) a square in Z[i] (equivalently in Q(i))"""
    if a == 0 and b == 0: return True
    N = a*a + b*b; c = isqrt(N)
    if c*c != N: return False
    for s2 in ((c + a)//2,):
        if (c + a) % 2: return False
        s = isqrt(s2); t2 = (c - a)//2; t = isqrt(t2)
        if s*s != s2 or t*t != t2: return False
        return (s*s - t*t, 2*s*t) == (a, b) or (s*s - t*t, -2*s*t) == (a, b)
    return False

def qsquare(r):
    r = Fr(r)
    return r > 0 and isqrt(r.numerator)**2 == r.numerator and isqrt(r.denominator)**2 == r.denominator

def certify_ratio(m, n, forge=None):
    """For the curve E+ at (m, n) (E- at (m, -n) is E+ at (m, -n)): build the model with sigma_exact.model, find P, R,
    Q, T on it, and check that R, Q, T are independent in E(Q)/2E(Q) by the maps alpha (x - e) and x - theta.
    Exact.  Returns the dimension of the span of the classes found (3 means rank >= 2)."""
    sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
    from sigma_exact import model, Curve
    M = model(m, n, 2); C = Curve(M['a2'], M['a4'], M['a6']); A, z0, S = M['A'], M['z0'], M['S']
    a2, a4, a6 = C.a2, C.a4, C.a6
    cub = lambda X: X**3 + a2*X*X + a4*X + a6
    e = Fr(A) / (Fr(4*m*n, S) - z0)
    assert cub(e) == 0
    # the quadratic factor and theta (must lie in Q(i))
    # X^3 + a2 X^2 + a4 X + a6 = (X - e)(X^2 + uX + v)
    u = a2 + e; v = a4 + e*u
    assert cub(Fr(0)) == -e*v
    D = u*u - 4*v
    assert D < 0
    sD = Fr(isqrt((-D).numerator), isqrt((-D).denominator)); assert sD*sD == -D
    th = (-u/2, sD/2)                                   # theta = -u/2 + i sqrt(-D)/2
    # P (t = 1, z = 2) and R
    XP = Fr(A) / (2 - z0); y2 = cub(XP); assert qsquare(y2)
    P = (XP, Fr(isqrt(y2.numerator), isqrt(y2.denominator)))
    xR = Fr(-4*m*n*S*(m*m - n*n)); yR2 = cub(xR); assert qsquare(yR2)
    R = (xR, Fr(isqrt(yR2.numerator), isqrt(yR2.denominator)))
    twoR = C.add(R, R)
    assert twoR[0] == P[0]
    pts = {'R': R}
    if forge == 'R=P': pts['R'] = P
    # Q: X = 0, needs a6 a square
    if forge == 'Q-square-type' or qsquare(Fr(a6)):
        if not qsquare(Fr(a6)):
            raise ValueError('a6 is not a square: there is no point with X = 0')
        pts['Q'] = (Fr(0), Fr(isqrt(a6)))
    T = (e, Fr(0))
    fpe = 3*e*e + 2*a2*e + a4
    # classes: alpha in Q*/Q*^2 as the squarefree part; beta = x - theta in Q(i)*/Q(i)*^2 tested by products
    def alpha(Pt):
        if Pt == T: return fpe
        return Pt[0] - e
    def beta(Pt):
        if Pt == T: return (e + u/2, -sD/2)                 # e - theta
        return (Pt[0] + u/2, -sD/2)
    def gmul(x, y): return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])
    names = list(pts) + ['T']
    allp = dict(pts); allp['T'] = T
    use = dict((k, k) for k in names)
    if forge == 'T=R': use['T'] = 'R'          # a forged dependence: must be detected
    if forge == 'Q=R' and 'Q' in use: use['Q'] = 'R'
    import itertools
    dim_ok = True
    for mask in itertools.product((0, 1), repeat=len(names)):
        if not any(mask): continue
        al = Fr(1); be = (Fr(1), Fr(0))
        for nm, bit in zip(names, mask):
            if bit:
                al *= alpha(allp[use[nm]]); be = gmul(be, beta(allp[use[nm]]))
        # a product of classes is trivial in E/2E only if alpha is a rational square and beta a square in Q(i)
        L = 1
        for q in (be[0].denominator, be[1].denominator): L = L * q // gcd(L, q)
        bi = (int(be[0]*L*L), int(be[1]*L*L))
        if qsquare(al) and gauss_is_square(*bi): dim_ok = False
    return len(names) if dim_ok else None


def twice_square_ratios(M):
    out = []
    for m in range(2, M + 1):
        for n in range(1, m):
            if gcd(m, n) == 1:
                S = m*m + n*n
                if S % 2 == 0 and isqrt(S//2)**2 == S//2: out.append((m, n))
    return out

def part_D(nrandom=200, seed=20261124):
    """the conclusion checked curve by curve, from the model of sigma_exact.py (independent of Part A's formulas)"""
    tw = twice_square_ratios(1000)
    assert len(tw) == 127
    sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
    from sigma_exact import model
    for (m, n) in tw:
        Mm, Mp = model(m, n, -2), model(m, -n, 2)   # E- at (m, n), built directly, is E+ at (m, -n)
        assert (Mm['a2'], Mm['a4'], Mm['a6']) == (Mp['a2'], Mp['a4'], Mp['a6'])
        for nn in (n, -n):                      # E+ at (m, n) and E- at (m, n) = E+ at (m, -n)
            assert certify_ratio(m, nn) == 3, (m, nn)
    say('D1  the 127 ratios m <= 1000 with m^2+n^2 twice a square (those of sigma_twice.py), both curves',
        '(model(m, n, -2) = model(m, -n, 2) checked on each):',
        'R, Q, T independent in E(Q)/2E(Q) on all 254 (rank >= 2, agreeing with sigma_twice.py)')
    # random twice-square ratios with m up to 10^6, from m = p^2 + 2pq - q^2, n = |-p^2 + 2pq + q^2|
    rnd = random.Random(seed); seen = set(); cnt = 0
    while cnt < nrandom:
        p_ = rnd.randrange(2, 1000); q_ = rnd.randrange(1, p_)
        if gcd(p_, q_) != 1 or (p_ + q_) % 2 == 0: continue
        m_, n_ = p_*p_ + 2*p_*q_ - q_*q_, abs(-p_*p_ + 2*p_*q_ + q_*q_)
        if m_ < n_: m_, n_ = n_, m_
        if n_ == 0 or m_ > 10**6 or (m_, n_) in seen: continue
        seen.add((m_, n_)); cnt += 1
        S = m_*m_ + n_*n_; assert S % 2 == 0 and isqrt(S//2)**2 == S//2 and gcd(m_, n_) == 1
        for nn in (n_, -n_):
            assert certify_ratio(m_, nn) == 3, (m_, nn)
    say('D2  %d random twice-square ratios with m <= 10^6 (largest m %d), both curves: R, Q, T independent'
        % (cnt, max(x[0] for x in seen)))
    # the square type: no Q; R and T are independent (rank >= 1), the theorem claims nothing more
    sq = [(m, n) for m in range(2, 200) for n in range(1, m) if gcd(m, n) == 1 and isqrt(m*m + n*n)**2 == m*m + n*n]
    for (m, n) in sq:
        for nn in (n, -n):
            assert certify_ratio(m, nn) == 2, (m, nn)
    say('D3  the %d square-type ratios with m < 200, both curves: only R and T (no rational point with X = 0);'
        ' they are independent in E(Q)/2E(Q): rank >= 1 (Theorem 1), nothing more is claimed' % len(sq))
    # other ratios: R, T independent as Theorem 1 says
    oth = [(m, n) for m in range(2, 60) for n in range(1, m) if gcd(m, n) == 1]
    for (m, n) in oth:
        for nn in (n, -n):
            assert certify_ratio(m, nn) in (2, 3)
    say('D4  all %d coprime ratios with m < 60, both curves: R, T independent in E(Q)/2E(Q)' % len(oth))
    return True

def part_D_pari():
    """guide only (floating point inside PARI's ellrank): ranks of a few curves of each type"""
    try:
        import cypari2
    except ImportError:
        say('D5  (cypari2 not available: skipped)'); return
    pari = cypari2.Pari(); pari.allocatemem(10**9)
    sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
    from sigma_exact import model
    rows = []
    for (m, n) in [(3, 2), (4, 1), (4, 3), (5, 2), (8, 3), (12, 5), (15, 8), (24, 7), (21, 20), (35, 12),
                   (7, 1), (17, 7), (23, 7), (41, 1)]:
        for nn in (n, -n):
            M = model(m, nn, 2)
            E = pari.ellinit([0, M['a2'], 0, M['a4'], M['a6']])
            S = m*m + n*n
            xR = -4*m*nn*S*(m*m - n*n)
            pts = [[xR, pari.ellordinate(E, xR)[0]]]
            if isqrt(M['a6'])**2 == M['a6']: pts.append([0, isqrt(M['a6'])])
            r = pari.ellrank(E, 0, pts)
            typ = 'twice' if (S % 2 == 0 and isqrt(S//2)**2 == S//2) else ('square' if isqrt(S)**2 == S else 'other')
            rows.append('(%d,%d) %s [%d,%d]' % (m, nn, typ, int(r[0]), int(r[1])))
            assert int(r[1]) >= (2 if typ == 'twice' else 1) and int(r[0]) >= (2 if typ == 'twice' else 1)
    say('D5  PARI ellrank [lower, upper] with R (and Q) supplied (guide only; consistent with the bounds proved here):',
        '; '.join(rows))
    # Lemma B's curve and the curve of Lemma D, by PARI (independent of the code above)
    rB = pari.ellrank(pari.ellinit([0, 14, 0, 1, 0]))
    assert int(rB[0]) == 0 and int(rB[1]) == 0
    # Frobenius of J(D2) (PARI hyperellcharpoly) against the point counts of E1 over the residue fields of Q(sqrt2):
    # J(D2) is isogenous to the restriction of scalars of E1, so at p split the charpoly is the product of the two
    # (x^2 - a x + p), and at p inert it is P(x^2) with P the charpoly of Frobenius of E1 over F_(p^2)
    d = (2, 1); A = fint(fmul((8, 0), d)); B = fint(fmul((-8, 0), fmul(d, d)))
    fD = '2*(x^2+1)*(x^4-2*x^3+2*x^2+2*x+1)'
    for pp in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
        cp = pari('hyperellcharpoly(Mod(1,%d)*(%s))' % (pp, fD))
        cnt = count_points_sqrt2(A, B, pp)
        if len(cnt) == 2:
            pr = pari('(x^2 - (%d)*x + %d)*(x^2 - (%d)*x + %d)' % (pp + 1 - cnt[0], pp, pp + 1 - cnt[1], pp))
        else:
            q = pp*pp; a = q + 1 - cnt[0]
            pr = pari('x^4 - (%d)*x^2 + %d' % (a, q))
        assert cp == pr, (pp, cp, pr)
    # the Q(sqrt2) descent on the conjugate of E1 (must also give rank 0) and on E0: y^2 = v(v^2 + 8v - 8), whose rank
    # over Q(sqrt2) is rank E0(Q) + rank E0^(2)(Q), computed by PARI
    dc = (2, -1); Ac = fint(fmul((8, 0), dc)); Bc = fint(fmul((-8, 0), fmul(dc, dc)))
    kc1 = descent_sqrt2(Ac, Bc, 'conj E1 '); kc2 = descent_sqrt2((-2*Ac[0], -2*Ac[1]), fint(fsub(fmul(Ac, Ac), fmul((4, 0), Bc))), "conj E1'")
    assert len(kc1)*len(kc2) == 4
    k01 = descent_sqrt2((8, 0), (-8, 0), 'E0 '); k02 = descent_sqrt2((-16, 0), (96, 0), "E0'")
    r0 = int(pari.ellrank(pari.ellinit([0, 8, 0, -8, 0]))[1]) + int(pari.ellrank(pari.ellinit([0, 16, 0, -32, 0]))[1])
    assert (len(k01)*len(k02)).bit_length() - 1 - 2 == r0 == 1
    say('D7  the Q(sqrt2) descent gives rank 0 for the conjugate of E1 too, and the bound 1 for y^2 = v(v^2 + 8v - 8) over',
        'Q(sqrt2), which is rank 1 + 0 of that curve and its twist by 2 over Q by PARI')
    say('D6  PARI (guide, independent of this script): rank of y^2 = x^3 + 14x^2 + x is [0, 0]; the Frobenius',
        'polynomials of J(D2) at 5, 7, ..., 31 (hyperellcharpoly) are those of the restriction of scalars of E1 computed',
        'from the point counts of Part C, so E1 is the quotient it should be')

def part_E():
    """negative tests: each forged input must be refused or must make the independence test fail"""
    n_ok = 0
    # 1. a forged point Q on a square-type curve is refused
    for (m, n) in [(4, 3), (12, 5), (15, 8), (24, 7)]:
        try:
            certify_ratio(m, n, forge='Q-square-type'); raise AssertionError('forged Q accepted')
        except ValueError:
            n_ok += 1
    # 2. R replaced by P = 2R (trivial in E/2E), T's class replaced by R's, Q's class replaced by R's
    for (m, n) in twice_square_ratios(200):
        for f in ('R=P', 'T=R', 'Q=R'):
            assert certify_ratio(m, n, forge=f) is None, (m, n, f); n_ok += 1
    # 3. the local solubility tests keep the classes of actual points: y^2 = x^3 - 6x over Q(sqrt2), the point
    #    (-2, 2) of infinite order has alpha = -2 = -1 mod squares; and over Q, y^2 = x^3 - 25x, (-4, 6): alpha = -1
    kept = descent_sqrt2((0, 0), (-6, 0), 'test')
    assert sclass((-2, 0)) in kept; n_ok += 1
    Z = QO(None)
    for L in (Local(Z, 'rat', 2), Local(Z, 'rat', 3), Local(Z, 'rat', 5)):
        assert locally_soluble(L, [(-1, 0), (0, 0), (0, 0), (0, 0), (25, 0)]); n_ok += 1
    # 4. a forged descent claim: pretending E1 had the class sqrt2 in its image must contradict the local test
    d = (2, 1); A = fint(fmul((8, 0), d)); B = fint(fmul((-8, 0), fmul(d, d)))
    Lr = Local(R2, 'ram', None)
    Bd = fint(fdiv(B, SQ2))
    assert not locally_soluble(Lr, [SQ2, (0, 0), A, (0, 0), Bd]); n_ok += 1
    say('E   negative tests: %d forged inputs, all refused' % n_ok)
    engine_crosscheck()
    return n_ok

def brute_local(L, F, B):
    """independent check of locally_soluble: try (M, e) with small integer coordinates (|a|, |b| <= B for a + b w),
    not both divisible by pi, and test the exact value F(M, e) with is_square; a hit proves solubility"""
    R = L.R
    rng = range(-B, B + 1)
    elts = [(a, b) for a in rng for b in (rng if R.D is not None else [0])]
    for M in elts:
        for e in elts:
            if L.val(M) > 0 and L.val(e) > 0: continue
            v = (0, 0)
            for i, c in enumerate(F):
                v = R.add(v, R.mul(c, R.mul(R.pw(M, 4 - i), R.pw(e, i))))
            if R.is_zero(v) or L.is_square(v): return True
    return None

def engine_crosscheck(seed=1):
    """random quartics over Q_2, Q_3, Q_5 and over Q(sqrt2) at sqrt2 and 3: whenever the brute force finds a point
    the engine must say soluble, and the two must agree on every sample (a brute force cannot prove insolubility,
    so agreement there is evidence, not proof)"""
    rnd = random.Random(seed); Z = QO(None)
    sol = ins = 0
    cases = [(Local(Z, 'rat', 2), 60, 0), (Local(Z, 'rat', 3), 60, 0), (Local(Z, 'rat', 5), 60, 0),
             (Local(R2, 'ram', None), 5, 1), (Local(R2, 'inert', 3), 5, 1)]
    for L, K, irr in cases:
        for k in range(40):
            rc = lambda h: (rnd.randint(-h, h), rnd.randint(-h, h) if irr else 0)
            if k % 2:   # a general quartic
                F = [rc(12) for _ in range(5)]
            else:       # a homogeneous space of a 2-isogeny descent, most of them insoluble somewhere
                F = [rc(12), (0, 0), rc(30), (0, 0), rc(40)]
            if R2.is_zero(F[0]) or R2.is_zero(F[4]): continue
            a = locally_soluble(L, F); b = brute_local(L, F, K)
            assert (b is True) == a, (L.kind, L.p, F, a, b)
            sol += a; ins += (not a)
    say('E2  local solubility engine against a brute-force search on %d random quartics (%d soluble, %d not):'
        ' agreement on all' % (sol + ins, sol, ins))


def part_A_dup():
    """the duplication identity used for Q: on y^2 = (x - e1)(x - e2)(x - e3),
    x(2P) - e1 = ((x - e1)^2 - (e1 - e2)(e1 - e3))^2 / (4 y^2), so x(2P) - e1 is a square in Q(e1)"""
    from sympy import symbols, expand, cancel
    x, e1, e2, e3 = symbols('x e1 e2 e3')
    f = (x - e1)*(x - e2)*(x - e3)
    fp = f.diff(x)
    lam2 = fp**2/(4*f)
    s1 = e1 + e2 + e3
    x2 = lam2 - (-s1) - 2*x          # x(2P) = lambda^2 - a2 - 2x with a2 = -(e1+e2+e3)
    assert cancel(x2 - e1 - ((x - e1)**2 - (e1 - e2)*(e1 - e3))**2/(4*f)) == 0
    say('A12 x(2P) - e1 = ((x - e1)^2 - (e1 - e2)(e1 - e3))^2 / (4y^2) on y^2 = (x - e1)(x - e2)(x - e3)')


def main():
    if '--negative' in sys.argv:
        part_E(); return
    say('sigma_generic.py: uniform rank bounds for the curves E+, E- of F0 F8 GL (UNIFORM_RANK.md); exact arithmetic')
    say('-- Part A: identities in general m, n (sympy)')
    part_A(); part_A_dup()
    say('-- Part B: Lemma B (descent over Q)')
    part_B()
    say('-- Part C: Lemma D (descent over Q(sqrt2))')
    part_C()
    say('-- Part D: the conclusion checked curve by curve')
    part_D()
    if '--pari' in sys.argv: part_D_pari()
    say('-- Part E: negative tests')
    part_E()
    say('sigma_generic OK')
    out = __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),
                                      'sigma_generic.txt')
    if '--write' in sys.argv:
        open(out, 'w').write('\n'.join(LOG) + '\n')

if __name__ == '__main__':
    main()
