#!/usr/bin/env python3
"""L(E,1) > 0 for the curve of (440,279) and (719,161), with the conductor, the signs a_p and the root number computed
here in exact arithmetic, not read from PARI (item B1 of AGENT_TASKS.md; written by an AI model, not read by a person).

What lball.py takes from PARI (its line 22: ellminimalmodel, ellglobalred, ellrootno, ellan), this program computes itself:

 1. Conductor.  Tate's algorithm (Silverman, Advanced Topics, IV.9.4; the order of the steps and the coordinate changes
    follow Cremona's implementation as in Sage's ell_local_data._tate) in integer arithmetic at every prime p | Delta,
    2 and 3 included: a p-minimal model, the Kodaira type, f_p, the Tamagawa number c_p, and for multiplicative p whether
    it is split.  Every division is asserted exact.  f_p is computed by Ogg's formula f_p = v(Delta) + 1 - (number of
    components) (Saito at 2) from the Kodaira type, as Tate's algorithm does, so asserting that formula would be true by
    construction (an earlier version did); instead the program asserts f_p = 0 at good and 1 at multiplicative p, and
    at additive p the bounds 2 <= f_p <= 2 (p >= 5), 5 (p = 3), 8 (p = 2) (Lockhart-Rosen-Silverman 1993, Brumer-
    Kramer 1994).  The conductor itself is compared with PARI by --tatetest and --census, and below.  The prime factors of Delta are found by trial division (sympy.factorint for a cofactor above
    10^6, as a finder only), the product is checked, and each factor is proved prime by trial division (for this curve
    all are below 10^3).
 2. a_n.  a_p = p + 1 - #E(F_p) at good p, by counting points: #E(F_p) = p + 1 + sum_x chi(4x^3 + b2 x^2 + 2 b4 x + b6),
    chi the Legendre symbol read off the table of squares mod p (numpy int64, all values below p^2 < 2^63: exact); for a
    model y^2 = x(x - e2)(x - e3), as here, a_p = -sum_x chi(x) chi(x - e2) chi(x - e3) from the same table (checked
    against the first count for every p < 3000).
    At bad p: a_p = +1 split, -1 non-split, 0 additive, from step 1.  a_n by multiplicativity and
    a_{p^(k+1)} = a_p a_{p^k} - p a_{p^(k-1)} (good p), a_{p^k} = a_p^k (bad p).  Hasse |a_p| <= 2 sqrt(p) is asserted.
 3. Root number, two ways.
    (a) Local: w = -prod_p w_p, w_p = -1 split, +1 non-split, additive p >= 5 by Rohrlich (Compositio 87 (1993),
        Prop. 2; as restated in Rizzo, Compositio 136 (2003) 1-23, section 1, item (5)), additive potentially
        multiplicative at 2 and 3 by Rohrlich as restated in Rizzo's item (4) ((-1/3) at 3; at 2, w_2 = -c6' mod 4),
        and in potentially good reduction at p = 2 one row of Rizzo's Table III (it completes Halberstadt, C. R. Acad.
        Sci. 326 (1998) 1047-1052): the row v_2(c4, c6, Delta) = (4, >=7, 6) with c4' = 3 mod 4, Kodaira III,
        f_2 = 5, where w_2 = +1 iff c4' - 4 c_{6,7} = 7 or 11 mod 16 (x' is x without its powers of 2; c_{6,7} is read
        as c6/2^7, integral in this row).  The table was read only through a text extraction of the PDF (a direct
        download was refused by the network proxy), twice: the first reading left the sign illegible, the second shows
        "-".  --rowtest: "-" agrees with PARI on all 3000 test curves of the row, "+" fails on all 678 with v2(c6) = 7;
        the 1166 with v2(c6) = 8, where 4 c6/2^7 = 8 mod 16, also agree, which supports the normalization c6/2^7.
        For this curve v_2(c6) = 10, so 4 c6/2^7 = 0 mod 16 and neither the sign nor the normalization matters: the
        condition is c4' = 7 or 11 mod 16, and c4' = 3 mod 16 gives w_2 = -1.
        So the sign of the c6 term in this row is PARI-calibrated (chosen because it agrees with PARI), not derived
        from Rizzo's paper: where that term is non-zero mod 16 (v_2(c6) = 7) this route is not independent of PARI.
        No curve used depends on it: this curve has v_2(c6) = 10, and the 50 closing curves of --census that fall in
        this row all have v_2(c6) >= 9 (checked on 26 November 2026 by a separate script, not committed).
        Any other additive case at 2 or 3 is reported as "w_p not implemented (would trust PARI)".
    (b) Global, from the functional equation itself (modularity): if f = sum a_n q^n is the newform of level N (Wiles,
        Breuil-Conrad-Diamond-Taylor; level = conductor, Carayol), then f(i/(y sqrt N)) = w y^2 f(i y/sqrt N) for y > 0;
        at y = 1 this is f(i/sqrt N) = w f(i/sqrt N).  So if the enclosure of F = sum a_n exp(-2 pi n/sqrt N) excludes 0,
        w = +1.  This uses only N (step 1) and the a_n (step 2).
 4. L(E,1).  For y > 0, L(E,1) = sum a_n/n (exp(-2 pi n y/sqrt N) + w exp(-2 pi n/(y sqrt N))) (split the Mellin integral
    at y and apply 3(b)).  It is enclosed at y = 1 (the sum of lball.py) and at y = 6/5, the "two cut-offs".  Arithmetic:
    fixed point with 2^-PREC, every quantity an integer interval [lo, hi] rounded outward: pi by Machin's formula with its
    truncation and rounding bounded, sqrt N by isqrt, exp(-t) by its alternating Taylor series (partial sums bracket it),
    q^n by repeated multiplication rounded down / up, each term a_n q^n/n rounded down / up.  Tail: |a_n| <= d(n) sqrt n
    <= 2n (Hasse and the recursion), so sum_{n>M} |a_n|/n q^n <= 2 q^(M+1)/(1-q) and sum_{n>M} |a_n| q^n
    <= 2 q^(M+1)((M+1) - M q)/(1-q)^2.  No floating point enters a certified bound (floats only choose M and print).
 5. Guard.  The two enclosures of L(E,1) must intersect and the Fricke residual f(i/(y sqrt N)) - w y^2 f(i y/sqrt N) at
    y = 6/5 must contain 0.  --negative runs evaluate() itself (both tests) with a forged conductor (N/2, 2N, N/7,
    3N) or sign (w = -1), tails below 10^-3, and with one forged a_p (a_2, a_3, a_7, a_521, a_719), tails below 10^-9,
    and requires that every forgery is refused.  The guard is a check, not part of the proof: the proof is N (1), the a_n
    (2), w = +1 (3b, and 3a), and the lower end of the enclosure at y = 1 (4).

L(E,1) > 0 gives rank 0 by modularity and then Kolyvagin (with Gross-Zagier and a nonvanishing twist: Bump-Friedberg-
Hoffstein, Murty-Murty), or by Kato's theorem, which also needs modularity.
Derivation used in 3(b) and 4: Lambda(s) = N^(s/2) (2 pi)^(-s) Gamma(s) L(E,s) = int_0^oo f(i t/sqrt N) t^(s-1) dt; the
Fricke involution gives f(-1/(N tau)) = eps N tau^2 f(tau) with w = -eps, i.e. f(i/(t sqrt N)) = w t^2 f(i t/sqrt N);
splitting the integral at t = y and substituting t -> 1/t on (0, y) gives at s = 1 the formula of 4.  Trusted: those theorems (with Atkin-Lehner: f is an eigenform of the Fricke involution W_N with eigenvalue eps = +-1,
and w = -eps), Tate's algorithm and Rohrlich's and Rizzo's formulas as published (the local product 3(a) is a
corroboration of 3(b), which alone is used), Python's
integers and numpy's int64.  PARI is used only for the cross-checks printed in their own section (and by --rowtest,
--tatetest); with --nopari it is not imported at all.

usage: lball_exact.py [--write] [--nopari]   the curve of (440,279); output lball_exact.txt (written only with --write)
       lball_exact.py --negative [--write]   the forged inputs; output lball_exact_negative.txt
       lball_exact.py --all [--write]        both, with one table of a_p (the two outputs one after the other)
       lball_exact.py --local a1,a2,a3,a4,a6  Tate's algorithm and the local root numbers only, for any curve
       lball_exact.py --census second_123.jsonl [--write]  N and w exactly for all 89 closing curves, against PARI
       lball_exact.py --tatetest K [--write] Tate's algorithm and w_p against PARI on K random models (a test; lball_exact_tatetest.txt)
       lball_exact.py --rowtest K [--write]  the row of Rizzo's Table III against PARI's ellrootno(E, 2) (lball_exact_rowtest.txt)
CI (proposed, job quick):  python3 lball_exact.py --all | diff - <(cat lball_exact.txt lball_exact_negative.txt)"""
import sys, math, time

PREC = 192                                   # bits of the fixed-point grid
ONE = 1 << PREC
CURVE = (0, -491040, 0, -13400146081, 0)     # y^2 = x(x - 719^2)(x + 161^2), the twist by -1 of the Frey curve of (440,279)
LABEL = '(440,279) and its partner (719,161): y^2 = x(x - 719^2)(x + 161^2)'
Y2 = (6, 5)                                  # the second cut-off y = 6/5
TARGET = 10**-9                              # tail of each L-value sum below this


# ---------------------------------------------------------------- integers
def val(x, p):
    assert x != 0
    k = 0
    while x % p == 0: x //= p; k += 1
    return k

def is_prime_td(n):
    """primality by trial division (a proof); only used for n < 10^14"""
    if n < 2: return False
    if n % 2 == 0: return n == 2
    r = math.isqrt(n); d = 3
    while d <= r:
        if n % d == 0: return False
        d += 2
    return True

def factor_certified(n):
    """{p: e} with prod p^e = |n|, each p proved prime by trial division.  The factors are found by trial division to 10^6
    and, for what is left, by sympy.factorint (a finder only: the product and the primality are checked here)."""
    n = abs(n); m = n; f = {}; d = 2
    while d * d <= m and d < 10**6:
        while m % d == 0: f[d] = f.get(d, 0) + 1; m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        import sympy
        for p, e in sympy.factorint(m).items(): f[int(p)] = f.get(int(p), 0) + int(e)
    assert math.prod(p**e for p, e in f.items()) == n
    for p in f:
        if not (p < 10**14 and is_prime_td(p)): raise ValueError('factor %d not proved prime here' % p)
    return f

def legendre(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


# ---------------------------------------------------------------- Weierstrass models (exact)
def b_inv(a):
    a1, a2, a3, a4, a6 = a
    b2 = a1*a1 + 4*a2; b4 = 2*a4 + a1*a3; b6 = a3*a3 + 4*a6
    b8 = a1*a1*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3*a3 - a4*a4
    return b2, b4, b6, b8

def c_inv(a):
    b2, b4, b6, b8 = b_inv(a)
    c4 = b2*b2 - 24*b4; c6 = -b2**3 + 36*b2*b4 - 216*b6
    D = -b2*b2*b8 - 8*b4**3 - 27*b6*b6 + 9*b2*b4*b6
    assert c4**3 - c6**2 == 1728*D and 4*b8 == b2*b6 - b4*b4
    return c4, c6, D

def rst(a, r, s, t):
    """x = x' + r, y = y' + s x' + t (u = 1); Silverman, Table III.1.2"""
    a1, a2, a3, a4, a6 = a
    n = (a1 + 2*s, a2 - s*a1 + 3*r - s*s, a3 + r*a1 + 2*t,
         a4 - s*a3 + 2*r*a2 - (t + r*s)*a1 + 3*r*r - 2*s*t,
         a6 + r*a4 + r*r*a2 + r**3 - t*a3 - t*t - r*t*a1)
    assert c_inv(n)[2] == c_inv(a)[2]
    return n

def unscale(a, p):
    """a_i -> a_i / p^i (the model x = p^2 x', y = p^3 y'); every division exact"""
    out = []
    for ai, i in zip(a, (1, 2, 3, 4, 6)):
        assert ai % p**i == 0
        out.append(ai // p**i)
    return tuple(out)

def ediv(x, d):
    assert x % d == 0, 'inexact division'
    return x // d


# ---------------------------------------------------------------- Tate's algorithm at one prime
def _roots_quad(A, B, C, p):
    """does A T^2 + B T + C have a root in F_p"""
    A, B, C = A % p, B % p, C % p
    if p <= 3 or A == 0:
        return any((A*T*T + B*T + C) % p == 0 for T in range(p)) if p <= 3 else (B != 0 or C == 0)
    return legendre(B*B - 4*A*C, p) >= 0

def _polymulmod(f, g, m, p):
    """f*g mod (m monic), coefficient lists low to high, over F_p"""
    r = [0] * (len(f) + len(g) - 1)
    for i, x in enumerate(f):
        if x:
            for j, y in enumerate(g): r[i + j] = (r[i + j] + x*y) % p
    d = len(m) - 1
    for k in range(len(r) - 1, d - 1, -1):
        c = r[k]
        if c:
            for j in range(d + 1): r[k - d + j] = (r[k - d + j] - c*m[j]) % p
    r = r[:d] if len(r) > d else r + [0]*(d - len(r))
    return r

def _polygcd_deg(f, g, p):
    def trim(h):
        h = [x % p for x in h]
        while h and h[-1] == 0: h.pop()
        return h
    f, g = trim(f), trim(g)
    while g:
        inv = pow(g[-1], -1, p); f = f[:]
        while len(f) >= len(g):
            c = f[-1] * inv % p; s = len(f) - len(g)
            for j in range(len(g)): f[s + j] = (f[s + j] - c*g[j]) % p
            f = trim(f)
            if not f: break
        f, g = g, f
    return len(f) - 1

def _ncubicroots(b, c, d, p):
    """number of distinct roots in F_p of T^3 + b T^2 + c T + d"""
    if p < 1000: return sum(1 for T in range(p) if (T**3 + b*T*T + c*T + d) % p == 0)
    m = [d % p, c % p, b % p, 1]; x = [0, 1, 0]; res = [1, 0, 0]; e = p
    while e:
        if e & 1: res = _polymulmod(res, x, m, p)
        x = _polymulmod(x, x, m, p); e >>= 1
    h = res[:]; h[1] = (h[1] - 1) % p                  # x^p - x mod m
    return _polygcd_deg(m, h, p) if any(h) else 3

def tate(a, p):
    """Tate's algorithm at p.  Returns dict: model (p-minimal, integral), kod (string), f, c, red, ap, vD, split."""
    a = tuple(a)
    while True:
        c4, c6, D = c_inv(a)
        vD = val(D, p)
        if vD == 0:
            return dict(model=a, kod='I0', f=0, c=1, red='good', ap=None, vD=0)
        a1, a2, a3, a4, a6 = a
        b2, b4, b6, b8 = b_inv(a)
        half = pow(2, -1, p) if p != 2 else None
        # move the singular point to (0, 0)
        if p == 2:
            if b2 % 2 == 0:
                r = a4 % 2; t = (((r + a2)*r + a4)*r + a6) % 2
            else:
                r = a3 % 2; t = (a4 + r*r) % 2            # a1 is odd here, so 1/a1 = 1 mod 2
        elif p == 3:
            r = (-b6) % 3 if b2 % 3 == 0 else (-pow(b2, -1, 3) * b4) % 3
            t = (a1*r + a3) % 3
        else:
            r = (-pow(12, -1, p) * b2) % p if c4 % p == 0 else (-pow(12*c4, -1, p) * (c6 + b2*c4)) % p
            t = (-half * (a1*r + a3)) % p
        a = rst(a, r, 0, t); a1, a2, a3, a4, a6 = a
        b2, b4, b6, b8 = b_inv(a)
        assert a3 % p == 0 and a4 % p == 0 and a6 % p == 0
        trace = 'singular point moved to (0,0) by (r, s, t) = (%d, 0, %d): a = %s; v(b2) %s, v(a6) %s, v(b8) %s, v(b6) %s' % (
            r, t, list(a), *(val(z, p) if z else 'inf' for z in (b2, a6, b8, b6)))
        if b2 % p != 0:                                   # multiplicative
            split = _roots_quad(1, a1, -a2, p)
            cp = vD if split else (2 if vD % 2 == 0 else 1)
            return dict(model=a, kod='I%d' % vD, f=1, c=cp, red='split' if split else 'nonsplit', ap=1 if split else -1, vD=vD)
        # additive
        if a6 % p**2 != 0:
            return dict(model=a, kod='II', f=vD, c=1, red='additive', ap=0, vD=vD, trace=trace + '; p^2 does not divide a6: II')
        if b8 % p**3 != 0:
            return dict(model=a, kod='III', f=vD - 1, c=2, red='additive', ap=0, vD=vD, trace=trace + '; p^2 | a6, p^3 does not divide b8: III, f = v(Delta) - 1')
        if b6 % p**3 != 0:
            cp = 3 if _roots_quad(1, ediv(a3, p), -ediv(a6, p*p), p) else 1
            return dict(model=a, kod='IV', f=vD - 2, c=cp, red='additive', ap=0, vD=vD)
        # change coordinates so that p | a1, a2; p^2 | a3, a4; p^3 | a6
        if p == 2: s = a2 % 2; t = 2 * (ediv(a6, 4) % 2)
        elif p == 3: s = a1; t = a3
        else: s = -a1 * half; t = -a3 * half             # not reduced mod p: a3 + 2t = a3 (1 - 2 half) is divisible by p^2
        a = rst(a, 0, s, t); a1, a2, a3, a4, a6 = a
        assert a1 % p == 0 and a2 % p == 0 and a3 % p**2 == 0 and a4 % p**2 == 0 and a6 % p**3 == 0
        b = ediv(a2, p); c = ediv(a4, p*p); d = ediv(a6, p**3)
        w = 27*d*d - b*b*c*c + 4*b**3*d - 18*b*c*d + 4*c**3   # minus the discriminant of T^3 + bT^2 + cT + d
        x = 3*c - b*b
        sw = 1 if w % p != 0 else (2 if x % p != 0 else 3)
        if sw == 1:                                       # I0*
            return dict(model=a, kod='I0*', f=vD - 4, c=1 + _ncubicroots(b, c, d, p), red='additive', ap=0, vD=vD)
        if sw == 2:                                       # I_m*, m > 0: move the double root to 0
            if p == 2: r = c % 2
            elif p == 3: r = (c * pow(b, -1, 3)) % 3
            else: r = ((b*c - 9*d) * pow(2*x, -1, p)) % p
            a = rst(a, p*r, 0, 0); a1, a2, a3, a4, a6 = a
            ix = iy = 3; mx = my = p*p
            while True:
                a2t = ediv(a2, p); a3t = ediv(a3, my); a4t = ediv(ediv(a4, p), mx); a6t = ediv(ediv(a6, mx), my)
                if (a3t*a3t + 4*a6t) % p == 0:
                    t = my * (a6t % 2) if p == 2 else my * ((-a3t * half) % p)
                    a = rst(a, 0, 0, t); a1, a2, a3, a4, a6 = a
                    my *= p; iy += 1
                    a2t = ediv(a2, p); a3t = ediv(a3, my); a4t = ediv(ediv(a4, p), mx); a6t = ediv(ediv(a6, mx), my)
                    if (a4t*a4t - 4*a6t*a2t) % p == 0:
                        r = mx * ((a6t * pow(a2t, -1, 2)) % 2) if p == 2 else mx * ((-a4t * pow(2*a2t, -1, p)) % p)
                        a = rst(a, r, 0, 0); a1, a2, a3, a4, a6 = a
                        mx *= p; ix += 1
                    else:
                        cp = 4 if _roots_quad(a2t, a4t, a6t, p) else 2; break
                else:
                    cp = 4 if _roots_quad(1, a3t, -a6t, p) else 2; break
            return dict(model=a, kod='I%d*' % (ix + iy - 5), f=vD - ix - iy + 1, c=cp, red='additive', ap=0, vD=vD)
        # sw == 3: triple root, move it to 0
        if p == 2: r = b % 2
        elif p == 3: r = (-d) % 3
        else: r = (-b * pow(3, -1, p)) % p
        a = rst(a, p*r, 0, 0); a1, a2, a3, a4, a6 = a
        a3t = ediv(a3, p*p); a6t = ediv(a6, p**4)
        if (a3t*a3t + 4*a6t) % p != 0:                   # IV*
            cp = 3 if _roots_quad(1, a3t, -a6t, p) else 1
            return dict(model=a, kod='IV*', f=vD - 6, c=cp, red='additive', ap=0, vD=vD)
        t = -p*p * (a6t % 2) if p == 2 else p*p * ((-a3t * half) % p)
        a = rst(a, 0, 0, t); a1, a2, a3, a4, a6 = a
        if a4 % p**4 != 0:
            return dict(model=a, kod='III*', f=vD - 7, c=2, red='additive', ap=0, vD=vD)
        if a6 % p**6 != 0:
            return dict(model=a, kod='II*', f=vD - 8, c=1, red='additive', ap=0, vD=vD)
        a = unscale(a, p)                                 # not minimal: divide out and start again


# ---------------------------------------------------------------- local root numbers
def rizzo_row_2(cm):
    """Rizzo, Table III, row v2(c4,c6,Delta) = (4, >=7, 6), c4' = 3 mod 4 (Kodaira III, f_2 = 5): w_2 = +1 iff
    c4' - 4*c6/2^7 = 7, 11 mod 16 (the sign is PARI-calibrated: it is the one that agrees with PARI in --rowtest, not
    read reliably from the paper; it matters only when v2(c6) = 7; see the docstring).
    Returns None outside this row."""
    c4, c6, D = cm
    if val(c4, 2) != 4 or c6 == 0 or val(c6, 2) < 7 or val(D, 2) != 6: return None
    c4u = c4 >> 4
    if c4u % 4 != 3: return None
    return 1 if (c4u - 4*ediv(c6, 2**7)) % 16 in (7, 11) else -1

def local_root_number(p, T):
    """(w_p, how) from the output T of tate(); w_p None if not implemented here"""
    if T['red'] == 'good': return 1, 'good'
    if T['red'] == 'split': return -1, 'split multiplicative'
    if T['red'] == 'nonsplit': return 1, 'non-split multiplicative'
    cm = c_inv(T['model']); c4, c6, D = cm; vD = val(D, p)
    if p >= 5:
        if c4 == 0 or 3*val(c4, p) >= vD:                     # v(j) >= 0: potentially good
            e = 12 // math.gcd(12, vD)
            w = {1: 1, 2: legendre(-1, p), 6: legendre(-1, p), 3: legendre(-3, p), 4: legendre(-2, p)}[e]
            return w, 'additive, potentially good, e = %d (Rohrlich 1993; Rizzo 2003, section 1 (5))' % e
        return legendre(-1, p), 'additive, potentially multiplicative: (-1/p) (Rohrlich 1993; Rizzo 2003, section 1 (4))'
    if c4 != 0 and 3*val(c4, p) < vD:                       # v(j) < 0: additive, potentially multiplicative, p = 2 or 3
        if p == 3: return -1, 'additive, potentially multiplicative: (-1/3) (Rohrlich 1993; Rizzo 2003, section 1 (4))'
        c6u = c6 >> val(c6, 2)
        return (1 if c6u % 4 == 3 else -1), "additive, potentially multiplicative: w_2 = -c6' mod 4, c6' = %d mod 4 (Rohrlich 1993; Rizzo 2003, section 1 (4))" % (c6u % 4)
    if p == 2:
        w = rizzo_row_2(cm)
        if w is not None:
            return w, "Rizzo 2003, Table III, row (4, >=7, 6), c4' = 3 mod 4: c4' = %d mod 16, 4 c6/2^7 = %d mod 16" % ((c4 >> 4) % 16, (4*ediv(c6, 2**7)) % 16)
    return None, 'additive at %d, v(c4, c6, Delta) = (%s, %s, %d): not implemented (would trust PARI)' % (
        p, val(c4, p) if c4 else 'inf', val(c6, p) if c6 else 'inf', vD)


def n_components(kod):
    """number of irreducible components of the special fibre of the Neron model (geometric), from the Kodaira symbol"""
    fixed = {'I0': 1, 'II': 1, 'III': 2, 'IV': 3, 'I0*': 5, 'IV*': 7, 'III*': 8, 'II*': 9}
    if kod in fixed: return fixed[kod]
    if kod.endswith('*'): return 5 + int(kod[1:-1])
    return int(kod[1:])

def local_data(a):
    c4, c6, D = c_inv(a)
    fac = factor_certified(D)
    rows = {}
    for p in sorted(fac):
        T = tate(a, p); T['w'], T['w_how'] = local_root_number(p, T); rows[p] = T
    N = 1
    for p, T in rows.items():
        N *= p**T['f']
        # f_p comes from Ogg's formula in tate(), so that formula is not re-asserted; the known range of f_p is:
        assert T['f'] == {'good': 0, 'split': 1, 'nonsplit': 1}.get(T['red'], T['f']), ('f_p', p, T)
        assert T['red'] != 'additive' or 2 <= T['f'] <= {2: 8, 3: 5}.get(p, 2), ('f_p out of range', p, T)
    ws = [T['w'] for T in rows.values()]
    w = None if None in ws else -math.prod(ws)
    return dict(fac=fac, rows=rows, N=N, w_local=w, D=D, c4=c4, c6=c6)


# ---------------------------------------------------------------- a_n by counting points
def primes_upto(M):
    import numpy as np
    s = np.ones(M + 1, bool); s[:2] = False
    for i in range(2, math.isqrt(M) + 1):
        if s[i]: s[i*i::i] = False
    return [int(p) for p in np.nonzero(s)[0]]

def ap_count(a, p):
    """p + 1 - #E(F_p) for a model with good reduction at p"""
    import numpy as np
    if p == 2:
        a1, a2, a3, a4, a6 = a
        n = 1 + sum(1 for x in range(2) for y in range(2) if (y*y + a1*x*y + a3*y - x**3 - a2*x*x - a4*x - a6) % 2 == 0)
        return 3 - n
    b2, b4, b6, b8 = b_inv(a)
    x = np.arange(p, dtype=np.int64)
    g = ((4*x + b2 % p) % p * x + (2*b4) % p) % p
    g = (g * x + b6 % p) % p
    sq = np.zeros(p, np.int64); sq[(x*x) % p] = 1          # sq[g] = 1 iff g is a square mod p (0 included)
    chi = 2*sq[g] - 1; chi[g == 0] = 0
    return -int(chi.sum())

def integer_roots(a):
    """(e1, e2, e3) if the model is y^2 = (x - e1)(x - e2)(x - e3) with integers e_i and e1 = 0, else None"""
    a1, a2, a3, a4, a6 = a
    if a1 or a3 or a6: return None
    d = a2*a2 - 4*a4; r = math.isqrt(d) if d >= 0 else -1
    if r < 0 or r*r != d or (a2 + r) % 2: return None
    e = (0, (-a2 + r)//2, (-a2 - r)//2)
    assert all(x**3 + a2*x*x + a4*x == 0 for x in e)
    return e

def ap_count_roots(e, p):
    """the same count for y^2 = x(x - e2)(x - e3): a_p = -sum_x chi(x) chi(x - e2) chi(x - e3), chi read off one table"""
    import numpy as np
    x = np.arange((p + 1)//2, dtype=np.int64)
    chi = -np.ones(p, np.int8); chi[(x*x) % p] = 1; chi[0] = 0
    s = chi * np.roll(chi, e[1] % p) * np.roll(chi, e[2] % p)      # np.roll(chi, k)[x] = chi[x - k mod p]
    return -int(s.sum(dtype=np.int64))

def ap_table(a, LD, M):
    """a_p for p <= M: counted at good p (on the p-minimal model of Tate's algorithm if p | Delta), +-1 or 0 at bad p.
    For a model y^2 = x(x - e2)(x - e3) the count uses ap_count_roots, checked against ap_count for every p < 3000."""
    bad = LD['rows']; ap = {}; e = integer_roots(a)
    for p in primes_upto(M):
        if p in bad and bad[p]['red'] != 'good': ap[p] = bad[p]['ap']
        else:
            model = bad[p]['model'] if p in bad else a     # p | Delta of the given model but good: count on the minimal model
            if e is not None and p > 2 and p not in bad:
                ap[p] = ap_count_roots(e, p)
                if p < 3000: assert ap[p] == ap_count(model, p)
            else:
                ap[p] = ap_count(model, p)
            assert ap[p]**2 <= 4*p, ('Hasse', p, ap[p])
    return ap

def an_from_ap(ap, LD, M):
    """a_1..a_M (index 0 unused) by multiplicativity and the Hecke recursion"""
    bad = LD['rows']
    spf = list(range(M + 1))
    for p in ap:
        if p*p > M: break
        for k in range(p*p, M + 1, p):
            if spf[k] == k: spf[k] = p
    an = [0]*(M + 1); an[1] = 1
    for n in range(2, M + 1):
        p = spf[n]; m = n; k = 0
        while m % p == 0: m //= p; k += 1
        if not (p in bad and bad[p]['red'] != 'good'):
            x0, x1 = 1, ap[p]
            for _ in range(k - 1): x0, x1 = x1, ap[p]*x1 - p*x0
            apk = x1
        else:
            apk = ap[p]**k
        an[n] = apk * an[m]
    return an


# ---------------------------------------------------------------- outward-rounded fixed point (integers scaled by 2^PREC)
def cdiv(x, d): return -((-x) // d)
def mul_lo(x, y): return (x*y) >> PREC                    # x, y >= 0
def mul_hi(x, y): return -((-(x*y)) >> PREC)

def atan_inv(k, Q):
    """arctan(1/k) * 2^Q as an integer interval"""
    s = 0; j = 0; err = 0; pw = k
    while True:
        term = (1 << Q) // ((2*j + 1) * pw)
        if term == 0: break
        s += term if j % 2 == 0 else -term; err += 1; j += 1; pw *= k*k
    return s - err - 1, s + err + 1                         # rounding of each term < 1, truncation < first omitted term <= 1

PI50 = 314159265358979323846264338327950288419716939937510       # pi * 10^50, truncated (a check only)

def pi_interval():
    Q = PREC + 20
    a5, b5 = atan_inv(5, Q); a2, b2 = atan_inv(239, Q)
    lo = 16*a5 - 4*b2; hi = 16*b5 - 4*a2
    lo, hi = lo >> 20, cdiv(hi, 1 << 20)
    assert hi - lo < 2**20 and lo * 10**50 < (PI50 + 1) << PREC and hi * 10**50 > PI50 << PREC   # agrees with 50 digits
    return lo, hi

def exp_neg(tlo, thi):
    """[lo, hi] containing exp(-t) for every t in [tlo, thi], 0 <= tlo <= thi < ONE"""
    assert 0 <= tlo <= thi < ONE
    def partial(t, lower):
        Tlo = Thi = ONE; s_lo = s_hi = ONE; k = 0
        while True:
            k += 1
            Tlo = (Tlo * t) // (k * ONE); Thi = cdiv(Thi * t, k * ONE)
            if k % 2: s_lo -= Thi; s_hi -= Tlo
            else: s_lo += Tlo; s_hi += Thi
            if Thi == 0 or (k % 2 == (1 if lower else 0) and Thi < 2): break
        return s_lo if lower else s_hi
    lo = partial(thi, True); hi = partial(tlo, False)
    # partial sums ending with a negative term are below exp(-t), ending with a positive term above (t < 1); a few ulp margin
    return lo - 4, hi + 4

def rate(N, ynum, yden, PI):
    """[lo, hi] containing 2 pi y / sqrt(N) for y = ynum/yden"""
    s = math.isqrt(N << (2*PREC))                          # s <= sqrt(N) 2^PREC < s + 1
    lo = (2*PI[0]*ynum << PREC) // (yden*(s + 1)); hi = cdiv(2*PI[1]*ynum << PREC, yden*s)
    return lo, hi

def q_interval(N, ynum, yden, PI):
    tlo, thi = rate(N, ynum, yden, PI)
    return exp_neg(tlo, thi)

def series(an, q, divide, M):
    """enclosure of sum_{n<=M} c_n q^n, c_n = a_n/n (divide) or a_n, plus the tail bound for n > M"""
    qlo, qhi = q; plo = phi = ONE; lo = hi = 0
    for n in range(1, M + 1):
        plo = mul_lo(plo, qlo); phi = mul_hi(phi, qhi)
        c = an[n]
        if c == 0: continue
        d = n if divide else 1
        if c > 0: lo += (c*plo) // d; hi += cdiv(c*phi, d)
        else: lo += (c*phi) // d; hi += cdiv(c*plo, d)
    qM1 = mul_hi(phi, qhi)                                  # upper bound for q^(M+1)
    om = ONE - qhi                                          # lower bound for 1 - q
    assert om > 0
    if divide: tail = cdiv(2*qM1*ONE, om)
    else: tail = cdiv(2*qM1*((M + 1)*ONE - M*qlo) * ONE, om*om)
    return lo - tail, hi + tail

def terms_needed(N, y, divide, target=TARGET):
    x = 2*math.pi*y/math.sqrt(N); M = 10
    while True:
        tail = 2*math.exp(-x*(M + 1))/(1 - math.exp(-x)) if divide else 2*math.exp(-x*(M + 1))*(M + 1)/(1 - math.exp(-x))**2
        if tail < target/4: return M
        M = int(M*1.05) + 10

def fmt(x, dig=10):
    """an integer interval end as a decimal, for printing only"""
    return '%.*f' % (dig, x / ONE)

def lower_str(lo, dig=9):
    """a decimal <= lo / 2^PREC"""
    if lo < 0: return '-' + upper_str(-lo, dig)
    q, r = divmod((lo * 10**dig) // ONE, 10**dig); return '%d.%0*d' % (q, dig, r)
def upper_str(hi, dig=9):
    """a decimal >= hi / 2^PREC"""
    if hi < 0: return '-' + lower_str(-hi, dig)
    q, r = divmod(cdiv(hi * 10**dig, ONE), 10**dig); return '%d.%0*d' % (q, dig, r)
def iv(I, dig=9): return '[%s, %s]' % (lower_str(I[0], dig), upper_str(I[1], dig))


# ---------------------------------------------------------------- L(E,1) and the guard
def evaluate(an, N, w, PI, Mmax, tL=TARGET, tF=0.1):
    """the two enclosures of L(E,1), the Fricke sum F at y = 1 and the residual at y = 6/5, all from (an, N, w)"""
    out = {}
    M1 = terms_needed(N, 1, True, tL); assert M1 <= Mmax
    q1 = q_interval(N, 1, 1, PI)
    S1 = series(an, q1, True, M1)
    out['L_y1'] = (S1[0] + w*S1[0], S1[1] + w*S1[1]) if w == 1 else (S1[0] - S1[1], S1[1] - S1[0])
    yn, yd = Y2
    Mb = terms_needed(N, yd/yn, True, tL); assert Mb <= Mmax
    qa = q_interval(N, yn, yd, PI); qb = q_interval(N, yd, yn, PI)
    Sa = series(an, qa, True, Mb); Sb = series(an, qb, True, Mb)
    out['L_y2'] = (Sa[0] + (Sb[0] if w == 1 else -Sb[1]), Sa[1] + (Sb[1] if w == 1 else -Sb[0]))
    MF = terms_needed(N, 1, False, tF); assert MF <= Mmax
    out['F'] = series(an, q1, False, MF)
    # residual R = f(i/(y sqrt N)) - w y^2 f(i y / sqrt N) at y = 6/5:  f(i/(y sqrt N)) = sum a_n qb^n, f(i y/sqrt N) = sum a_n qa^n
    MR = terms_needed(N, yd/yn, False, tF); assert MR <= Mmax
    Fa = series(an, qa, False, MR); Fb = series(an, qb, False, MR)
    k = yn*yn; kd = yd*yd                                    # y^2 = k/kd
    if w == 1: R = (Fb[0] - cdiv(k*Fa[1], kd), Fb[1] - (k*Fa[0]) // kd)
    else: R = (Fb[0] + (k*Fa[0]) // kd, Fb[1] + cdiv(k*Fa[1], kd))
    out['R'] = R
    out['M'] = (M1, Mb, MF, MR)
    out['overlap'] = max(out['L_y1'][0], out['L_y2'][0]) <= min(out['L_y1'][1], out['L_y2'][1])
    out['R_has_0'] = R[0] <= 0 <= R[1]
    out['guard'] = out['overlap'] and out['R_has_0']
    return out

def max_terms(N, tL=TARGET, tF=0.1):
    return max(terms_needed(N, Y2[1]/Y2[0], True, tL), terms_needed(N, Y2[1]/Y2[0], False, tF), terms_needed(N, 1, False, tF))


# ---------------------------------------------------------------- reports
def local_report(a, LD, lines):
    lines.append('model a = %s' % (list(a),))
    lines.append('Delta = %s = %s' % (LD['D'], ' * '.join('%d^%d' % (p, e) if e > 1 else str(p) for p, e in sorted(LD['fac'].items()))))
    lines.append('Tate\'s algorithm (exact):')
    for p, T in LD['rows'].items():
        cm = c_inv(T['model'])
        lines.append('  p = %d: v(Delta) = %d, v(c4, c6, Delta) of the p-minimal model = (%s, %s, %d), Kodaira %s, f_p = %d, c_p = %d, %s, a_p = %s, w_p = %s (%s)' % (
            p, T['vD'], val(cm[0], p) if cm[0] else 'inf', val(cm[1], p) if cm[1] else 'inf', val(cm[2], p), T['kod'], T['f'], T['c'],
            T['red'], T['ap'], T['w'], T['w_how']))
        if 'trace' in T: lines.append('    (at %d: %s)' % (p, T['trace']))
    lines.append('conductor N = %d = %s' % (LD['N'], ' * '.join('%d^%d' % (p, T['f']) if T['f'] > 1 else str(p) for p, T in LD['rows'].items() if T['f'])))
    lines.append('root number from the local root numbers, w = -prod w_p: %s' % (LD['w_local'] if LD['w_local'] is not None else 'not determined here (an additive case at 2 or 3 not implemented)'))

def pari_crosscheck(a, LD, an, lines, nterms=None):
    # every a_n used (n <= Mmax) is compared with PARI's ellan; an earlier version stopped at n = 20000, while a wrong
    # a_p above about 1.5*10^5 would not move the guard
    try:
        import cypari2
    except ImportError:
        lines.append('PARI cross-check: cypari2 not available, skipped'); return True
    pari = cypari2.Pari(); E = pari.ellinit(list(a))
    gr = pari.ellglobalred(E); Np = int(gr[0]); wp = int(pari.ellrootno(E)); ok = True
    lines.append('PARI cross-check (not used above): conductor %d (%s), root number %d (%s)' % (
        Np, 'agrees' if Np == LD['N'] else 'DISAGREES', wp, 'agrees' if wp == LD['w_local'] else 'DISAGREES'))
    ok &= Np == LD['N'] and wp == LD['w_local']
    for p, T in LD['rows'].items():
        lr = pari.elllocalred(E, p); fp = int(lr[0]); cp = int(lr[3]); wpp = int(pari.ellrootno(E, p))
        same = fp == T['f'] and cp == T['c'] and (T['w'] is None or wpp == T['w'])
        ok &= same
        lines.append('  p = %d: PARI f_p = %d, c_p = %d, w_p = %d: %s' % (p, fp, cp, wpp, 'agrees' if same else 'DISAGREES'))
    m = len(an) - 1 if nterms is None else min(nterms, len(an) - 1); pa = [int(x) for x in pari.ellan(E, m)]
    same = pa == an[1:m + 1]; ok &= same
    lines.append('  a_n for n <= %d: %s with PARI ellan' % (m, 'agree' if same else 'DISAGREE'))
    return ok

def main_curve(write, nopari, ap=None):
    t0 = time.time(); lines = ['lball_exact.py: ' + LABEL, '']
    a = CURVE; LD = local_data(a); local_report(a, LD, lines)
    N = LD['N']; PI = pi_interval()
    Mmax = max_terms(N)
    if ap is None: ap = ap_table(a, LD, Mmax)
    ap = {p: x for p, x in ap.items() if p <= Mmax}; an = an_from_ap(ap, LD, Mmax)
    lines.append('a_p for the %d primes p <= %d: by counting points at good p (Hasse checked for each), from Tate\'s algorithm at bad p' % (len(ap), Mmax))
    ev = evaluate(an, N, 1, PI, Mmax)
    F = ev['F']
    lines.append('')
    lines.append('Fricke sign: F = sum a_n exp(-2 pi n / sqrt N) in %s  (%d terms, tail included)' % (iv(F, 3), ev['M'][2]))
    wF = 1 if F[0] > 0 or F[1] < 0 else None
    lines.append('  F != 0, so w = +1 (f = sum a_n q^n is the newform of level N; f(i/sqrt N) = w f(i/sqrt N), so w = -1 would force F = 0)' if wF else '  F not shown nonzero: w NOT determined this way')
    lines.append('root numbers: local %s, Fricke %s: %s' % (LD['w_local'], wF, 'agree' if wF == LD['w_local'] == 1 else 'DO NOT AGREE'))
    lines.append('')
    lines.append('L(E,1), cut-off y = 1   (%d terms): %s' % (ev['M'][0], iv(ev['L_y1'])))
    lines.append('L(E,1), cut-off y = 6/5 (%d terms): %s' % (ev['M'][1], iv(ev['L_y2'])))
    lines.append('guard: the two enclosures intersect: %s; Fricke residual at y = 6/5 in %s (%d terms) contains 0: %s' % (
        ev['overlap'], iv(ev['R'], 3), ev['M'][3], ev['R_has_0']))
    lo = ev['L_y1'][0]
    positive = lo > 0 and ev['guard'] and wF == 1
    lines.append('certified (the enclosure at y = 1, with N, a_n and w as above): L(E,1) >= %s > 0: %s' % (lower_str(lo), positive))
    assert positive and ev['guard'] and wF == 1 and LD['w_local'] == 1
    assert LD['N'] == 1929934048                     # the value quoted in SECOND_DESCENT.md, now computed here
    lines.append('')
    if not nopari:
        ok = pari_crosscheck(a, LD, an, lines)
        assert ok, 'PARI cross-check disagrees'
    out = '\n'.join(lines) + '\n'
    sys.stdout.write(out)
    sys.stderr.write('[%.0f s]\n' % (time.time() - t0))
    if write: open('lball_exact.txt', 'w').write(out)

NEG_TL = 1e-3      # a forged N or w moves the two values apart by more than 10^-3; a forged a_p at a large p by ~10^-6,
                   # so those are run with the full precision TARGET (and the true N, so no more terms are needed)

def neg_terms(N):
    return max(max(max_terms(n, NEG_TL) for n in (N // 2, 2*N, N // 7, 3*N)), max_terms(N))

def negative(write, ap=None):
    """forged N, forged w, forged a_p: the guard (the two cut-offs) must refuse each"""
    t0 = time.time(); a = CURVE; LD = local_data(a); N = LD['N']; PI = pi_interval()
    forged = [('N/2 (f_2 = 4)', N // 2, 1, None, NEG_TL), ('2N (f_2 = 6)', 2*N, 1, None, NEG_TL),
              ('N/7 (7 left out)', N // 7, 1, None, NEG_TL), ('3N (3 added)', 3*N, 1, None, NEG_TL), ('w = -1', N, -1, None, NEG_TL),
              ('a_2 = +1 (additive read as split)', N, 1, {2: 1}, TARGET), ('a_3 one more', N, 1, 'a3', TARGET),
              ('a_7 = +1 (non-split read as split)', N, 1, {7: 1}, TARGET), ('a_521 = -1 (split read as non-split)', N, 1, {521: -1}, TARGET),
              ('a_719 = -1 (split read as non-split)', N, 1, {719: -1}, TARGET)]
    Mmax = neg_terms(N)
    if ap is None: ap = ap_table(a, LD, Mmax)
    ap0 = {p: x for p, x in ap.items() if p <= Mmax}; an0 = an_from_ap(ap0, LD, Mmax)
    lines = ['lball_exact.py --negative: the guard (L(E,1) at the cut-offs y = 1 and y = 6/5) with forged inputs; each must be refused', '']
    # every forged input goes through evaluate(), the function main_curve() certifies with, and is refused by its guard
    # (the two enclosures meet and the Fricke residual at y = 6/5 contains 0).  Before 26 November 2026 this test used a
    # copy of the overlap test, without the residual.
    def two(an, Nf, wf, tl):
        ev = evaluate(an, Nf, wf, PI, Mmax, tL=tl)
        return ev['L_y1'], ev['L_y2'], ev['guard'], ev
    def why(ev):
        return ', '.join(x for x, y in (('enclosures disjoint', not ev['overlap']), ('residual excludes 0', not ev['R_has_0'])) if y)
    for tl in (NEG_TL, TARGET):
        L1, L2, ok, ev = two(an0, N, 1, tl); dg = 4 if tl == NEG_TL else 10
        lines.append('true inputs, tails below %g: N = %d, w = +1: L(y=1) %s, L(y=6/5) %s, residual %s: guard %s' % (tl, N, iv(L1, dg), iv(L2, dg), iv(ev['R'], 4), 'passes' if ok else 'FAILS'))
        assert ok
    refused = 0
    for name, Nf, wf, ov, tl in forged:
        if ov == 'a3': ov = {3: ap0[3] + 1}
        if ov is None: an = an0
        else: apf = dict(ap0); apf.update(ov); an = an_from_ap(apf, LD, Mmax)
        L1, L2, ok, ev = two(an, Nf, wf, tl); dg = 4 if tl == NEG_TL else 10
        lines.append('forged %-38s (tails below %g) L(y=1) %s, L(y=6/5) %s, residual %s: guard %s' % (name + ':', tl, iv(L1, dg), iv(L2, dg), iv(ev['R'], 4),
                     'passes (NOT REFUSED)' if ok else 'refuses (%s)' % why(ev)))
        refused += not ok
    lines.append('negative tests: %d forged inputs, %d refused' % (len(forged), refused))
    out = '\n'.join(lines) + '\n'; sys.stdout.write(out); sys.stderr.write('[%.0f s]\n' % (time.time() - t0))
    assert refused == len(forged)
    if write: open('lball_exact_negative.txt', 'w').write(out)

def both(write, nopari):
    """main_curve and negative with one table of a_p"""
    a = CURVE; LD = local_data(a); N = LD['N']
    ap = ap_table(a, LD, max(max_terms(N), neg_terms(N)))
    main_curve(write, nopari, ap); negative(write, ap)

def local_only(a):
    LD = local_data(a); lines = []; local_report(a, LD, lines)
    N = LD['N']; M = terms_needed(N, 1, True)
    lines.append('terms the series of L(E,1) would need for a tail below 1e-9: about %.2g' % M)
    try:
        import cypari2
        pari = cypari2.Pari(); E = pari.ellinit(list(a))
        lines.append('PARI cross-check (not used above): conductor %s, root number %s' % (pari.ellglobalred(E)[0], pari.ellrootno(E)))
    except ImportError: pass
    print('\n'.join(lines))

def integral_model(ai):
    """a rational model (strings like '11085/8') made integral by x = x'/u^2, y = y'/u^3: a_i -> a_i u^i, u = lcm of the
    denominators (an isomorphism over Q, so N and w are unchanged)"""
    from fractions import Fraction
    f = [Fraction(x) for x in ai]; u = 1
    for x in f: u = u * x.denominator // math.gcd(u, x.denominator)
    out = [x * u**i for x, i in zip(f, (1, 2, 3, 4, 6))]
    assert all(x.denominator == 1 for x in out)
    return tuple(int(x) for x in out)

def census(path, write=False):
    """N and w by Tate's algorithm and the local root numbers for every closing curve in second_123.jsonl (the curves of
    lproof.py / lball.py), compared with PARI's ellglobalred and ellrootno: a check of the inputs those programs take
    from PARI (their a_n and the sums are not redone here)."""
    import json, cypari2
    pari = cypari2.Pari(); rows = []
    for l in open(path):
        r = json.loads(l)
        if r['status'] != 'closed': continue
        c = [c for c in r['curves'] if c['curve'] == r['by']][0]
        rows.append((tuple(r['pair']), integral_model(c['rank0_curve'])))
    lines = []; agree = undetermined = disagree = 0; small = small_w = 0
    for pair, a in rows:
        E = pari.ellinit(list(a)); Np = int(pari.ellglobalred(E)[0]); wp = int(pari.ellrootno(E))
        try: LD = local_data(a)
        except ValueError as ex:
            lines.append('%d %d: %s' % (pair[0], pair[1], ex)); undetermined += 1; continue
        missing = [p for p, T in LD['rows'].items() if T['w'] is None]
        okN = LD['N'] == Np; okw = LD['w_local'] == wp
        if okN and (okw or LD['w_local'] is None): agree += 1
        else: disagree += 1
        if LD['w_local'] is None: undetermined += 1
        small += Np <= 10**14; small_w += Np <= 10**14 and LD['w_local'] is not None
        lines.append('%d %d: N = %d (%s), w = %s (%s)%s' % (pair[0], pair[1], LD['N'], 'agrees with PARI' if okN else 'DISAGREES with PARI %d' % Np,
            LD['w_local'], 'agrees' if okw else ('PARI %d' % wp if LD['w_local'] is None else 'DISAGREES with PARI %d' % wp),
            '' if not missing else '; w_p not implemented at %s' % missing))
    lines.append('census: %d closing curves; N exact for all that factored, equal to PARI: %s; w determined exactly for %d, '
                 'equal to PARI for all of them: %s; of the %d with N <= 1e14 (the ones lball.py sums), w exact for %d' % (
        len(rows), disagree == 0, len(rows) - undetermined, disagree == 0, small, small_w))
    out = '\n'.join(lines) + '\n'; sys.stdout.write(out)
    assert disagree == 0
    if write: open('lball_exact_census.txt', 'w').write(out)

KOD_PARI = {1: 'I0', 2: 'II', 3: 'III', 4: 'IV', -1: 'I0*', -2: 'II*', -3: 'III*', -4: 'IV*'}
def kod_from_pari(k):
    k = int(k)
    if k in KOD_PARI: return KOD_PARI[k]
    return 'I%d' % (k - 4) if k > 4 else 'I%d*' % (-k - 4)

def tatetest(K, write=False):
    """Tate's algorithm against PARI on K random models (small coefficients, every fourth one replaced by a quadratic twist,
    then scaled by 2^i 3^j to force non-minimal models, and half of them moved by a random (r, s, t))"""
    import random, cypari2
    pari = cypari2.Pari(); rnd = random.Random(20261125); n = bad = 0; cases = {}
    while n < K:
        a = [rnd.randint(-30, 30) for _ in range(5)]
        if n % 4 == 3:                                        # a quadratic twist: additive, potentially multiplicative cases
            c4, c6, D = c_inv(a)
            if D == 0: continue
            d = rnd.choice([-1, 2, -2, 3, -3, 6, -6, 5, -5, 7, -7, 10, 15, -21])
            a = [0, 0, 0, -27*c4*d*d, -54*c6*d**3]
        u = rnd.choice([1, 1, 2, 3, 4, 6, 8, 9])
        a = [x * u**i for x, i in zip(a, (1, 2, 3, 4, 6))]
        if rnd.random() < 0.5: a = list(rst(tuple(a), rnd.randint(-9, 9), rnd.randint(-3, 3), rnd.randint(-9, 9)))
        D = c_inv(a)[2]
        if D == 0: continue
        try: fac = factor_certified(D)
        except ValueError: continue
        E = pari.ellinit(a); Em = pari.ellminimalmodel(E); n += 1
        for p in fac:
            T = tate(a, p); lr = pari.elllocalred(E, p)
            if (T['f'], T['c'], T['kod']) != (int(lr[0]), int(lr[3]), kod_from_pari(lr[1])):
                bad += 1; print('DISAGREE', a, p, T['f'], T['c'], T['kod'], lr)
            if T['red'] in ('split', 'nonsplit') and T['ap'] != int(pari.ellap(Em, p)): bad += 1; print('a_p DISAGREE', a, p)
            if T['red'] == 'good' and ap_count(T['model'], p) != int(pari.ellap(Em, p)): bad += 1; print('a_p (good) DISAGREE', a, p)
            w, how = local_root_number(p, T)
            key = ('p = %d' % p if p < 5 else 'p >= 5') + ', ' + (how.split(':')[0].split(' (')[0].split(', e =')[0] if w is not None else 'w_p not implemented')
            cases[key] = cases.get(key, 0) + 1
            if w is not None and w != int(pari.ellrootno(E, p)): bad += 1; print('w_p DISAGREE', a, p, w, how)
    lines = ['tatetest: %d random models, every prime of Delta: f_p, c_p, Kodaira, a_p and w_p compared with PARI: %d disagreements' % (n, bad)]
    for k in sorted(cases): lines.append('  %6d  %s' % (cases[k], k))
    out = '\n'.join(lines) + '\n'; sys.stdout.write(out)
    if write: open('lball_exact_tatetest.txt', 'w').write(out)

def rowtest(K, write=False):
    """the row of Rizzo's Table III, both sign readings, against PARI's local root number at 2, by v_2(c6) = 7, 8, >= 9"""
    import random, cypari2
    pari = cypari2.Pari(); rnd = random.Random(7); n = 0; tally = {}
    while n < K:
        A = rnd.randrange(1, 10**6, 2); B = rnd.randrange(1, 10**6, 2)
        if rnd.random() < 0.5: A, B = A*A, B*B
        a = (0, B - A, 0, -A*B, 0)
        if rnd.random() < 0.3: a = (rnd.randrange(0, 4)*2, a[1], rnd.randrange(0, 4)*4, a[3], rnd.randrange(0, 8)*8)
        if c_inv(a)[2] == 0: continue
        T = tate(a, 2); cm = c_inv(T['model'])
        if rizzo_row_2(cm) is None: continue
        c4, c6, D = cm; n += 1
        wp = int(pari.ellrootno(pari.ellinit(list(a)), 2))
        c4u = c4 >> 4; t = 4*ediv(c6, 2**7); k = min(val(c6, 2), 9)
        row = tally.setdefault(k, [0, 0, 0])
        row[0] += 1; row[1] += wp == (1 if (c4u + t) % 16 in (7, 11) else -1); row[2] += wp == (1 if (c4u - t) % 16 in (7, 11) else -1)
        assert T['kod'] == 'III' and T['f'] == 5 and int(pari.elllocalred(pari.ellinit(list(a)), 2)[0]) == 5
    lines = []
    for k, (m, ap, am) in sorted(tally.items()):
        lines.append('rowtest: v2(c6) %s %d: %d curves in the row (4, >=7, 6), c4\' = 3 mod 4; PARI w_2 agrees with "c4\' + 4 c6/2^7" on %d, with "c4\' - 4 c6/2^7" on %d' % ('=' if k < 9 else '>=', k, m, ap, am))
    out = '\n'.join(lines) + '\n'; sys.stdout.write(out)
    if write: open('lball_exact_rowtest.txt', 'w').write(out)

if __name__ == '__main__':
    args = sys.argv[1:]
    if '--local' in args: local_only(tuple(int(x) for x in args[args.index('--local') + 1].split(',')))
    elif '--tatetest' in args: tatetest(int(args[args.index('--tatetest') + 1]), '--write' in args)
    elif '--rowtest' in args: rowtest(int(args[args.index('--rowtest') + 1]), '--write' in args)
    elif '--negative' in args: negative('--write' in args)
    elif '--all' in args: both('--write' in args, '--nopari' in args)
    elif '--census' in args: census(args[args.index('--census') + 1], '--write' in args)
    else: main_curve('--write' in args, '--nopari' in args)
