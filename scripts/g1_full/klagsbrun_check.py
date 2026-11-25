#!/usr/bin/env python3
"""Is Theorem S1 of SELMER_SQUARE.md an instance of Klagsbrun's lower-bound obstruction?  (slice v5, run 3)

    python3 klagsbrun_check.py [--max M] [--pari] [--negative]

Written by an AI model (Claude); not read by a person.  Experiment; nothing here is counted.

Klagsbrun (Trans. AMS 369 (2017) 3357-3388, arXiv:1201.5408, Lemma 4.1; and Math. Res. Lett., arXiv:1201.5407,
Remark 3.2 and Prop. 4.1): for E/K with E(K)[2] = <T>, M = K(E[2]), M' = K(E'[2]) (E' = E/<T>),
    E has a cyclic 4-isogeny over K                 iff M' = K,
    E has one over M (kernel containing T), not over K   iff M = M'  (M' != K).
The lower bound d_2(E^F/K) >= r_2(K) is proved under M = M'.

Here, for E+- = sigma_exact.model(m, n, +-2), moved to y^2 = x(x^2 + a x + b):

PART A (sympy, general m, n; the model is rebuilt symbolically exactly as sigma_exact.model builds it):
  A1  e = -2 S^3 is a root of the cubic; a = 3e + a2, b = 3e^2 + 2 a2 e + a4 are the closed forms of selmer_square.py;
  A2  b' := a^2 - 4b = -c^2,  c = 8 m n (m - n)(m + n)^3, so M = Q(E[2]) = Q(i);
  A3  b = 4 S^2 H H-,  H H- = S^4 - 4 m^2 n^2 (m^2 - n^2)^2,  H = m^2 (m - n)^2 + n^2 (m + n)^2  (> 0);
  A4  E': X^2 + a' X + b' (a' = -2a) has discriminant 16 b, so M' = Q(E'[2]) = Q(sqrt(H H-));
  A5  halving: on y^2 = x(x^2 + Ax + B), x(2P) = (x^2 - B)^2 / (4 y^2), so 2P = T iff x(P)^2 = B;
      the subgroup <P> (2P = T) is defined over a field L iff x(P) in L iff sqrt(B) in L;
  A6  the other 2-torsion points over Q(i): e1, e2 = (-a +- i c)/2 (roots of x^2 + a x + b), e1 e2 = b; the point
      P with 2P = (e1, 0) has x(P) = e1 + sqrt((e1 - e2)(e1 - 0)), so <P> is defined over Q(i) iff
      w := (e1 - e2) e1 = i c e1 is a square in Q(i); N_{Q(i)/Q}(w) = c^2 b (sympy), so w is a square only if
      H H- is a rational square;
  A7  for q in Q*: q is a square in Q(i) iff q or -q is a square in Q (since Q(i)*^2 meets Q* in Q*^2 u -Q*^2).
  Consequence (written proof, every identity checked here):
      cyclic 4-isogeny over Q(i), kernel containing T     iff b in Q*^2 u -Q*^2  iff H H- in Q*^2  (as H H- > 0);
      cyclic 4-isogeny over Q(i), kernel containing e1/e2  only if H H- in Q*^2;
      cyclic 4-isogeny over Q                             iff M' = Q iff H H- in Q*^2;
      M = M' (Klagsbrun's case)                           iff -H H- in Q*^2: never (H H- > 0), at every ratio.
  So E+- has a cyclic 4-isogeny over Q(i) (any kernel) only if H(m, n) H(m, -n) is a rational square.

PART B (exact, curve by curve, from sigma_exact.model itself, not from the closed forms): every coprime square-type
ratio (m^2 + n^2 = k^2) and twice-square ratio (m^2 + n^2 = 2 k^2) with m <= M (default 1000), both curves: the
cubic's root, a, b, b' recomputed; b' = -c^2; b, -b, H H-, 16 b, and c^2 b are not rational squares (isqrt); so no
cyclic 4-isogeny over Q(i) with any kernel, and none over Q.  (For the square type, H(m, +-n) not a square is L0.7
of SELMER_SQUARE.md, through Lemma D of UNIFORM_RANK.md; for the twice-square type it is Lemma D/B of UNIFORM_RANK.md
used in Theorem 2.  The isqrt checks here re-establish it on these finite lists without either lemma.)

--pari (guide only, not used): PARI's ellisomat over nfinit(t^2 + 1) on the first curves of each type: the 2-power
isogeny class over Q(i) has exactly 4 curves, E adjacent to the other three by degree 2 and with no degree-4 entry
(E' -> E -> E1 is the cyclic 4-isogeny of E', expected, as sqrt(b') = i c is in Q(i); but Q(i) != Q(E'[2])).

--negative: the same test must report a cyclic 4-isogeny on
  N1  y^2 = x(x^2 + 3x + 1)       (b = 1 square: a cyclic 4-isogeny over Q),
  N2  y^2 = x(x^2 + 4x + 2)       (M = M' = Q(sqrt 2): Klagsbrun's case, kernel containing T),
  N3  y^2 = x(x^2 - 288x - 31752) (M = Q(sqrt 2), M' = Q(sqrt -2): kernel containing e1 over M, not T; with --pari
                                    ellisomat over Q(sqrt 2) shows degree-4 entries from E),
  and a forged "square-type" input with H H- replaced by a square must be refused.
"""
import os, sys
from fractions import Fraction as Fr
from math import gcd, isqrt

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import sigma_exact as SE

def say(*a):
    print(*a, flush=True)

def issq(q):
    q = Fr(q)
    if q < 0: return False
    a, b = q.numerator, q.denominator
    return isqrt(a)**2 == a and isqrt(b)**2 == b

def H(m, n): return m**4 - 2*m**3*n + 2*m*m*n*n + 2*m*n**3 + n**4

# ------------------------------------------------------------------ A. symbolic
def part_A():
    import sympy as sp
    m, n = sp.symbols('m n', real=True)
    z, X, Y, x, y, A_, B_ = sp.symbols('z X Y x y A B')
    S = m*m + n*n
    t = sp.symbols('t')
    f = sp.expand((m*m*t*t + n*n) * (n*n*t*t + m*m) * (2*S*t*t - 8*m*n*t + 2*S))
    co = sp.Poly(f, t).all_coeffs()[::-1]           # low degree first, as sigma_exact.model
    assert len(co) == 7 and all(sp.expand(co[k] - co[6 - k]) == 0 for k in range(7))
    a3, a4c, a5, a6c = co[3], co[4], co[5], co[6]
    g = a3 + a4c*z + a5*(z*z - 2) + a6c*(z**3 - 3*z)  # f(t) = t^3 g(t + 1/t)
    assert sp.expand(sp.cancel(f - t**3 * g.subs(z, t + 1/t))) == 0
    out = {}
    for s in (2, -2):
        q = sp.expand((z + s) * g); z0 = -s
        d = [q]
        for _ in range(4): d.append(sp.diff(d[-1], z))
        v = [sp.expand(p.subs(z, z0)) for p in d]
        A = v[1]
        a2, a4, a6 = sp.expand(v[2] / 2), sp.expand(A * v[3] / 6), sp.expand(A*A * v[4] / 24)
        # the map (z, Y) -> (A/(z - z0), A Y/(z - z0)^2): q(z) (A/(z - z0)^2)^2 = cubic(A/(z - z0))
        Xz = A / (z - z0)
        assert sp.simplify(sp.cancel(q * (A/(z - z0)**2)**2 - (Xz**3 + a2*Xz**2 + a4*Xz + a6))) == 0
        nn = n if s == 2 else -n
        Sx = m*m + nn*nn
        e = -2 * Sx**3
        assert sp.expand(e**3 + a2*e*e + a4*e + a6) == 0, 'A1: e is not a root'
        a = sp.expand(3*e + a2); b = sp.expand(3*e*e + 2*a2*e + a4)
        # cubic(X) = (X - e)(X^2 + c1 X + c0), and with x = X - e: x(x^2 + a x + b)
        assert sp.expand((x + e)**3 + a2*(x + e)**2 + a4*(x + e) + a6 - x*(x*x + a*x + b)) == 0
        assert sp.expand(a + 4*(m**6 - m**4*nn**2 - 8*m**3*nn**3 - m*m*nn**4 + nn**6)) == 0, 'A1: a'
        Hp, Hm = H(m, nn), H(m, -nn)
        assert sp.expand(b - 4*Sx**2*Hp*Hm) == 0, 'A1/A3: b'
        c = 8*m*nn*(m - nn)*(m + nn)**3
        bp = sp.expand(a*a - 4*b)
        assert sp.expand(bp + c*c) == 0, "A2: b' = -c^2"
        assert sp.expand(Hp*Hm - (Sx**4 - 4*m*m*nn*nn*(m*m - nn*nn)**2)) == 0, 'A3: H H-'
        assert sp.expand(Hp - (m*m*(m - nn)**2 + nn*nn*(m + nn)**2)) == 0, 'A3: H sum of two squares'
        ap = -2*a
        assert sp.expand(ap*ap - 4*bp - 16*b) == 0, "A4: disc of E' 2-division"
        out[s] = (a, b, c)
    say("A1  model rebuilt symbolically as sigma_exact.model; e = -2S^3 is a root; a, b = closed forms (E+ and E-)")
    say("A2  b' = a^2 - 4b = -(8mn(m - n)(m + n)^3)^2: Q(E[2]) = Q(i)")
    say("A3  b = 4 S^2 H(m,n) H(m,-n);  H H- = S^4 - 4m^2n^2(m^2 - n^2)^2;  H = m^2(m - n)^2 + n^2(m + n)^2 > 0")
    say("A4  E' 2-division X^2 - 2aX + b' has discriminant 16b: Q(E'[2]) = Q(sqrt(H H-))")
    # A5: duplication on y^2 = x(x^2 + A x + B)
    lam = (3*x*x + 2*A_*x + B_) / (2*y)
    x2 = lam**2 - A_ - 2*x
    x2s = sp.cancel(x2.subs(y, sp.sqrt(x*(x*x + A_*x + B_))))
    assert sp.simplify(sp.cancel(x2s - (x*x - B_)**2 / (4*x*(x*x + A_*x + B_)))) == 0
    say("A5  x(2P) = (x^2 - B)^2 / (4 y^2) on y^2 = x(x^2 + Ax + B): 2P = T iff x(P)^2 = B")
    # A6: other 2-torsion over Q(i), and the norm of w = (e1 - e2) e1
    I = sp.I
    for s in (2, -2):
        a, b, c = out[s]
        e1 = (-a + I*c) / 2; e2 = (-a - I*c) / 2
        assert sp.expand(e1*e1 + a*e1 + b) == 0 and sp.expand(e1*e2 - b) == 0
        w = sp.expand((e1 - e2) * e1)
        wc = sp.expand(sp.conjugate(w))
        assert sp.expand(w * wc - c*c*b) == 0, 'A6: norm'
        # translate (e1, 0) to the origin: y^2 = x'(x' + e1 - e2)(x' + e1); its "B" is (e1 - e2) e1
        assert sp.expand((x + e1)*((x + e1)**2 + a*(x + e1) + b) - x*(x + e1 - e2)*(x + e1)) == 0
    say("A6  e1, e2 = (-a +- ic)/2; translating (e1, 0) to 0 gives B = (e1 - e2) e1 = w; N(w) = c^2 b == H H- (mod squares)")
    say("A7  q in Q* is a square in Q(i) iff q in Q*^2 u -Q*^2; a square of Q(i) has a rational square norm")
    say("PART A conclusion: E+- has a cyclic 4-isogeny over Q(i) (any kernel) or over Q only if H(m,n)H(m,-n) is a")
    say("  rational square; Klagsbrun's case M = M' needs -H H- in Q*^2, which never holds (H H- > 0).")
    return True

# ------------------------------------------------------------------ B. exact, curve by curve
def ratios(M, kind):
    out = []
    for m in range(2, M + 1):
        for n in range(1, m):
            if gcd(m, n) != 1: continue
            S = m*m + n*n
            if kind == 'square' and isqrt(S)**2 == S: out.append((m, n))
            if kind == 'twice' and S % 2 == 0 and isqrt(S // 2)**2 == S // 2: out.append((m, n))
    return out

def shifted(m, n, s):
    Md = SE.model(m, n, s)
    a2, a4, a6 = Md['a2'], Md['a4'], Md['a6']
    roots = []
    # the integer roots of the cubic divide a6 (a6 != 0); we find them as the rational root -2S^3 and check that the
    # quadratic cofactor has no rational root (negative discriminant)
    S = m*m + n*n
    e = -2 * S**3
    assert e**3 + a2*e*e + a4*e + a6 == 0
    c1 = a2 + e; c0 = a4 + e*c1
    assert a6 == -e*c0 and c1*c1 - 4*c0 < 0
    a = 3*e + a2; b = 3*e*e + 2*a2*e + a4
    assert (a, b) == (2*e + c1, e*e + c1*e + c0)   # x = X - e
    return a, b

def four_isogeny_test(a, b):
    """For y^2 = x(x^2 + a x + b) over Q with x^2 + ax + b irreducible, M = Q(sqrt(a^2 - 4b)): return the list of
    reasons a cyclic 4-isogeny over M could exist (empty list = none, by A5-A7 generalised to M = Q(sqrt D))."""
    D = a*a - 4*b
    assert D != 0 and not issq(D)
    # squarefree class of D
    def sqclass_eq(u, v):         # u, v nonzero rationals: same class in Q*/Q*^2
        return issq(Fr(u) * Fr(v))
    reasons = []
    # q in Q* is a square in Q(sqrt D) iff q in Q*^2 u D Q*^2
    if issq(b): reasons.append('b square: cyclic 4-isogeny over Q')
    if sqclass_eq(b, D): reasons.append('b in D Q*^2: M = M\' (kernel containing T, over M)')
    # kernel containing e1: w = (e1 - e2) e1 with e1 = (-a + sqrt D)/2: w = (D - a sqrt D)/2, N(w) = (D^2 - a^2 D)/4
    Nw = Fr(D*D - a*a*D, 4)        # = -D b
    assert Nw == -D * b
    if issq(Nw):
        # then decide whether w itself is a square in Q(sqrt D): w = u + v sqrt D, sqrt(w) = p + q sqrt D with
        # p^2 + D q^2 = u, 2 p q = v; p^2 = (u + sqrt(N w))/2 or (u - sqrt(N w))/2 must be a rational square
        u, v = Fr(D, 2), Fr(-a, 2)
        r = Fr(isqrt(Nw.numerator), isqrt(Nw.denominator))
        for p2 in ((u + r) / 2, (u - r) / 2):
            if p2 != 0 and issq(p2):
                p = Fr(isqrt(p2.numerator), isqrt(p2.denominator)); qq = v / (2*p)
                assert p*p + D*qq*qq == u and 2*p*qq == v
                reasons.append('w = (e1 - e2) e1 is a square in M: kernel containing e1, over M')
                break
    return reasons

def part_B(M):
    tot = 0
    for kind in ('square', 'twice'):
        R = ratios(M, kind)
        cnt = 0
        for (m, n) in R:
            for s in (2, -2):
                a, b = shifted(m, n, s)
                nn = n if s == 2 else -n
                S = m*m + n*n
                c = 8*m*nn*(m - nn)*(m + nn)**3
                assert a*a - 4*b == -c*c, "b' != -c^2"
                Hp, Hm = H(m, nn), H(m, -nn)
                assert b == 4*S*S*Hp*Hm and Hp > 0 and Hm > 0
                assert not issq(Hp*Hm), ('H H- square', m, n, s)
                assert not issq(b) and not issq(-b) and not issq(16*b) and not issq(c*c*b)
                assert four_isogeny_test(a, b) == [], (m, n, s)
                cnt += 1
        say(f"B  {kind}-type ratios m <= {M}: {len(R)} ratios, {cnt} curves; on every one: b' = -c^2 (Q(E[2]) = Q(i)),"
            f" b, -b, H H-, c^2 b not squares; four_isogeny_test empty")
        tot += cnt
    return tot

def part_pari(M):
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(2*10**9)
    say("PARI (guide only): ellisomat(E, 2) over nfinit(t^2 + 1)")
    for kind in ('square', 'twice'):
        for (m, n) in ratios(M, kind)[:6]:
            for s in (2, -2):
                a, b = shifted(m, n, s)
                r = pari(f'K=nfinit(t^2+1); E=ellinit([0,{a},0,{b},0],K); [L,Mx]=ellisomat(E,2); [#L, Mx[1,]]')
                say(f"   {kind} ({m},{n}) {'E+' if s == 2 else 'E-'}: class size, degrees from E = {r}")
    for name, (a, b) in (('N2', (4, 2)), ('N3', (-288, -31752))):
        r = pari(f'K=nfinit(t^2-2); E=ellinit([0,{a},0,{b},0],K); [L,Mx]=ellisomat(E,2); [#L, Mx[1,]]')
        r0 = pari(f'E=ellinit([0,{a},0,{b},0]); [L,Mx]=ellisomat(E,2); [#L, Mx[1,]]')
        say(f"   {name} y^2 = x(x^2 + {a}x + {b}) over Q(sqrt 2): {r};  over Q: {r0}")

def negative():
    ok = 0
    r = four_isogeny_test(3, 1); assert r and 'over Q' in r[0]; ok += 1
    say(f"N1  y^2 = x(x^2 + 3x + 1): {r}")
    r = four_isogeny_test(4, 2); assert any("M = M'" in z for z in r); ok += 1
    say(f"N2  y^2 = x(x^2 + 4x + 2): {r}")
    r = four_isogeny_test(-288, -31752); assert r == ['w = (e1 - e2) e1 is a square in M: kernel containing e1, over M']; ok += 1
    say(f"N3  y^2 = x(x^2 - 288x - 31752): {r}  (M = Q(sqrt 2), M' = Q(sqrt -2): Klagsbrun's M = M' fails, yet a")
    say("    cyclic 4-isogeny over M exists with kernel containing e1; his Lemma 4.1(ii) needs the kernel to contain T)")
    # a curve of the family's shape (b' = -c^2, so M = Q(i)) with b a square: y^2 = x(x^2 + 6x + 25), 36 - 100 = -8^2
    r = four_isogeny_test(6, 25); assert r and 'over Q' in r[0]; ok += 1
    say(f"N4  y^2 = x(x^2 + 6x + 25) (b' = -8^2 as in the family, b = 5^2): {r}")
    # a forged family curve: (4,3) E+ with H H- replaced by 7^2 must fail the per-curve assertions of part B
    m, n = 4, 3
    a, b = shifted(m, n, 2)
    S = 25; c = 8*m*n*(m - n)*(m + n)**3
    try:
        bf = 4*S*S*49
        assert a*a - 4*bf == -c*c and not issq(bf)
        refused = False
    except AssertionError:
        refused = True
    assert refused; ok += 1
    say("N5  (4,3) E+ with b forged to 4 S^2 7^2: refused by the per-curve assertions")
    say(f"negative: {ok}/5 refused or flagged as they must be")

def main():
    args = sys.argv[1:]
    M = int(args[args.index('--max') + 1]) if '--max' in args else 1000
    if '--negative' in args:
        negative(); return
    part_A()
    part_B(M)
    negative()
    if '--pari' in args: part_pari(M)
    say("ALL CHECKS PASSED")

if __name__ == '__main__':
    main()
