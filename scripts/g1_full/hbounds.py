"""Two-sided bounds between the naive height h(x(P)) and the canonical height (Sage's normalisation, twice
Silverman's) for the rational points of an integral Weierstrass model, as in height_certify.py:
    h - hhat <= BUP = sum_p (alpha_p + ord_p(D/Dmin)/6) log p + log+|j|/4 + log+|b2/12| + log 2* + 1.946
    hhat - h <= U   = log|D|/6 + log+|j|/6 + log+|b2/12| + log 2* + 2.14
(Cremona, Prickett, Siksek 2006: equation (1), Proposition 8 and Table 1, Theorem 15 quoting Silverman 1990;
log 2* is present when b2 != 0).  Exact rational arithmetic; the results are returned as floats rounded up."""
from fractions import Fraction as Fr
from sage.all__sagemath_schemes import ZZ

def atanh_bounds(z, terms=60):
    s = Fr(0); p = z; z2 = z * z
    for i in range(terms):
        s += p / (2 * i + 1); p *= z2
    return s, s + p / ((2 * terms + 1) * (1 - z2))
LN2 = tuple(2 * x for x in atanh_bounds(Fr(1, 3)))
def ln_bounds(q, bits=96):
    q = Fr(q); assert q > 0
    k = q.numerator.bit_length() - q.denominator.bit_length()
    y = q / Fr(2) ** k
    if y < 1: y *= 2; k -= 1
    one = 1 << bits
    ylo = Fr((y.numerator * one) // y.denominator, one); yhi = Fr(-((-y.numerator * one) // y.denominator), one)
    lo = 2 * atanh_bounds((ylo - 1) / (ylo + 1))[0]; hi = 2 * atanh_bounds((yhi - 1) / (yhi + 1))[1]
    return (lo + k * (LN2[0] if k >= 0 else LN2[1]), hi + k * (LN2[1] if k >= 0 else LN2[0]))
def lnplus_hi(q): return ln_bounds(q)[1] if q > 1 else Fr(0)

def alpha(code, c):
    """Table 1 of Cremona-Prickett-Siksek; code = PARI's code of the Kodaira type, c = Tamagawa number."""
    if c == 1: return Fr(0)
    if code >= 5:                                   # I_m, m = code - 4
        m = code - 4
        if m % 2 == 0 and c in (2, m): return Fr(m, 4)
        if m % 2 == 1 and c == m: return Fr(m * m - 1, 4 * m)
    if code == 3 and c == 2: return Fr(1, 2)        # III
    if code == 4 and c == 3: return Fr(2, 3)        # IV
    if code == -1 and c in (2, 4): return Fr(1)     # I0*
    if code <= -5:                                  # I_m*, m = -code - 4
        m = -code - 4
        if c == 2: return Fr(1)
        if c == 4: return Fr(m + 4, 4)
    if code == -4 and c == 3: return Fr(4, 3)       # IV*
    if code == -3 and c == 2: return Fr(3, 2)       # III*
    raise ValueError((code, c))

def height_bounds_exact(ainv, E):
    a1, a2, a3, a4, a6 = [int(v) for v in ainv]
    b2 = a1*a1 + 4*a2; b4 = 2*a4 + a1*a3; b6 = a3*a3 + 4*a6; b8 = a1*a1*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3*a3 - a4*a4
    D = -b2*b2*b8 - 8*b4**3 - 27*b6*b6 + 9*b2*b4*b6; c4 = b2*b2 - 24*b4; j = Fr(c4**3, D)
    assert int(E.discriminant()) == D and [int(v) for v in E.ainvs()] == [a1, a2, a3, a4, a6]
    na = Fr(0); prod = 1
    for p, e in ZZ(D).factor():
        p = int(p); e = int(e); ld = E.local_data(p); emin = int(ld.discriminant_valuation())
        assert (e - emin) % 12 == 0 and e >= emin
        na += (alpha(int(ld.kodaira_symbol()._pari_code()), int(ld.tamagawa_number())) + Fr(e - emin, 6)) * ln_bounds(p)[1]
        prod *= p ** e
    assert prod == abs(D)
    tail = lnplus_hi(Fr(abs(b2), 12)) + (LN2[1] if b2 else 0)
    BUP = na + lnplus_hi(abs(j)) / 4 + tail + Fr(1946, 1000)
    U = ln_bounds(abs(D))[1] / 6 + lnplus_hi(abs(j)) / 6 + tail + Fr(214, 100)
    return BUP, U
def height_bounds(ainv, E):
    BUP, U = height_bounds_exact(ainv, E)
    return float(BUP) + 1e-9, float(U) + 1e-9
