#!/usr/bin/env python3
"""A rational enclosure of the canonical height of the point G of LAST_PAIR.md (the normalisation of PARI and Sage).

On the minimal model of E1, for every rational point Q (Cremona, Prickett and Siksek, J. Number Theory 116 (2006): the
difference h(x(Q)) - hhat(Q) is a sum of local terms; at a prime p it lies between 0 and alpha_p log p with alpha_p from their
Table 1; at the real place Silverman's bounds as quoted in their Theorem 15):
    h - hhat <= BUP := sum_p alpha_p log p + log+|j|/4 + log+|b2/12| + log 2* + 1.946,
    hhat - h <= U   := log|D|/6 + log+|j|/6 + log+|b2/12| + log 2* + 2.14,
where h(x) = log max(|numerator|, denominator) and log 2* is log 2 if b2 != 0.  With Q = 2^k G computed exactly from the
duplication formula for x,
    (h(x(Q)) - BUP) / 4^k <= hhat(G) <= (h(x(Q)) + U) / 4^k.
All arithmetic is exact (fractions; logarithms enclosed in rational intervals), as in scripts/g1_full/height_certify.py, from
which the three helper functions are copied.  Kodaira types and Tamagawa numbers are Sage's (Tate's algorithm).  The statements
of Cremona, Prickett and Siksek were read from the text of their paper with an AI model.
usage: height_enclosure.py m,n d [k]        (as a module: enclose(Em, xG, k) for a minimal model Em)"""
import sys
from fractions import Fraction as Fr

def atanh_bounds(z, terms=60):
    s = Fr(0); p = z; z2 = z*z
    for i in range(terms):
        s += p/(2*i + 1); p *= z2
    return s, s + p/((2*terms + 1)*(1 - z2))
LN2 = tuple(2*x for x in atanh_bounds(Fr(1, 3)))
def ln_bounds(q, bits=96):
    q = Fr(q); assert q > 0
    k = q.numerator.bit_length() - q.denominator.bit_length()
    y = q/Fr(2)**k
    if y < 1: y *= 2; k -= 1
    one = 1 << bits
    ylo = Fr((y.numerator*one)//y.denominator, one); yhi = Fr(-((-y.numerator*one)//y.denominator), one)
    lo = 2*atanh_bounds((ylo - 1)/(ylo + 1))[0]; hi = 2*atanh_bounds((yhi - 1)/(yhi + 1))[1]
    return (lo + k*(LN2[0] if k >= 0 else LN2[1]), hi + k*(LN2[1] if k >= 0 else LN2[0]))
def lnplus_hi(q): return ln_bounds(q)[1] if q > 1 else Fr(0)
def alpha(code, c):
    """Table 1 of Cremona-Prickett-Siksek; code = PARI's code of the Kodaira type, c = Tamagawa number."""
    if c == 1: return Fr(0)
    if code >= 5:
        m = code - 4
        if m % 2 == 0 and c in (2, m): return Fr(m, 4)
        if m % 2 == 1 and c == m: return Fr(m*m - 1, 4*m)
    if code == 3 and c == 2: return Fr(1, 2)
    if code == 4 and c == 3: return Fr(2, 3)
    if code == -1 and c in (2, 4): return Fr(1)
    if code <= -5:
        m = -code - 4
        if c == 2: return Fr(1)
        if c == 4: return Fr(m + 4, 4)
    if code == -4 and c == 3: return Fr(4, 3)
    if code == -3 and c == 2: return Fr(3, 2)
    raise ValueError((code, c))

def enclose(Em, xG, k=7):
    """(lo, hi, BUP, U): rationals with lo <= hhat(G) <= hi, for the point of x-coordinate xG on the minimal model Em."""
    a1, a2, a3, a4, a6 = (int(a) for a in Em.a_invariants())
    b2 = a1*a1 + 4*a2; b4 = 2*a4 + a1*a3; b6 = a3*a3 + 4*a6; b8 = a1*a1*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3*a3 - a4*a4
    D = -b2*b2*b8 - 8*b4**3 - 27*b6*b6 + 9*b2*b4*b6; c4 = b2*b2 - 24*b4; j = Fr(c4**3, D)
    assert int(Em.discriminant()) == D and Em.is_minimal()
    na = Fr(0)
    for ld in Em.local_data():
        p = int(ld.prime().gen()); assert int(ld.discriminant_valuation()) == Em.discriminant().valuation(p)
        na += alpha(int(ld.kodaira_symbol()._pari_code()), int(ld.tamagawa_number()))*ln_bounds(p)[1]
    tail = lnplus_hi(Fr(abs(b2), 12)) + (LN2[1] if b2 else 0)
    BUP = na + lnplus_hi(abs(j))/4 + tail + Fr(1946, 1000)
    U = ln_bounds(abs(D))[1]/6 + lnplus_hi(abs(j))/6 + tail + Fr(214, 100)
    x = Fr(int(xG.numerator()), int(xG.denominator())) if hasattr(xG, 'numerator') and callable(xG.numerator) else Fr(xG)
    for _ in range(k):
        num = x**4 - b4*x*x - 2*b6*x - b8; den = 4*x**3 + b2*x*x + 2*b4*x + b6
        assert den != 0; x = num/den
    hx = ln_bounds(max(abs(x.numerator), x.denominator))
    return (hx[0] - BUP)/4**k, (hx[1] + U)/4**k, BUP, U

if __name__ == '__main__':
    from second_descent import *
    m, n = (int(v) for v in sys.argv[1].split(',')); d = int(sys.argv[2]); k = int(sys.argv[3]) if len(sys.argv) > 3 else 7
    Q1, r1 = family(m, n)['+13']; E1, A1 = weierstrass(d*Q1, r1); EK = E1.change_ring(K); conj = [a for a in K.automorphisms() if a(I) == -I][0]
    vd = d*Q1.change_ring(K)(-I); wd = 1/(-I - r1); Pd = EK(A1*wd, A1*vd.sqrt()*wd*wd); Gs = Pd + EK(conj(Pd[0]), conj(Pd[1])); G = E1(QQ(Gs[0]), QQ(Gs[1]))
    Em = E1.minimal_model(); Gm = E1.isomorphism_to(Em)(G)
    lo, hi, BUP, U = enclose(Em, QQ(Gm[0]), k)
    print('minimal model %s; x(G) = %s' % (list(Em.a_invariants()), Gm[0]))
    print('h - hhat <= %.3f and hhat - h <= %.3f for every rational point (rounded up)' % (float(BUP) + 0.0005, float(U) + 0.0005))
    print('from 2^%d G, in exact arithmetic: %d/10^6 <= hhat(G) <= %d/10^6' % (k, int(lo*10**6), -int(-hi*10**6 // 1)))
    print('floating-point values for comparison: PARI %.6f, Sage %.6f' % (float(pari.ellheight(pari.ellinit([pari(str(a)) for a in Em.a_invariants()]), [pari(str(Gm[0])), pari(str(Gm[1]))])) if 'pari' in globals() else float('nan'), float(Gm.height())))
