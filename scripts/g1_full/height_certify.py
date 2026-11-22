#!/usr/bin/env python3
"""Result 4: certified enclosures of the generator heights, and the listings rechecked with them.

For every line of dm_evidence.txt, in exact rational arithmetic (logarithms enclosed in rational intervals):

  1. Bounds for the difference between the naive height h(x(Q)) and the canonical height hhat(Q) (Sage's
     normalisation), for all rational points Q of the integral model used, from Cremona, Prickett and Siksek,
     J. Number Theory 116 (2006): h - hhat is the sum of local terms Psi_v (their equation (1)); at a finite
     place p, 0 <= Psi_p <= (alpha_p + ord_p(D/Dmin)/6) log p with alpha_p from their Table 1 (Proposition 8);
     at the real place, Silverman's bound as quoted in their Theorem 15. So
        h - hhat <= BUP := sum_p (alpha_p + ord_p(D/Dmin)/6) log p + log+|j|/4 + log+|b2/12| + log 2* + 1.946,
        hhat - h <= U   := log|D|/6 + log+|j|/6 + log+|b2/12| + log 2* + 2.14.
     Kodaira types, Tamagawa numbers and minimal discriminants come from Sage's local_data (Tate's algorithm).
     Checked: BUP <= B + 0.025 and U <= B, where B is the recorded Silverman bound used by dm.py, dm2.py and
     dm3.py on both sides. So |h - hhat| <= BS := B + 0.025.
  2. An enclosure of hhat(G): with Q = 2^k G computed exactly from the duplication formula for x,
        (h(x(Q)) - BUP) / 4^k <= hhat(G) <= (h(x(Q)) + U) / 4^k.
     Checked: the recorded hG lies in the enclosure.
  3. The bound X on |a| of the Demjanenko-Manin argument, recomputed with BS and the lower end of the
     enclosure. Checked: X < amax + 1, so the listed multiples |a| <= amax cover every case.
     For E0L the constants C1 (recorded to two decimals, increased by 0.01 here) and rho are taken from the
     evidence; the derivation of that inequality is not rechecked here.

What this does not certify: the rank, that G generates, and the Demjanenko-Manin derivations themselves.
The statements of Cremona, Prickett and Siksek were read from the text of their paper with an AI model.
Needs passagemath (sage.all__sagemath_schemes) for local_data only; under a minute."""
import os, sys
from fractions import Fraction as Fr
from sage.all__sagemath_schemes import EllipticCurve, QQ, ZZ

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib.util
here = os.path.dirname(os.path.abspath(__file__))

# ---- rational enclosures of logarithms (as in height_sides.py) ----
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

def alpha(code, c, m):
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

def model(curve, m, n, f):
    S = m*m + n*n
    if curve == 'E04': return [0, 2*S*S, 0, 4*S*S*m*m*n*n, 0]
    if curve == 'D1':
        c3 = 4*S*S*m*m*n*n; c2 = -16*S*m**3*n**3; c1 = 4*S*S*(m*m - n*n)**2; c0 = -16*S*m*n*(m*m - n*n)**2
        return [0, c1, 0, c0*c2, c0*c0*c3]
    return [int(x) for x in f['model'].split(',')]

count = 0; widest = Fr(0); least = None; maxdB = None
for line in open(os.path.join(here, 'dm_evidence.txt')):
    if line.startswith('#') or not line.strip(): continue
    w = line.split(); curve = w[2]; m, n = (int(x) for x in w[1].split(','))
    f = dict(x.split('=', 1) for x in w[3:] if '=' in x)
    a1, a2, a3, a4, a6 = model(curve, m, n, f)
    b2 = a1*a1 + 4*a2; b4 = 2*a4 + a1*a3; b6 = a3*a3 + 4*a6; b8 = a1*a1*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3*a3 - a4*a4
    D = -b2*b2*b8 - 8*b4**3 - 27*b6*b6 + 9*b2*b4*b6; c4 = b2*b2 - 24*b4; j = Fr(c4**3, D)
    E = EllipticCurve(QQ, [a1, a2, a3, a4, a6]); assert int(E.discriminant()) == D
    # 1. the two bounds
    na = Fr(0); prod = 1
    for p, e in ZZ(D).factor():
        p = int(p); e = int(e); ld = E.local_data(p); emin = int(ld.discriminant_valuation())
        assert (e - emin) % 12 == 0 and e >= emin
        na += (alpha(int(ld.kodaira_symbol()._pari_code()), int(ld.tamagawa_number()), None) + Fr(e - emin, 6)) * ln_bounds(p)[1]
        prod *= p ** e
    assert prod == abs(D)
    tail = lnplus_hi(Fr(abs(b2), 12)) + (LN2[1] if b2 else 0)
    BUP = na + lnplus_hi(abs(j)) / 4 + tail + Fr(1946, 1000)
    U = ln_bounds(abs(D))[1] / 6 + lnplus_hi(abs(j)) / 6 + tail + Fr(214, 100)
    B = Fr(f['B']); BS = B + Fr(25, 1000)
    assert BUP <= BS and U <= B, (w[0], float(BUP), float(U), float(B))
    # 2. enclosure of hhat(G) from Q = 2^k G
    x = Fr(f['G']); hG = Fr(f['hG']); k = 0
    while 4 ** (k + 1) * hG <= 120000:
        num = x**4 - b4*x*x - 2*b6*x - b8; den = 4*x**3 + b2*x*x + 2*b4*x + b6
        assert den != 0; x = num / den; k += 1
    hx = ln_bounds(max(abs(x.numerator), x.denominator)) if x != 0 else (Fr(0), Fr(0))
    lo = (hx[0] - BUP) / 4**k; hi = (hx[1] + U) / 4**k
    assert lo <= hG <= hi and lo > 0, (w[0], float(lo), float(hG), float(hi))
    # 3. the listing
    amax = int(f['amax'])
    if curve == 'E04': X = (2 * BS + 2 * ln_bounds(2 * (m*m + n*n) * n*n)[1]) / lo
    elif curve == 'D1': X = 2 * BS / lo
    else:
        K = Fr(f['C1']) + Fr(1, 100) + LN2[1] + 2 * BS
        X = (2 * K + 4 * (ln_bounds(max(m, n))[1] + LN2[1])) / (2 * lo) + Fr(abs(int(f['rho'])), 2)
    assert X < amax + 1, (w[0], float(X), amax)
    slack = Fr(amax + 1) / X - 1
    count += 1; widest = max(widest, (hi - lo) / lo); least = slack if least is None else min(least, slack)
    maxdB = BUP - B if maxdB is None else max(maxdB, BUP - B)
    if '-v' in sys.argv: print(w[0], curve, 'k=%d' % k, 'hhat(G) in [%.6f, %.6f]' % (lo, hi), 'recorded %.4f' % hG, 'X < %.3f' % X, 'amax %d' % amax, flush=True)
assert count == 122
print('1. lines: %d; h - hhat <= B + 0.025 and hhat - h <= B in every line (largest BUP - B: %.4f)' % (count, maxdB))
print('2. hhat(G) enclosed for every line; widest enclosure, relative to its lower end: %.3f%%; the recorded hG lies in every one' % (100 * widest))
print('3. with these enclosures the bound on |a| is below amax + 1 in every line; least room: %.1f%%' % (100 * least))
print('height_certify OK')
