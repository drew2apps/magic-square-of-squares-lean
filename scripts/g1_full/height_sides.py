#!/usr/bin/env python3
"""Result 4: Silverman's height bound on both sides, and the margins of the listings.

dm.py, dm2.py and dm3.py use |h(x(Q)) - hhat(Q)| <= B for every rational point Q, with B = Sage's
silverman_height_bound() (computed by eclib) and hhat = Sage's canonical height (twice Silverman's).
eclib documents B as a bound for h - hhat only.

Sources. J. E. Cremona, M. Prickett, S. Siksek, Height difference bounds for elliptic curves over number
fields, J. Number Theory 116 (2006), in the normalisation used here and for any Weierstrass equation with
integral coefficients (minimal or not):
  - Theorem 1 and Proposition 8: at a finite place the local difference (naive minus canonical) is >= 0,
    so only the real place contributes to an upper bound for hhat - h;
  - Theorem 15, quoted there from J. H. Silverman, The difference between the Weil height and the canonical
    height on elliptic curves, Math. Comp. 55 (1990), p. 737: at the archimedean place
        naive - canonical >= -( log+|D|/6 + log+|j|/6 + log+|b2/12| + log 2* + 2.14 ),
    with D the discriminant, b2 = a1^2 + 4 a2, and 2* = 2 if b2 != 0, 1 if b2 = 0.
So, D being a nonzero integer,
    hhat - h <= U := log|D|/6 + log+|j|/6 + log+|b2/12| + [log 2 if b2 != 0] + 2.14.
eclib's bound for the other side (libsrc/htconst.cc, with Bremner's constant 1.922 instead of 1.946) is
    h - hhat <= B := h(j)/12 + log|D|/6 + log+|j|/6 + log+|b2/12| + [log 2 if b2 != 0] + 1.922,
so U = B - h(j)/12 + 0.218, and U <= B as soon as h(j) >= 2.616.
I read the paper of Cremona, Prickett and Siksek through an automated summary of its PDF and have not seen
Silverman's paper; the two statements above should be compared with the originals.

For the model of every line of dm_evidence.txt this script, in exact rational arithmetic (logarithms are
enclosed in rational intervals by the series for atanh; no floating point enters a comparison):
  1. encloses B and checks that the recorded value (three decimals) is within 0.001 of it; this also shows
     that B was computed for the model whose points the scripts list;
  2. checks U <= B - 2, and that h(j) >= log 32 (the numerator or denominator of j is at least 32);
  3. computes the bound X on |a| from the recorded hG, B (and kappa, or K and rho) and checks that the
     listed bound amax leaves room: X would have to grow by the reported fraction before a multiple could
     be missed. The heights hG are floating-point values computed by Sage, so X is a computed bound, proved
     only to the extent that hG is exact; this step measures the margin against the error of hG, it does
     not bound that error.
No dependencies."""
import os, sys
from fractions import Fraction as Fr

# ---- rational enclosures of logarithms ----
def atanh_bounds(z, terms=60):
    """z rational, 0 <= z <= 1/2. Returns (lo, hi) with lo <= atanh(z) <= hi."""
    s = Fr(0); p = z; z2 = z * z
    for i in range(terms):
        s += p / (2 * i + 1); p *= z2
    return s, s + p / ((2 * terms + 1) * (1 - z2))        # tail <= z^(2N+1) / ((2N+1)(1 - z^2))

LN2 = tuple(2 * x for x in atanh_bounds(Fr(1, 3)))        # ln 2 = 2 atanh(1/3)

def ln_bounds(q, bits=96):
    """q rational > 0. Returns rationals (lo, hi) with lo <= ln q <= hi."""
    q = Fr(q); assert q > 0
    k = q.numerator.bit_length() - q.denominator.bit_length()
    y = q / Fr(2) ** k                                    # 1/2 < y < 2
    if y < 1: y *= 2; k -= 1                              # now 1 <= y < 2
    one = 1 << bits
    ylo = Fr((y.numerator * one) // y.denominator, one)                 # ylo <= y <= yhi
    yhi = Fr(-((-y.numerator * one) // y.denominator), one)
    lo = 2 * atanh_bounds((ylo - 1) / (ylo + 1))[0]
    hi = 2 * atanh_bounds((yhi - 1) / (yhi + 1))[1]
    return (lo + k * (LN2[0] if k >= 0 else LN2[1]), hi + k * (LN2[1] if k >= 0 else LN2[0]))

def lnplus_bounds(q):
    return ln_bounds(q) if q > 1 else (Fr(0), Fr(0))

def add(*iv): return (sum(x[0] for x in iv), sum(x[1] for x in iv))
def scale(iv, c): return (iv[0] * c, iv[1] * c)           # c > 0

assert ln_bounds(2)[0] > Fr(693147, 10**6) and ln_bounds(2)[1] < Fr(693148, 10**6)
assert ln_bounds(Fr(10**40 + 1, 7))[1] - ln_bounds(Fr(10**40 + 1, 7))[0] < Fr(1, 10**20)

# ---- the models ----
def invariants(a1, a2, a3, a4, a6):
    b2 = a1*a1 + 4*a2; b4 = 2*a4 + a1*a3; b6 = a3*a3 + 4*a6
    b8 = a1*a1*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3*a3 - a4*a4
    return b2, b2*b2 - 24*b4, -b2*b2*b8 - 8*b4**3 - 27*b6*b6 + 9*b2*b4*b6          # b2, c4, discriminant

def model(curve, m, n, f):
    S = m*m + n*n
    if curve == 'E04': return (0, 2*S*S, 0, 4*S*S*m*m*n*n, 0)
    if curve == 'D1':
        c3 = 4*S*S*m*m*n*n; c2 = -16*S*m**3*n**3; c1 = 4*S*S*(m*m - n*n)**2; c0 = -16*S*m*n*(m*m - n*n)**2
        return (0, c1, 0, c0*c2, c0*c0*c3)
    return tuple(int(x) for x in f['model'].split(','))

here = os.path.dirname(os.path.abspath(__file__))
count = 0; worst_rec = Fr(0); least_gap = None; least_slack = None; least_j = None
for line in open(os.path.join(here, 'dm_evidence.txt')):
    if line.startswith('#') or not line.strip(): continue
    w = line.split(); curve = w[2]; m, n = (int(x) for x in w[1].split(','))
    f = dict(x.split('=', 1) for x in w[3:] if '=' in x)
    a = model(curve, m, n, f); b2, c4, D = invariants(*a); assert D != 0
    j = Fr(c4**3, D); jsize = max(abs(j.numerator), j.denominator)
    assert jsize >= 32, w[0]                                           # h(j) >= ln 32 > 3.4 > 2.616
    common = add(scale(ln_bounds(abs(D)), Fr(1, 6)), scale(lnplus_bounds(abs(j)), Fr(1, 6)),
                 lnplus_bounds(Fr(abs(b2), 12)), LN2 if b2 else (Fr(0), Fr(0)))
    B = add(common, scale(ln_bounds(jsize), Fr(1, 12)), (Fr(1922, 1000),) * 2)   # enclosure of eclib's bound
    U = add(common, (Fr(214, 100),) * 2)                                        # enclosure of the other side
    Brec = Fr(f['B'])
    dev = max(B[1] - Brec, Brec - B[0]); assert dev <= Fr(1, 1000), (w[0], float(dev))
    gap = Brec - U[1]; assert gap >= 2, (w[0], float(gap))
    # the margin of the listing
    hG = Fr(f['hG']); amax = int(f['amax'])
    if curve == 'E04': X = (2 * Brec + 2 * ln_bounds(2 * (m*m + n*n) * n*n)[1]) / hG
    elif curve == 'D1': X = 2 * Brec / hG
    else: X = (2 * Fr(f['K']) + 4 * (ln_bounds(max(m, n))[1] + LN2[1])) / (2 * hG) + Fr(abs(int(f['rho'])), 2)
    slack = Fr(amax + 1) / X - 1                                       # a multiple is missed only if X_true >= amax + 1
    assert slack >= Fr(2, 100), (w[0], float(slack))
    count += 1; worst_rec = max(worst_rec, dev)
    least_gap = gap if least_gap is None else min(least_gap, gap)
    least_slack = slack if least_slack is None else min(least_slack, slack)
    least_j = jsize if least_j is None else min(least_j, jsize)
    if '-v' in sys.argv: print(w[0], curve, 'B in [%.6f, %.6f]' % (B[0], B[1]), 'U <= %.4f' % U[1], 'X = %.3f' % X, 'amax = %d' % amax)
assert count == 121
print('1. lines: %d; recorded B within %.4f of the enclosure of eclib\'s formula in every line' % (count, worst_rec))
print('2. B - U >= %.3f in every line (U bounds hhat - h); smallest max(|num j|, den j) has %d digits' % (least_gap, len(str(least_j))))
print('3. listed amax against the bound X computed from the recorded heights: X would have to grow by at least %.1f%% in every line' % (100 * least_slack))
print('height_sides OK')
