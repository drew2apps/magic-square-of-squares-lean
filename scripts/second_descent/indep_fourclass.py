#!/usr/bin/env python3
"""
indep_fourclass.py -- how the class d of the quartics of the second descent is read off from the four rational
classes of LAST_PAIR.md section 12 (a check asked for by the review of slice s3; exact integers only).

The four forms of LAST_PAIR.md section 12 are alpha0, alpha8, (1+i) gammaL, (1+i) gammaU; at a solution each has
square norm, so its class in K*/K*^2 (K = Q(i)) is that of a rational number: d0, d8, dL, dU.  A form of square
norm has the same class as its conjugate.  From the identities of SecondDescentIdentities.lean:

    v^4 L1 L3 L5 L6 = alpha0 conj(alpha8) N(gamma4),  N(gamma4) in 2S Q*^2   ->  class 2S d0 d8
    v^4 q+ L1 L4   = -gammaL conj(gammaU) alpha0 alpha8                     ->  class 2 d0 d8 dL dU
    v^4 q+ L2 L3   = -gammaL conj(gammaU) conj(alpha0) conj(alpha8)         ->  class 2 d0 d8 dL dU

(gammaL conj(gammaU) = (1+i)gammaL * conj((1+i)gammaU) / ((1+i)(1-i)), and 1/2 has the class of 2, which is
NOT a square in K: 2 = -i (1+i)^2 has the class of i.)  The constants 2S and 2 depend only on the ratio, so a
solution whose four classes are those of a degenerate solution has, for each of the three quartics, the class d
of that degenerate solution.  This script checks the two formulas at both degenerate solutions of the nine open
ratios, with the value of the quartic computed directly from its coefficients (homogeneous form at (u, v), so
t = infinity is (u, v) = (1, 0)), and compares with the d of second_123.jsonl.  Negative check: the formula without
the factor 2 (the false version of an earlier draft of SECOND_DESCENT_CHAIN.md) must fail.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from indep_class import (gmul, gconj, gnorm, is_square_zi, rational_class, form_coeffs, form_value,
                         degenerate_points)

def quartic_coeffs(m, n, name):
    """Coefficients q0..q4 of the quartic in s (exact integers), by polynomial multiplication."""
    S, D = m * m + n * n, m * m - n * n
    lin = {"1": [m, -n], "2": [m, n], "3": [n, m], "4": [n, -m], "5": [1, -1], "6": [1, 1], "+": [S, -2 * D, S]}
    def mul(a, b):
        out = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] += x * y
        return out
    q = [1]
    for ch in name:
        q = mul(q, lin[ch])
    assert len(q) == 5
    return q

def quartic_value(q, u, v):
    """v^4 Q(-i u / v) = sum q_k (-i u)^k v^(4-k), as a Gaussian integer."""
    re, im = 0, 0
    powers = [(1, 0), (0, -1), (-1, 0), (0, 1), (1, 0)]      # (-i)^k
    for k in range(5):
        c = q[k] * u ** k * v ** (4 - k)
        re += c * powers[k][0]
        im += c * powers[k][1]
    return (re, im)

def same_class(x, r):
    """x (Gaussian integer) and the rational integer r have the same class in K*/K*^2."""
    return is_square_zi(gmul(x, (r, 0)))

PAIRS = [(195, 28), (223, 167), (440, 279), (719, 161), (609, 200), (809, 409), (861, 620), (627, 364), (991, 263)]
NAMES = {"1356": "L1L3L5L6", "+14": "q+L1L4", "+23": "q+L2L3"}

def main():
    rec = {}
    for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "second_123.jsonl")):
        r = json.loads(l)
        rec[tuple(r["pair"])] = {c["curve"]: c.get("d") for c in r["curves"]}
    ok = True
    wrong_formula_fails = 0
    for (m, n) in PAIRS:
        S = m * m + n * n
        co = form_coeffs(m, n)
        for (u, v) in degenerate_points(m, n):
            a0, a8, g4, gL, gU = [form_value(c, u, v) for c in co]
            four = []
            for x in (a0, a8, gmul((1, 1), gL), gmul((1, 1), gU)):
                d = rational_class(x)
                assert d is not None, ("a form of square norm without a rational class", m, n, u, v)
                four.append(d)
            d0, d8, dL, dU = four
            line = "%d,%d at (u,v) = (%d,%d): four classes (d0, d8, dL, dU) = %s;" % (m, n, u, v, tuple(four))
            for code, name in NAMES.items():
                val = quartic_value(quartic_coeffs(m, n, code), u, v)
                formula = 2 * S * d0 * d8 if code == "1356" else 2 * d0 * d8 * dL * dU
                good = same_class(val, formula)
                dr = int(rec[(m, n)][code]) if rec[(m, n)].get(code) else None
                good_rec = dr is not None and same_class(val, dr)
                ok &= good and good_rec
                if code != "1356" and not same_class(val, d0 * d8 * dL * dU):
                    wrong_formula_fails += 1
                line += " %s: formula %s, record d = %s %s;" % (name, "ok" if good else "FAILS", dr,
                                                                 "ok" if good_rec else "DIFFERS")
            print(line)
    # step 10 of SECOND_DESCENT_CHAIN.md: partner_forms_half assumes t - 1 != 0, so t = 1 for (719,161) is
    # excluded separately: there 2S N(gamma4) = 2S * 2 = 4S is not a square (S = 2 * 521^2)
    from math import isqrt
    m, n = 719, 161
    S = m * m + n * n
    for (u, v) in [(1, 1), (-1, 1)]:
        g4 = form_value(form_coeffs(m, n)[2], u, v)
        val = 2 * S * gnorm(g4)
        assert isqrt(val) ** 2 != val
    print("719,161 at t = 1 and t = -1: 2S N(gamma4) = %d is not a square, so neither is a solution" % (2 * S * 2))
    print("the false formula d0 d8 dL dU (no factor 2) fails at %d of %d evaluations" % (wrong_formula_fails, 2 * 2 * len(PAIRS)))
    assert wrong_formula_fails == 2 * 2 * len(PAIRS), "negative check: the formula without the factor 2 should fail"
    print("ALL: the class of each quartic is 2S d0 d8 (L1L3L5L6) and 2 d0 d8 dL dU (q+L1L4, q+L2L3), and equals the"
          " recorded d, at both degenerate solutions of all nine ratios: %s" % ok)
    assert ok

if __name__ == "__main__":
    main()
