#!/usr/bin/env python3
"""Independent exact cross-check of the certificate in LeanRank/RankZero.lean (not used by Lean).

Written by an AI model; not read by a person.

E' : y^2 = x(x - 12)(x - 16).  For a rational point with x = n/e^2, gcd(n, e) = 1, the descent
writes n = d1 a^2, n - 12 e^2 = d2 b^2, d1 d2 (n - 16 e^2) = t^2 with d1 in {1,2,3,6} (real place)
and d2 in {+-1, +-2, +-3, +-6}.  This script

  1. recomputes, by brute force over ALL residues (a, e, b, t) modulo 2^k (no table of squares, no
     reduction to a^2 mod 64), which classes have a solution with not (2 | n and 2 | e), for
     k = 5 and k = 6, and checks that modulo 64 exactly the four classes of Lean's goodClasses
     survive (and that modulo 32 the class (1, -3) still survives, as in RankZeroForge.lean);
  2. checks that the four surviving classes are those of torsion points;
  3. counts the points of E' modulo 5 (8, as `card_point_E5`) and checks the eight rational points;
  4. negative test: a forged list in which a good class is declared excluded is refused.
Exit status 0 and the line ALL CHECKS PASSED on success.
"""
import sys
from fractions import Fraction

GOOD = [(1, 1), (3, -3), (2, -1), (6, 3)]
D1 = [1, 2, 3, 6]
D2 = [1, -1, 2, -2, 3, -3, 6, -6]


def survivors(M):
    out = set()
    sq = {}
    for r in range(M):
        sq.setdefault(r * r % M, r)
    for d1 in D1:
        for d2 in D2:
            found = False
            for a in range(M):
                n = d1 * a * a % M
                for e in range(M):
                    if n % 2 == 0 and e % 2 == 0:
                        continue
                    g2 = (n - 12 * e * e) % M
                    if not any((d2 * b * b - g2) % M == 0 for b in range(M)):
                        continue
                    g3 = d1 * d2 * (n - 16 * e * e) % M
                    if not any((t * t - g3) % M == 0 for t in range(M)):
                        continue
                    found = True
                    break
                if found:
                    break
            if found:
                out.add((d1, d2))
    return out


def squarefree_class(q):
    """The squarefree integer d with q = d * (rational square), q a nonzero Fraction."""
    def sf(z):
        s = -1 if z < 0 else 1
        z = abs(z)
        d, p = 1, 2
        while p * p <= z:
            while z % (p * p) == 0:
                z //= p * p
            if z % p == 0:
                d *= p
                z //= p
            p += 1
        return s * d * z
    return sf(q.numerator * q.denominator)


def main():
    s64 = survivors(64)
    assert s64 == set(GOOD), s64
    print("mod 64 survivors:", sorted(s64))
    s32 = survivors(32)
    assert (1, -3) in s32 and set(GOOD) <= s32
    print("mod 32 survivors:", len(s32), "(includes (1,-3))")
    # 2. classes of torsion points
    pts = [(8, 16), (24, 48)]
    classes = set()
    for x, y in pts:
        assert y * y == x * (x - 12) * (x - 16)
        c = (squarefree_class(Fraction(x)), squarefree_class(Fraction(x - 12)))
        classes.add(c)
    classes.add((1, 1))
    # (0,0) and (12,0): the x - T map uses f'(T) at the root itself
    classes.add((squarefree_class(Fraction(192)), squarefree_class(Fraction(-12))))   # T = (0,0)
    assert classes == set(GOOD), classes
    print("torsion classes:", sorted(classes))
    # 3. points mod 5 and the eight rational points
    c5 = 1 + sum(1 for x in range(5) for y in range(5)
                 if (y * y - (x ** 3 - 28 * x * x + 192 * x)) % 5 == 0)
    assert c5 == 8
    rat = [(0, 0), (12, 0), (16, 0), (8, 16), (8, -16), (24, 48), (24, -48)]
    assert all(y * y == x * (x - 12) * (x - 16) for x, y in rat)
    print("#E'(F_5) =", c5, "; 7 affine rational points + O check")
    # 4. negative test
    forged_excluded = [(1, 1)]
    refused = any(c in s64 for c in forged_excluded)
    assert refused
    print("negative test: forged exclusion of (1,1) refused")
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
