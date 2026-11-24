#!/usr/bin/env python3
"""The two point counts used in Proposition T of PERFECT_TRIANGLE.md, by plain enumeration (Python only).

E: W^2 = (X - 1)(X - 2)(X + 2) has good reduction at an odd prime modulo which the roots 1, 2, -2 are distinct, as at 5
and 13 (checked below).  The number of points over F_p, with the point at infinity:
"""
for p in (5, 13):
    assert all(((a - b) % p) for a, b in ((1, 2), (1, -2), (2, -2))), 'bad reduction'
    n = 1 + sum(1 for X in range(p) for W in range(p) if (W*W - (X - 1)*(X - 2)*(X + 2)) % p == 0)
    print('#E(F_%d) = %d' % (p, n))
print('5 and 13 split in Q(i), so the residue fields of Q(i) above them are F_5 and F_13; the torsion of E over Q(i) injects into both groups, so its order divides gcd(8, 16) = 8, and E has 8 torsion points over Q already (pt_tors.txt).')
