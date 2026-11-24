#!/usr/bin/env python3
"""The law at 5 (LAST_PAIR.md, section 17), by a complete check modulo 25.  Python only; about a minute.

Statement: m, n coprime, 5 divides exactly one of m, n, m - n, m + n and 25 does not divide it; u, v not both divisible by 5;
d an integer such that d times each of the eight products P0 * P8^± * PL^± * PU^± is a Gaussian square.  Then d is not
= ±1 modulo 5.  (At a solution 5 does not divide d, by Theorem E, so d = ±2 modulo 5: a non-residue.)

The hypothesis and the property of the eight products can be read modulo 25, so the statement follows from a check of all
residues (m, n, u, v, d) modulo 25, which is what this program does; no norm has to be a square.  It also prints how often
a non-residue passes, and repeats the check in the form that a Lean proof would use (the points of the projective line
modulo 25 and the four norms squares modulo 25: `survivors` and `killed` of lean_cert.py).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import lean_cert as L
M = 25; sq, g = L.masks(M)
qr = [d for d in range(1, M) if d % 5 in (1, 4)]; nqr = [d for d in range(1, M) if d % 5 in (2, 3)]
pairs = tested = passing = pts_n = 0; kills_ok = True
for m in range(M):
    for n in range(M):
        if m % 5 == 0 and n % 5 == 0: continue
        if sum(x % 5 == 0 and x % 25 != 0 for x in (m, n, (m - n) % 25, (m + n) % 25)) != 1 or any(x % 25 == 0 for x in (m, n, (m - n) % 25, (m + n) % 25)): continue
        pairs += 1
        for u in range(M):
            for v in range(M):
                if u % 5 == 0 and v % 5 == 0: continue
                Qs = [L.Q(L.forms(m, n, u, v), f) for f in L.FLAGS]
                for d in qr:
                    tested += 1; passing += all(((d*X) % M, (d*Y) % M) in g for X, Y in Qs)
                pts_n += any(all(((d*X) % M, (d*Y) % M) in g for X, Y in Qs) for d in nqr)
        kills_ok &= all(L.killed(d, L.survivors(m, n, M, 5), M) for d in qr)
print('pairs (m, n) modulo 25 with 5 exactly dividing one of m, n, m-n, m+n: %d' % pairs)
print('(point, d = ±1 mod 5) tested: %d; passing: %d' % (tested, passing))
print('points at which some d = ±2 mod 5 passes: %d' % pts_n)
print('the same in the form of the Lean checker (projective points, norms squares): no d = ±1 mod 5 survives for any pair: %s' % kills_ok)
