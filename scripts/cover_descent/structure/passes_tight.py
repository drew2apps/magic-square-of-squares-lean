#!/usr/bin/env python3
"""Is the repaired predicate `Passes` (ClassTable.lean) exactly the table of Theorem C?  A test by exhaustion (numpy).

Lean proves one direction (`passes_of_row`): every class of the row passes modulo every power of p.  This program checks
the other direction at a finite power, which no Lean statement covers yet: a class outside the row does not pass.  For
an odd prime p, a power M = p^K and a ratio (m, n) with p | m, it goes through every point (u, v) modulo M, not both
divisible by p, and keeps the classes d in {1, e, p, p e} (e a non-residue) for which the definition holds there:

    ClassAt M (p^3) d m n u v:  each of the four norms is w^2 modulo M for some w not divisible by p^3, and
                                d times each of the eight products is (a + b i)^2 modulo M for some a, b.

`Passes p k d m n` with k + 1 = K is the existence of such a point.  It can only get weaker as K decreases (reduce a
point modulo a smaller power), so a class outside the row that passes at K would also pass at every smaller power.
Expected, from the row of the table: p = 3 (mod 8): the units only; p = 7 (mod 8): all four; p = 5: e only if 5 || m,
1 and e if 25 | m.  Written after a review of 7b4ff7a.

usage: passes_tight.py            the cases below; about three minutes
       passes_tight.py --weak     the same with the first definition of Passes (no roots of the norms), as a control
                                  that the test can fail: there the class p passes at p = 3 (mod 8)
"""
import sys
import numpy as np

def tables(p, K, weak=False):
    M = p ** K; a = np.arange(M, dtype=np.int64)
    root = np.zeros(M, dtype=bool); root[(a * a % M)[a % p ** 3 != 0]] = True       # w^2 with p^3 not dividing w
    if weak: root[:] = True                                                         # the first definition: no roots
    G = np.zeros((M, M), dtype=bool)                                                # (a + b i)^2 = (a^2 - b^2, 2ab)
    for b in range(M):
        G[(a * a - b * b) % M, (2 * a * b) % M] = True
    return M, root, G

def gm(x1, y1, x2, y2, M): return (x1 * x2 - y1 * y2) % M, (x1 * y2 + y1 * x2) % M

def passing(p, K, m, n, reps, M, root, G):
    out = set()
    U = np.arange(M, dtype=np.int64)
    for v in range(M):
        u = U if v % p else U[U % p != 0]
        P = [((m * v) % M, (n * u) % M), ((n * v) % M, (m * u) % M),
             (((m + n) * (u - v)) % M, ((m - n) * (u + v)) % M), (((m - n) * (u - v)) % M, ((m + n) * (u + v)) % M)]
        ok = np.ones(len(u), dtype=bool)
        for x, y in P: ok &= root[(x * x + y * y) % M]
        if not ok.any(): continue
        P = [(np.broadcast_to(x, u.shape)[ok], np.broadcast_to(y, u.shape)[ok]) for x, y in P]
        prods = []
        for c8 in (1, -1):
            for cL in (1, -1):
                for cU in (1, -1):
                    x, y = gm(*P[0], P[1][0], c8 * P[1][1], M)
                    x, y = gm(x, y, P[2][0], cL * P[2][1], M)
                    prods.append(gm(x, y, P[3][0], cU * P[3][1], M))
        for d in reps:
            if d in out: continue
            good = np.ones(len(prods[0][0]), dtype=bool)
            for x, y in prods: good &= G[(d * x) % M, (d * y) % M]
            if good.any(): out.add(d)
        if len(out) == len(reps): break
    return out

def row(p, m, d):
    """the condition Row p m d of ClassTable.lean"""
    if d % p:
        return p != 5 or d % 5 in (2, 3) or m % 25 == 0
    return (d // p) % p != 0 and p % 8 in (1, 7)

CASES = [(3, 7, [(3, 1), (3, 2), (6, 1), (9, 1), (9, 2), (18, 5), (27, 1)]),
         (5, 5, [(5, 1), (5, 2), (10, 1), (15, 2), (25, 1), (25, 2), (50, 3), (125, 1)]),
         (7, 4, [(7, 1), (7, 3), (14, 1), (49, 1), (49, 2)]),
         (11, 4, [(11, 1), (11, 2), (121, 1)])]
if __name__ == '__main__':
    weak = '--weak' in sys.argv; bad_out = bad_in = total = 0
    if weak: CASES = [(p, K, r) for p, K, r in CASES if p % 8 == 3]
    for p, K, ratios in CASES:
        e = next(r for r in range(2, p) if pow(r, (p - 1) // 2, p) == p - 1); reps = [1, e, p, p * e]
        M, root, G = tables(p, K, weak)
        for m, n in ratios:
            got = passing(p, K, m, n, reps, M, root, G); want = {d for d in reps if row(p, m, d)}
            total += 1; bad_out += len(got - want); bad_in += len(want - got)
            print('p = %d, modulo %d^%d, ratio %d : %d: classes that pass %s; the row of the table %s%s'
                  % (p, p, K, m, n, sorted(got), sorted(want), '' if got == want else '  DIFFERENT'))
    print('%scases %d; classes outside the row that pass: %d; classes of the row that do not pass: %d' % ('first definition (no roots), ' if weak else '', total, bad_out, bad_in))
