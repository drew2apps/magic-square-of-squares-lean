#!/usr/bin/env python3
"""The first descent restricted to a subset of the five forms.

For a subset T of {F0, F8, G4, GL, GU} only the tuples (delta_j), j in T, are kept, with the same local sets.  If no tuple
is left, the forms of T cannot all be squares at the same rational t: a statement about the hyperelliptic curve
Z^2 = prod_{j in T} f_j(t), on which a solution gives a rational point with every f_j a square.
For a subset of three forms this curve has genus 2, and Magma's TwoCoverDescent (the two-cover descent of Bruin and Stoll)
computes its fake 2-Selmer set.  The tuples left here correspond to the elements of that set at which each of the three norms
is a square.  magma_twocover_raw.txt has Magma's answers; twocover_check.py compares the two counts.

The two-cover descent works modulo Q*: two tuples that differ by the same rational factor in every form (for instance by 2,
whose class is that of i) are one element there.  So the number to compare is the number of tuples up to such a factor.

usage: subsets.py m,n [m,n ...]      prints, for every subset of 2 to 5 forms, the number of tuples left (primes up to 200 and the bad
                                     primes up to 20000), followed by /k when they form only k classes up to a common rational factor"""
import sys, itertools, numpy as np
import gauss_descent as gd, scan
from real_points import is_gsquare, is_sq
def sub(m, n, pmax=200):
    S = m*m + n*n; A = m**4 - 2*m**3*n + 2*m*m*n*n + 2*m*n**3 + n**4; B = m**4 + 2*m**3*n + 2*m*m*n*n - 2*m*n**3 + n**4
    bad = set(gd.factor(2*m*n)) | set(gd.factor(m*m - n*n)) | set(gd.factor(S)) | scan.small_factors(A) | scan.small_factors(B)
    D = gd.dsets(m, n); subsets = [T for k in (2, 3, 4, 5) for T in itertools.combinations(range(5), k)]
    masks = {T: np.ones(tuple(len(D[j]) for j in T), dtype=bool) for T in subsets}; dead = {}
    for p in sorted(set([q for q in bad if q <= 20000] + [q for q in gd.PR if q <= pmax])):
        k = 5 if p > 50 else max(5, 8 if p == 2 else 6)
        L, W = gd.local_set(m, n, p, k); cd = [[L.cls(x) for x in d] for d in D]
        idxs = []
        for w in W:
            idx = [[i for i, c in enumerate(cd[j]) if gd.match(c, w[j])] for j in range(5)]
            idxs.append(idx)
        for T in subsets:
            if T in dead: continue
            al = np.zeros(masks[T].shape, dtype=bool)
            for idx in idxs:
                if all(idx[j] for j in T): al[np.ix_(*[idx[j] for j in T])] = True
            masks[T] &= al
            if not masks[T].any(): dead[T] = p
    orbits = {}
    for T in subsets:                                    # tuples up to a common rational factor (the two-cover descent works modulo Q*)
        if T in dead or masks[T].sum() > 64: continue
        tup = [tuple(D[j][i] for j, i in zip(T, idx)) for idx in zip(*np.nonzero(masks[T]))]
        def same(x, y):
            r = [gd.gmul(a, b) for a, b in zip(x, y)]                  # the quotients, modulo squares
            return is_sq(gd.gnorm(r[0])) and all(is_gsquare(gd.gmul(r[0], rk)) for rk in r[1:])
        reps = []
        for t in tup:
            if not any(same(t, q) for q in reps): reps.append(t)
        orbits[T] = len(reps)
    return dead, {T: int(masks[T].sum()) for T in subsets if T not in dead}, orbits
names = ['F0', 'F8', 'G4', 'GL', 'GU']
if __name__ == '__main__':
    for a in sys.argv[1:]:
        m, n = (int(v) for v in a.split(','))
        dead, alive, orbits = sub(m, n)
        cnt = {T: 0 for T in dead}; cnt.update(alive)
        print('(%d,%d):' % (m, n), ' '.join('%s=%d%s' % ('.'.join(names[j] for j in T), c, '/%d' % orbits[T] if T in orbits and orbits[T] != c else '') for T, c in sorted(cnt.items(), key=lambda kv: (len(kv[0]), kv[0]))), flush=True)
