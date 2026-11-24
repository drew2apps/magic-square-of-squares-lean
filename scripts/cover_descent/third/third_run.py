#!/usr/bin/env python3
"""Run third.py on a list of ratios: the classes possible at every prime up to B (default 61), as an upper bound
(exhaustion modulo p^N; N as in third.depth for p <= 61, N = 1 above) and as a certified lower bound (explicit local
points).  If the two final sets differ, the depth N is raised by one at the primes up to 17 (while p^N stays below 5000) and the ratio is run again (at most twice).
usage: third_run.py "m,n m,n ..." [B]  |  third_run.py @file [B]     (a file of pairs m,n separated by blanks)"""
import sys
import numpy as np
import third as bf

B = int(sys.argv[2]) if len(sys.argv) > 2 else 61
PRIMES = [q for q in range(2, B + 1) if all(q % r for r in range(2, int(q**0.5) + 1))]

def run(m, n, extra):
    cl = bf.cands(bf.cover_primes(m, n)); up = np.ones(len(cl), bool); lo = np.ones(len(cl), bool); exact = 0
    for p in PRIMES:
        N = bf.depth(p, m, n) if p <= 61 else 1
        for _ in range(extra):
            if p <= 17 and p**(N + 1) <= 5000: N += 1
        keep, aux = bf.cover_upper(m, n, p, N, cl); low = bf.cover_lower(m, n, p, N, keep, aux, cl)
        up &= keep; lo &= low; exact += bool((keep == low).all())
    return sorted(cl[j] for j in np.nonzero(up)[0]), sorted(cl[j] for j in np.nonzero(lo)[0]), len(cl), exact

arg = sys.argv[1]
pairs = open(arg[1:]).read().split() if arg.startswith('@') else arg.split()
for pr in pairs:
    m, n = (int(x) for x in pr.split(','))
    for extra in (0, 1, 2):
        up, lo, ncl, exact = run(m, n, extra)
        if up == lo: break
    print('(%d,%d): %d classes; possible at every prime up to %d: upper bound %s, certified lower bound %s; %s; the local set is exact at %d of %d primes; known classes %s%s'
          % (m, n, ncl, B, up, lo, 'EQUAL' if up == lo else 'DIFFERENT', exact, len(PRIMES), sorted(bf.known_cover_classes(m, n)), '' if extra == 0 else '; depth raised %d time(s) at the primes up to 17' % extra), flush=True)
