#!/usr/bin/env python3
"""Does the four-class descent need the primes of good reduction?  (An experiment written after a review of 7b4ff7a.)

The one-class criterion of COMPLETE.md, section 7, needs the good primes (the rows R, N, E of Theorem C').  For the
four-class descent of LAST_PAIR.md, section 12, the question decides whether a law at the bad primes alone, like
Theorem C, could explain its exclusions.  This program runs `descend` of tuple_descent.py on ratios of a list with the
bad primes only: 2, the primes of mn(m^2 - n^2) (the supports of the classes), and the primes below 2000 of m^2 + n^2 and
of the octic (larger ones are left out because the disc test cuts the line into p discs, so the result can only be
weaker than with them).  No good prime is used.  Ratios whose tuple space is larger than CAP are skipped.

usage: tuple_badprimes.py LIST [CAP]     LIST: one ratio m,n at the start of each line; CAP default 2^18
"""
import sys, os
from math import prod
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tuple_descent as T
from sympy import factorint

pairs = [tuple(map(int, l.split()[0].split(','))) for l in open(sys.argv[1]) if l.strip()]
cap = int(sys.argv[2]) if len(sys.argv) > 2 else 2 ** 18
tested = emptied = skipped = 0; left = []
for m, n in pairs:
    forms, sup, known = T.magic(m, n)
    if prod(2 ** len(P) for P in sup) > cap: skipped += 1; continue
    S = m * m + n * n; W = m ** 8 + 14 * m ** 4 * n ** 4 + n ** 8
    extra = sorted(q for q in set(factorint(S)) | set(factorint(W)) if 2 < q < 2000)
    rest, used, sizes = T.descend(forms, sup, known, pmax=2, extra_primes=extra)
    tested += 1; emptied += not rest
    if rest: left.append('(%d,%d):%d' % (m, n, len(rest)))
print('ratios tested (tuple space at most %d): %d; skipped: %d' % (cap, tested, skipped))
print('emptied by the bad primes alone: %d; not emptied (tuples left): %d' % (emptied, len(left)))
print('not emptied: ' + ' '.join(left))
