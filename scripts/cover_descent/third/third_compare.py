#!/usr/bin/env python3
"""Compare third.py with cover_descent.py (used as a black box), prime by prime, at the 18 primes up to 61.

For a ratio and a prime p:  L = the classes with an explicit local point at p (third.cover_lower: certified),
U = the classes that the exhaustion modulo p^N keeps (third.cover_upper),  R = the classes that the descent allows at p
(cover_descent.Local.allowed, with the depth that cover_descent.run uses).  The true set of classes possible at p lies
between L and U.  The descent passes the WITNESS test at (ratio, p) if L is inside R: it rejects no class for which a local point was
found (at most 400 candidate points are tried for each class, so L may be smaller than the true set).  The STRONGER
inclusion, U inside R, proves that it rejects no class with a local point at that pair.
usage: third_compare.py "m,n ..." | @file [@file ...]"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import cover_descent as cd
import third as bf

pairs = []
for arg in sys.argv[1:]: pairs += open(arg[1:]).read().split() if arg.startswith('@') else arg.split()
n = sound = exact = bracket = weaker = sharper = strong = equalU = other = 0; unsound = []
for pr in pairs:
    m, k = (int(x) for x in pr.split(',')); cl = bf.cands(bf.cover_primes(m, k))
    for p in bf.PRIMES61:
        N = bf.depth(p, m, k); keep, aux = bf.cover_upper(m, k, p, N, cl); low = bf.cover_lower(m, k, p, N, keep, aux, cl)
        U = set(cl[j] for j in np.nonzero(keep)[0]); L = set(cl[j] for j in np.nonzero(low)[0])
        loc = cd.Local(m, k, p, 14 if p == 2 else (8 if p < 50 else 5)); al = loc.allowed()
        R = set(cl) if al == 'all' else set(d for d in cl if loc.dclass(d) in al)
        n += 1; assert L <= U
        if L <= R: sound += 1
        else: unsound.append((m, k, p, sorted(L - R)[:4]))
        exact += (R == L); bracket += (L == U); weaker += (R > U); sharper += (R < U); strong += (U <= R); equalU += (R == U)
        other += not (R == U or R > U or R < U)
print('pairs (ratio, prime): %d (%d ratios, the 18 primes up to 61)' % (n, len(pairs)))
print('witness test (every class with an explicit local point found is allowed by the descent): %d of %d' % (sound, n))
print('the descent allows exactly the certified set: %d; third.py brackets the truth exactly (upper = lower): %d' % (exact, bracket))
print('the descent allows strictly more than the exhaustion keeps: %d; strictly less (and still all certified classes): %d' % (weaker, sharper))
print('the stronger inclusion (every class kept by the exhaustion is allowed by the descent, which proves that no class with a local point is rejected at that pair): %d of %d; descent = exhaustion: %d; neither set contains the other: %d' % (strong, n, equalU, other))
print('UNSOUND at: %s' % unsound if unsound else 'no pair where the descent rejects a class that has an explicit local point')
