#!/usr/bin/env python3
"""Sanity check of the genus: point counts of the five-form curve over F_p against the Weil bound
|N - (p+1)| <= 98 sqrt(p) (genus 49), for all ratios of good reduction at a few primes."""
import math, sys
import numpy as np
worst = 0.0
for p in [int(x) for x in sys.argv[1:]] or [101, 199, 401, 809, 1201]:
    chi = np.full(p, -1, dtype=np.int64); x = np.arange(p, dtype=np.int64); chi[(x * x) % p] = 1; chi[0] = 0
    a = np.concatenate([x, [1]]); b = np.concatenate([np.ones(p, dtype=np.int64), [0]])
    a2, b2, ab = a * a % p, b * b % p, a * b % p
    w = 0.0; cnt = 0
    for m in range(1, p):
        if (m * m - 1) % p == 0 or (m * m + 1) % p == 0: continue
        m2 = m * m % p
        g4 = 2 * (m2 + 1) % p * ((a2 + b2) % p) % p; xx = 8 * m % p * ab % p
        N = int(((1 + chi[(m2 * b2 + a2) % p]) * (1 + chi[(b2 + m2 * a2) % p]) * (1 + chi[g4])
                 * (1 + chi[(g4 - xx) % p]) * (1 + chi[(g4 + xx) % p])).sum())
        w = max(w, abs(N - (p + 1)) / (98 * math.sqrt(p))); cnt += 1
    worst = max(worst, w)
    print(f'p={p}: {cnt} ratios of good reduction, worst |N-(p+1)|/(98 sqrt p) = {w:.3f}')
print('largest ratio:', round(worst, 3), '| first prime with p + 1 - 98 sqrt(p) > 0:',
      next(q for q in range(9000, 10000) if all(q % d for d in range(2, 100)) and q + 1 - 98 * math.sqrt(q) > 0))
