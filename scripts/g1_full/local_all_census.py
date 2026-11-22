#!/usr/bin/env python3
"""Compares local_all.py with the two conditions of Result 1b on every coprime pair 1 <= n < m <= 1000.

For each pair: (i) is it excluded by Result 1b, that is, some prime = 5 mod 8 divides m^2 + n^2 to an odd
power, or m : n lies in one of the 484 classes of modp_classes.txt; (ii) does local_all.report find a prime
of insolubility, among the primes below 700 and the primes dividing 2mn(m^2-n^2)(m^2+n^2)(m^8+14m^4n^4+n^8).
(The bound 700 is passed explicitly here; the default of local_all.py on the command line is 20000.)
Usage: python3 local_all_census.py SHARD NSHARDS   (pairs with m = SHARD mod NSHARDS; writes a JSON file
and prints 'DONE pairs agreements differences undecided soluble'). The run of 22 November 2026 is recorded
in local_all_census.txt. Needs cypari2."""
import math, collections, sys, json, os
import local_all as L
here = os.path.dirname(os.path.abspath(__file__))
cls = collections.defaultdict(set)
for line in open(os.path.join(here, 'modp_classes.txt')):
    p, r = line.split()[:2]; cls[int(p)].add(int(r))
def odd5(N):
    x = N; d = 2
    while d * d <= x:
        if x % d == 0:
            k = 0
            while x % d == 0: x //= d; k += 1
            if d % 8 == 5 and k % 2: return True
        d += 1 if d == 2 else 2
    return x > 1 and x % 8 == 5
def inclass(m, n):
    for p, rs in cls.items():
        if n % p and (m * pow(n, -1, p)) % p in rs: return True
    return False
sh, ns = int(sys.argv[1]), int(sys.argv[2]); tot = agree = 0; diff = []; und_all = []; sol = []
for m in range(2, 1001):
    if m % ns != sh: continue
    for n in range(1, m):
        if math.gcd(m, n) != 1: continue
        tot += 1
        ab = odd5(m * m + n * n) or inclass(m, n)
        k, top, ins, und = L.report(m, n, BOUND=700)
        if und: und_all.append((m, n, und))
        loc = bool(ins)
        if not loc: sol.append((m, n))
        if ab == loc: agree += 1
        else: diff.append((m, n, ab, ins))
json.dump(dict(tot=tot, agree=agree, diff=diff, und=und_all, soluble=sol), open('local_all_census_%d.json' % sh, 'w'))
print('DONE', tot, agree, len(diff), len(und_all), len(sol), flush=True)
