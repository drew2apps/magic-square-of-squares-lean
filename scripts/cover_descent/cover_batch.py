#!/usr/bin/env python3
"""cover_descent.py on a list of ratios, one line per ratio.  Good primes up to 150 first; up to 1500 (or the third argument) where that leaves more than the known classes.
usage: cover_batch.py pairs.txt [jobs [largest good prime of the second pass]]"""
import sys
from multiprocessing import Pool
PMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 1500
from cover_descent import run
def one(a):
    m, n = (int(x) for x in a.split(','))
    P, known, cands, info = run(m, n, 150); kd = sorted(set(known.values())); used = 150
    if cands != kd: P, known, cands, info = run(m, n, PMAX); used = PMAX      # more primes only where the first pass leaves more
    return '(%d,%d): %d classes on the primes %s; good primes up to %d; known points: %s; left: %s; %s' % (m, n, 2**len(P), P, used, known or 'none', cands,
        ('no class is possible' if not cands else 'exactly the classes of the known points') if cands == kd else 'MORE THAN THE KNOWN POINTS')
if __name__ == '__main__':
    pairs = open(sys.argv[1]).read().split()
    with Pool(int(sys.argv[2]) if len(sys.argv) > 2 else 2) as pool:
        for line in pool.imap(one, pairs): print(line, flush=True)
