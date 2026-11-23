#!/usr/bin/env python3
"""The ratios left with 1000 < m <= 3000, from the three tallies, and what the route of LAST_PAIR.md does to the list.

Reads the lines "not closed" of wide_tally.txt, wide_2000_2500_tally.txt, wide_2500_3000_tally.txt (the list must be the
committed pairs_wide_open.txt), the ratios of route_wide.txt for which only (1, 1) is left, and adds the partner
(m + n, m - n) (reduced) of each when it is in the list.  usage: wide_open.py"""
import re
from math import gcd, isqrt
left = []
for fn in ('wide_tally.txt', 'wide_2000_2500_tally.txt', 'wide_2500_3000_tally.txt'):
    for l in open(fn):
        if l.startswith('not closed:'): left += [tuple(int(x) for x in a.split(',')) for a in l.split(':', 2)[2].split()]
assert left == [tuple(int(x) for x in a.split(',')) for a in open('pairs_wide_open.txt').read().split()], 'pairs_wide_open.txt is not the list of the tallies'
sq = [p for p in left if isqrt(p[0]**2 + p[1]**2)**2 == p[0]**2 + p[1]**2]
ok = 'generator: True; '; end = 'within the bound: [(1, 1)]; range covers the bound: True; fibre over a = 1 or b = 1: rational t = -1, 1'
done = [tuple(int(x) for x in re.match(r'\((\d+),(\d+)\)', l).groups()) for l in open('route_wide.txt') if ok in l and l.rstrip().endswith(end)]
assert all(p in sq for p in done)
print('ratios left by the two descents with 1000 < m <= 3000: %d; of square type: %d' % (len(left), len(sq)))
print('closed by the route (route_wide.txt): %d: %s' % (len(done), ' '.join('%d,%d' % p for p in done)))
closed = set(done); pairs = []
for m, n in done:
    a, b = m + n, m - n; g = gcd(a, b); a, b = a//g, b//g
    if (a, b) in left: closed.add((a, b)); pairs.append('%d,%d~%d,%d' % (m, n, a, b))
print('their partners in the list: %d: %s' % (len(pairs), ' '.join(pairs)))
print('closed with partners: %d; left: %d' % (len(closed), len(left) - len(closed)))
print('left: ' + ' '.join('%d,%d' % p for p in left if p not in closed))
