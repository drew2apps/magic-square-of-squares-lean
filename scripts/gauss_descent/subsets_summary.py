#!/usr/bin/env python3
"""What the first descent says when fewer than five of the forms are asked to be squares.

subsets_small.txt: the output of subsets.py for every coprime ratio m > n >= 1 up to a bound.  A magic square whose four edge
entries are squares has a ratio m:n and a rational t (Theorem 1 of PAPER.md), and its five diagonal entries are squares
exactly when the five forms are.  If the descent set of four of the forms is empty, no such square with that ratio has eight
square entries with those four among them; if the set of three is empty, none has those seven.  A set that is not empty
proves nothing.  The second part looks for rational points t = u/v at which two or more forms are squares: where there is
one, the descent must have left a tuple for that set of forms (it has no right to exclude a real point).  The points with
three or more forms, for a ratio with m^2 + n^2 neither a square nor twice a square, are printed with the edge roots of the
square they give.
usage: subsets_summary.py [file [bound for |u|, v]]"""
import sys, math
names = ['F0', 'F8', 'G4', 'GL', 'GU']
def issq(x): return x >= 0 and math.isqrt(x)**2 == x
def typ(m, n):
    S = m*m + n*n
    return 'square' if issq(S) else 'twice a square' if S % 2 == 0 and issq(S//2) else 'neither'
rows = {}
for l in open(sys.argv[1] if len(sys.argv) > 1 else 'subsets_small.txt'):
    pair, rest = l.split(':', 1); m, n = (int(v) for v in pair.strip('() ').split(','))
    rows[(m, n)] = {tuple(k.split('.')): int(v.split('/')[0]) for k, v in (x.split('=') for x in rest.split())}
M = max(m for m, n in rows); assert all((m, n) in rows for m in range(2, M + 1) for n in range(1, m) if math.gcd(m, n) == 1)
print('%d coprime ratios m > n >= 1 with m <= %d' % (len(rows), M))
for t in ('neither', 'square', 'twice a square'):
    R = [p for p in rows if typ(*p) == t]
    print('m^2 + n^2 %s: %d ratios; every set of k forms empty, k = 5, 4, 3, 2: %s' % (t, len(R), ', '.join(
        str(sum(all(v == 0 for T, v in rows[p].items() if len(T) == k) for p in R)) for k in (5, 4, 3, 2))))
left = {p: sorted('.'.join(T) for T, v in d.items() if len(T) == 3 and v) for p, d in rows.items() if typ(*p) == 'neither'}
left = {p: v for p, v in left.items() if v}
print('ratios of the first kind with a set of three forms that is not empty: %d' % len(left))
B = int(sys.argv[2]) if len(sys.argv) > 2 else 200; bad = 0; npts = 0; nsub = 0
for (m, n) in sorted(rows):
    S = m*m + n*n; found = {}
    for v in range(0, B + 1):
        for u in range(-B, B + 1):
            if math.gcd(u, v) != 1 or (v == 0 and u != 1): continue
            vals = [m*m*v*v + n*n*u*u, n*n*v*v + m*m*u*u, 2*S*(u*u + v*v), 2*S*(u*u + v*v) - 8*m*n*u*v, 2*S*(u*u + v*v) + 8*m*n*u*v]
            sq = tuple(names[j] for j in range(5) if issq(vals[j]))
            if len(sq) >= 2: found.setdefault(sq, (u, v)); npts += 1
    for sq, (u, v) in sorted(found.items()):
        sub = [T for T in rows[(m, n)] if set(T) <= set(sq)]; missing = [T for T in sub if rows[(m, n)][T] == 0]; bad += len(missing); nsub += len(sub)
        if typ(m, n) == 'neither' and len(sq) >= 3:
            e = sorted(abs(x) for x in (m*v + n*u, m*v - n*u, n*v + m*u, n*v - m*u))
            print('  (%d,%d), t = %d/%d: %s are squares; edge roots %s%s; the descent leaves a tuple for every subset of these: %s' % (
                m, n, u, v, ', '.join(sq), e, ' (two equal: a degenerate square)' if len(set(e)) < 4 else ' (DISTINCT)', not missing))
print('rational points with |u|, v <= %d at which two or more forms are squares: %d; subsets of forms that are squares at such a point: %d; of these with an empty descent set: %d' % (B, npts, nsub, bad))
assert bad == 0
