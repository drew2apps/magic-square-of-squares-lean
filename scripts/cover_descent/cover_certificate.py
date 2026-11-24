#!/usr/bin/env python3
"""The local eliminations of LAST_PAIR.md, section 11, prime by prime, from two programs (EXPERIMENT).

For each prime p up to a bound: the candidates d (products of the primes of 2 m n (m^2 - n^2)) that cover_descent.py allows
at p (discs, valuations, square classes) and those that cover_bruteforce.py keeps at p (every point modulo p^N, lists of
squares).  Each is a set that contains the class of every point over Q_p.  The line says whether the two sets are the same;
the list of the candidates REMOVED at p is printed (each one is a separate statement: this d is not the class of a point over
Q_p), and then the candidates left by all primes so far.
usage: cover_certificate.py m,n [largest prime, default 61] [largest modulus for the exhaustion, default 5000]"""
import sys
import cover_descent as cd, cover_bruteforce as bf

m, n = (int(x) for x in sys.argv[1].split(',')); pmax = int(sys.argv[2]) if len(sys.argv) > 2 else 61
Mmax = int(sys.argv[3]) if len(sys.argv) > 3 else 5000
P = sorted(bf.factor(2*m*n*(m*m - n*n))); cands = [1]
for q in P: cands += [c*q for c in cands]
cands.sort(); print('(%d,%d): candidates: the %d products of %s' % (m, n, len(cands), P))
left1 = set(cands); left2 = set(cands); same = True
for p in bf.primes_to(pmax):
    L = cd.Local(m, n, p, 14 if p == 2 else (8 if p < 50 else 5)); al = L.allowed()
    s1 = set(cands) if al == 'all' else set(d for d in cands if L.dclass(d) in al)
    N = 1
    while p**(N + 1) <= Mmax: N += 1
    s2 = bf.kept(m, n, cands, p, p**N)[1]
    left1 &= s1; left2 &= s2; same = same and s1 == s2
    rel = 'the same set' if s1 == s2 else ('the descent keeps fewer' if s1 < s2 else ('the exhaustion keeps fewer' if s2 < s1 else 'DIFFERENT SETS, neither inside the other'))
    print('prime %d: allowed by the descent: %d; kept by the exhaustion modulo %d^%d: %d; %s' % (p, len(s1), p, N, len(s2), rel))
    if len(s1 & s2) < len(cands):
        print('  removed at %d by both: %s' % (p, sorted(set(cands) - (s1 | s2))))
        print('  left after the primes up to %d: descent %d, exhaustion %d%s' % (p, len(left1), len(left2), ': ' + str(sorted(left1)) if left1 == left2 and len(left1) <= 8 else ''))
print('left by the descent: %s; left by the exhaustion: %s; the same set at every prime: %s' % (sorted(left1), sorted(left2), same))
