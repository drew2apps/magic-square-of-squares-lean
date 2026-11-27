"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
Claim 2 of H3 says good primes p = 3 mod 4 never remove anything.  At such p, Q_p(i) is unramified, every rational
unit d_j is a square in it, so the four-class set is everything or nothing; gamma_j(t) is a square in F_{p^2} iff its
norm F_j(t) is a square in F_p, and F_j has no zero on P^1(F_p).  So the set is everything iff some t in P^1(F_p) makes
F0, F8, GL, GU all nonzero squares, and empty otherwise.  This lists, for every prime p = 3 mod 4 up to 1160, the good
residues r = m/n mod p with no such t.  Rests on: that argument; enumeration in exact integers.
usage: p3mod4.py [--write]   (p3mod4.txt)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit, isprime, good
def Fs(m, n, u, v):
    return (m*m*v*v + n*n*u*u, n*n*v*v + m*m*u*u, 2*((m*u - n*v)**2 + (m*v - n*u)**2), 2*((m*u + n*v)**2 + (m*v + n*u)**2))
out = ['z8: good primes p = 3 mod 4, p <= 1160: good residues r = m/n mod p at which no t in P^1(F_p) makes all four forms nonzero squares (the four-class set is then empty). Written by an AI model; not read by a person.']
emp = []
for p in range(3, 1161):
    if p % 4 != 3 or not isprime(p): continue
    sq = set(x * x % p for x in range(1, p))
    pts = [(x, 1) for x in range(p)] + [(1, 0)]
    e = []; ng = 0
    for r in range(1, p):
        if not good(r, 1, p): continue
        ng += 1
        if not any(all(f % p in sq for f in Fs(r, 1, u, v)) for u, v in pts): e.append(r)
    if e: emp.append(p)
    out.append('p = %d (p mod 8 = %d): %d good residues, %d with no point %s' % (p, p % 8, ng, len(e), e if e else ''))
out.append('primes p = 3 mod 4 at which some good residue has an empty set: %s' % emp)
emit(out, 'p3mod4.txt')
