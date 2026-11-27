"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
G1's four-class descent (exact_local.descend_exact) on the 24 phantom ratios of H1 with EVERY bad prime (supports,
m^2 + n^2, the octic; PARI's factor with isprime) and the good primes up to 20,000; the tuples left against H1's lists.
Rests on G1's exact_local.py, the repository's tuple_descent.py (at p = 2 its one-sided depth-12 test).
usage: full_descent.py [--write]   (full_descent.txt; about 7 minutes)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit, G1, load_phantoms, pari, fac_proved
sys.path.insert(0, G1)
import exact_local as X
TD = X.TD
P = pari()
def factorint_proved(x):
    x = abs(int(x))
    if x <= 1: return {}
    f = P.factor(x); o = {int(q): int(e) for q, e in zip(f[0], f[1])}
    for q in o: assert P.isprime(q)
    return o
TD.factorint = factorint_proved
out = ['z8: the four-class descent with every bad prime and good primes <= 20000 on H1\'s 24 ratios. Written by an AI model; not read by a person.']
keep = []
for (m, n), (known, extra) in sorted(load_phantoms().items()):
    forms, sup, kn = TD.magic(m, n)
    assert sorted(kn) == sorted(known), ('known tuples differ', m, n)
    bad = set(q for S in sup for q in S) | set(fac_proved(P, m*m + n*n)) | set(fac_proved(P, m**8 + 14*m**4*n**4 + n**8))
    primes = sorted(bad | set(int(q) for q in P.primes([2, 20000])))
    S, used = X.descend_exact(forms, sup, kn, primes)
    assert set(kn) <= set(S)
    le = sorted(set(S) - set(kn))
    if le: keep.append((m, n))
    out.append('(%d,%d): %d tuples left; extra tuples left: %s' % (m, n, len(S), ('the %d of H1' % len(le)) if le == sorted(extra) else (le if le else 'none')))
out.append('ratios that keep extra tuples (phantoms): %d %s' % (len(keep), keep))
emit(out, 'full_descent.txt')
