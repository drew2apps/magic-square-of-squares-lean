"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
Claim 2 of H3, with the bundle's own exact test: G1's exact_local.allowed_at on tuple_descent.magic(m, n) at good primes
p = 7 mod 8 where p3mod4.py finds no point: the allowed set must be empty.  Rests on G1's exact_local.py.
usage: p3_methA.py [--write]   (p3_methA.txt; a few s)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit, good, G1
sys.path.insert(0, G1)
import exact_local as X
out = ['z8: good primes p = 3 mod 4 that empty the four-class set, by G1\'s exact test. Written by an AI model; not read by a person.']
for (m, n, p) in [(2, 1, 7), (9, 1, 7), (3, 1, 23), (5, 1, 47), (25, 1, 127), (14, 1, 151), (14 + 151 * 1000003, 1, 151)]:
    assert good(m, n, p)
    forms, sup, known = X.TD.magic(m, n)
    L, al = X.allowed_at(forms, p)
    out.append('(%d,%d) at p = %d (good): %d allowed label tuples' % (m, n, p, len(al)))
emit(out, 'p3_methA.txt')
