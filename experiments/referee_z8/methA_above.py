"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
Claim 1 of H3, with the bundle's own exact test: at p = 1181, 1201, 1249 (good primes p = 1 mod 4 above 1153), the
first four good residues r at which my_local.py finds a loss, recomputed by H3's method_A (G1's exact_local on
tuple_descent.magic(r, 1)).  Rests on G1's exact_local.py and the repository's tuple_descent.py.
usage: methA_above.py [--write]   (methA_above.txt; a few s)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit, good, H3, G1
os.environ.setdefault('G1_DIR', G1)
sys.path.insert(0, H3)
_a = sys.argv; sys.argv = [sys.argv[0]]
import h3_tables as H
sys.argv = _a
import my_local as M
out = ['z8: losses above 1153 recomputed with the bundle\'s method A. Written by an AI model; not read by a person.']
for p in (1181, 1201, 1249):
    ch = M.chitab(p); k = 0
    for r in range(1, p):
        if not good(r, 1, p): continue
        b = M.allowed(r, 1, p, ch)
        if b == (1 << 16) - 1: continue
        A = H.method_A(r, 1, p)
        out.append('p = %d, r = %d: my test %d vectors, method A %d vectors, %s; missing %s' % (p, r, bin(b).count('1'), bin(A).count('1'), 'equal' if A == b else 'DIFFER', [v for v in range(16) if not b >> v & 1]))
        k += 1
        if k == 4: break
emit(out, 'methA_above.txt')
