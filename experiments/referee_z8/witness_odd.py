"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
Explicit p-adic points for the extra tuples of the 7 true phantom pairs (14 ratios) at every odd prime p of m^2 + n^2
and of the octic outside the supports (all p = 1 mod 4; every d_j is a unit): a residue x in F_p at which no component
vanishes and chi(l_j^+(x)) = chi(l_j^-(x)) = chi(d_j) + c for every j (one c); any lift to Z_p is a point.  Random x
(seeded), at most 200,000 tries.  The only misses are at p = 5 (see witness5.py).  Rests on octic_cut.comps.
usage: witness_odd.py [--write]   (witness_odd.txt; about 16 s)"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit, load_phantoms, TRUE_PHANTOMS, pari, fac_proved
import octic_cut as O
P = pari(); R = load_phantoms(); rnd = random.Random(3); tot = 0; found = 0
out = ['z8: explicit points at the odd primes of m^2 + n^2 and of the octic (outside the supports) for the extra tuples of the 7 true phantom pairs. Written by an AI model; not read by a person.']
for (m, n) in TRUE_PHANTOMS:
    known, extra = R[(m, n)]
    supp = set([2]) | set(fac_proved(P, m)) | set(fac_proved(P, n)) | set(fac_proved(P, m - n)) | set(fac_proved(P, m + n))
    a = (m*m - m*n - n*n)**2 + 3*m*m*n*n; b = (m*m + m*n - n*n)**2 + 3*m*m*n*n
    ps = sorted((set(fac_proved(P, m*m + n*n)) | set(fac_proved(P, a)) | set(fac_proved(P, b))) - supp)
    miss = []
    for p in ps:
        assert p % 4 == 1
        I = O.sqrtm1(p, 1); L = O.comps(m, n, I, p)
        for T in extra:
            e = [O.chi(d, p) for d in T]; tot += 1; ok = None
            for _ in range(200000):
                x = rnd.randrange(p)
                vals = [(aa * x + bb) % p for (_, _, aa, bb) in L]
                if 0 in vals: continue
                ch = [O.chi(v, p) for v in vals]
                if all(ch[2*j] == ch[2*j+1] for j in range(4)):
                    c = (ch[0] + e[0]) % 2
                    if all((ch[2*j] + e[j]) % 2 == c for j in range(4)): ok = x; break
            found += ok is not None
            if ok is None: miss.append((p, T))
    out.append('(%d,%d): primes %s; extra tuples without a point found: %s' % (m, n, ps, miss))
out.append('(prime, extra tuple) cases: %d; with an explicit point: %d' % (tot, found))
emit(out, 'witness_odd.txt')
