"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
An exact argument for a cut at a bad prime p = 1 mod 4 that divides the octic m^8 + 14 m^4 n^4 + n^8 and not
2mn(m^2 - n^2) (so every d_j is a unit at p).  There two of the eight linear components l (of the Gaussian forms, in
Z_p x Z_p) have the same root x0 mod p.  If they belong to different forms j_a != j_b, no third component vanishes at
x0, and v_p(resultant) = 1 (root difference of valuation 1; iota mod p^3), then on the disc of x0 one of l_a, l_b
always has odd valuation, so the disc has no allowed point; everywhere else both are units with l_b = kappa l_a mod p.
A form's two components have the same class, so every allowed point has e_ja + e_jb = chi(kappa).  A tuple whose labels
(chi(d_j)) break one such equation is not allowed at p.  Rests on: that argument; exact integers.
usage: octic_cut.py [--write]   (the "open" pair (11875080,6497849) at p = 5462461 -> octic_cut.txt; under 1 s)
       (as a module: analyse(m, n, p, tuples))"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit
from sympy import isprime          # (a proof of primality is not needed here: p comes from PARI's factor with isprime)

def sqrtm1(p, k):
    g = 2
    while pow(g, (p - 1) // 2, p) != p - 1: g += 1
    z = pow(g, (p - 1) // 4, p); M = p ** k
    for _ in range(2 * k + 2): z = (z - (z * z + 1) * pow(2 * z, -1, M)) % M
    assert (z * z + 1) % M == 0
    return z

def chi(a, p): a %= p; assert a; return 0 if pow(a, (p - 1) // 2, p) == 1 else 1

def comps(m, n, I, M):
    """the eight components (j, s, A, B): A u + B v = gamma_j(s iota) (tuple_descent.magic's forms, (1+i) multiplied out)"""
    F = [((0, n), (m, 0)), ((0, -m), (n, 0)), ((m + n, m - n), (-(m + n), m - n)), ((m - n, m + n), (-(m - n), m + n))]
    out = []
    for j, (A, B) in enumerate(F):
        for s in (1, -1): out.append((j, s, (A[0] + s * I * A[1]) % M, (B[0] + s * I * B[1]) % M))
    return out

def analyse(m, n, p, tuples):
    """per root collision: (text, valid, [(tuple, 'ok' or 'VIOLATES')])"""
    k = 3; M = p ** k; I = sqrtm1(p, k); L = comps(m, n, I, M)
    def root(a, b):
        a %= p; b %= p
        return (1, 0) if a == 0 else ((-b * pow(a, -1, p)) % p, 1)
    R = [root(a, b) for (_, _, a, b) in L]
    coll = {}
    for i in range(8):
        for j in range(i + 1, 8):
            if R[i] == R[j]: coll.setdefault(R[i], []).append((i, j))
    res = []
    for x0, pairs in sorted(coll.items()):
        if len(set(q for pr in pairs for q in pr)) != 2: res.append(('more than two components vanish at %s' % (x0,), False, [])); continue
        (ia, ib), = pairs
        ja, sa, aa, ba = L[ia]; jb, sb, ab, bb = L[ib]
        if ja == jb: res.append(('the two components of one form collide', False, [])); continue
        resl = (aa * bb - ab * ba) % M; v = 0; t = resl
        while t % p == 0 and v < k: t //= p; v += 1
        kap = ab * pow(aa, -1, p) % p if aa % p else bb * pow(ba, -1, p) % p
        assert (ab - kap * aa) % p == 0 and (bb - kap * ba) % p == 0
        c = chi(kap, p)
        verdict = [(T, 'ok' if (chi(T[ja], p) + chi(T[jb], p)) % 2 == c else 'VIOLATES') for T in tuples]
        res.append(('forms %s and %s meet at %s, v_p(resultant) = %d: e_%s + e_%s = %d' % ('08LU'[ja], '08LU'[jb], x0, v, '08LU'[ja], '08LU'[jb], c), v == 1, verdict))
    return res

def check(m, n, p, known, extra):
    assert isprime(p) and p % 4 == 1
    for q in (2, m, n, m - n, m + n): assert q % p != 0, 'p divides 2mn(m^2 - n^2)'
    r = analyse(m, n, p, list(known) + list(extra))
    valid = all(v for _, v, _ in r)
    kok = all(x == 'ok' for _, v, vd in r if v for T, x in vd if T in known)
    cut = sorted(set(T for _, v, vd in r if v for T, x in vd if x == 'VIOLATES' and T in extra))
    return r, valid, kok, cut

if __name__ == '__main__':
    m, n, p = 11875080, 6497849, 5462461
    known = [(2, 1, 5377231, 18372929), (2, 1, 36745858, 10754462)]; extra = [(2, 1, 1, 98795483379599), (2, 1, 197590966759198, 2)]
    r, valid, kok, cut = check(m, n, p, known, extra)
    out = ['z8: the octic prime %d for (%d,%d) (octic mod p = %d, m^2 + n^2 mod p = %d). Written by an AI model; not read by a person.' % (p, m, n, (m**8 + 14*m**4*n**4 + n**8) % p, (m*m + n*n) % p)]
    for txt, v, vd in r: out.append('  %s (%s); %s' % (txt, 'argument applies' if v else 'argument does not apply', ', '.join('%s %s' % (T, x) for T, x in vd)))
    out.append('known tuples satisfy every equation: %s; extra tuples cut: %d of %d' % ('yes' if kok else 'NO', len(cut), len(extra)))
    emit(out, 'octic_cut.txt')
