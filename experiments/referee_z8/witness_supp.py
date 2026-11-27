"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
Explicit p-adic points at the odd SUPPORT primes (p | m n (m - n)(m + n), p odd) for the known and extra tuples of the 7
true phantom pairs (14 ratios).  Candidates t = x0 + p^a s (x0 in {0, 1, -1}, a = 0..8, s random, seeded) and
t = 1/(p^a s), as primitive (u, v); with K = 14, a value z = lam d_j gamma_j(u, v) is a square in Q_p(i) if
  p = 1 mod 4: both components z(iota), z(-iota) mod p^K have even valuation v < K - 1 and a QR unit part;
  p = 3 mod 4: z in Z[i]/p^K has even valuation v < K - 1 and its unit part has QR norm mod p
(lam in {1, eps} at p = 1 mod 4; at p = 3 mod 4 every unit of Z_p is a square in Q_p(i)).  Every lift is a point.
p = 2 is witness2.py's.  Rests on that test (tuple_descent's forms).
usage: witness_supp.py [--write]   (witness_supp.txt; a few s)"""
import sys, os, random
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit, load_phantoms, TRUE_PHANTOMS, pari, fac_proved
import octic_cut as O
K = 14
def vp(x, p):
    v = 0
    while x % p == 0 and v < K: x //= p; v += 1
    return v, x
def gforms(m, n, u, v): return [(m*v, n*u), (n*v, -m*u), ((m+n)*(u-v), (m-n)*(u+v)), ((m-n)*(u-v), (m+n)*(u+v))]
def ok1(z, p, I, M):
    for s in (1, -1):
        c = (z[0] + s * I * z[1]) % M
        if c == 0: return False
        v, w = vp(c, p)
        if v % 2 or v >= K - 1 or pow(w % p, (p - 1) // 2, p) != 1: return False
    return True
def ok3(z, p, M):
    a, b = z[0] % M, z[1] % M
    if a == 0 and b == 0: return False
    v = 0
    while a % p == 0 and b % p == 0 and v < K: a //= p; b //= p; v += 1
    if v % 2 or v >= K - 1: return False
    return pow((a*a + b*b) % p, (p - 1) // 2, p) == 1
P = pari(); R = load_phantoms(); rnd = random.Random(9); tot = 0; found = 0
out = ['z8: explicit points at the odd support primes for the tuples of the 7 true phantom pairs. Written by an AI model; not read by a person.']
for (m, n) in TRUE_PHANTOMS:
    known, extra = R[(m, n)]
    ps = sorted((set(fac_proved(P, m)) | set(fac_proved(P, n)) | set(fac_proved(P, m - n)) | set(fac_proved(P, m + n))) - {2})
    miss = []
    for p in ps:
        M = p ** K
        I = O.sqrtm1(p, K) if p % 4 == 1 else None
        eps = next(g for g in range(2, p) if pow(g, (p - 1) // 2, p) == p - 1)
        lams = [1, eps] if p % 4 == 1 else [1]
        for T in known + extra:
            hit = None; tot += 1
            for trial in range(60000):
                a = trial % 9; x0 = rnd.choice((0, 1, -1, None)); s = rnd.randrange(1, M)
                if x0 is None: u, v = 1, p ** a * s
                else: u, v = x0 + p ** a * s, 1
                if a == 0 and x0 is None: u, v = rnd.randrange(M), 1
                g = gcd(u, v); u //= g; v //= g
                F = gforms(m, n, u, v)
                for lam in lams:
                    zs = [(lam * d * F[j][0], lam * d * F[j][1]) for j, d in enumerate(T)]
                    if all((ok1(z, p, I, M) if p % 4 == 1 else ok3(z, p, M)) for z in zs): hit = True; break
                if hit: break
            if hit: found += 1
            else: miss.append((p, T))
    out.append('(%d,%d): odd support primes %s; tuples without a point found: %s' % (m, n, ps, miss))
out.append('(prime, tuple) cases: %d; with an explicit point: %d' % (tot, found))
emit(out, 'witness_supp.txt')
