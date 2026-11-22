#!/usr/bin/env python3
"""Everywhere-local solubility of the G1 system for a pair (m, n):
is there t in P^1(Q_p) with F0(t), F8(t), G4(t), GL(t), GU(t) all squares in Q_p?

p = 2: tree search on classes t = a mod 2^k (and 1/t = a mod 2^k, a even); on a class, a value with
valuation v <= k - 3 has a determined square class.
p odd: recursion on a list of quadratics f_i(u), u in Z_p. Write f_i = p^(c_i) h_i with h_i primitive.
  - u not a root mod p of any h_i: every f_i(u) is p^(c_i) times a unit, a square iff c_i is even and
    h_i(u) is a quadratic residue. Searched exhaustively for p <= 5000, by 3000 random u otherwise
    (if that fails for a large p the answer is 'undecided', which is reported and never counted as insoluble).
  - u = r mod p, r a root of some h_i (at most 10 residues): recurse on f_i(r + p u).
Start: t = u, and 1/t = p u.  A point with some f_i = 0 is found through nearby points when the root is simple.
Primes tested: all primes below 20000 and all primes dividing 2 m n (m^2 - n^2)(m^2 + n^2)(m^8 + 14 m^4 n^4 + n^8);
for the other primes see Result 1 of G1_FULL_SYSTEM.md. Needs cypari2 (factoring, square roots mod p).
Usage: python3 local_all.py 724,665 740,539 ...
"""
import sys, random
import cypari2
pari = cypari2.Pari()

def forms(m, n):
    S = m*m + n*n
    return [(n*n, 0, m*m), (m*m, 0, n*n), (2*S, 0, 2*S), (2*S, -8*m*n, 2*S), (2*S, 8*m*n, 2*S)]

def val(x, p):
    if x == 0: return 10**9
    v = 0
    while x % p == 0: x //= p; v += 1
    return v

def roots_mod(h, p):
    c2, c1, c0 = (x % p for x in h)
    if c2 == 0:
        return [(-c0 * pow(c1, -1, p)) % p] if c1 else []
    d = (c1*c1 - 4*c2*c0) % p
    inv = pow(2*c2, -1, p)
    if d == 0: return [(-c1 * inv) % p]
    if pow(d, (p - 1)//2, p) != 1: return []
    s = int(pari.Mod(d, p).sqrt().lift())
    return [((-c1 + s) * inv) % p, ((-c1 - s) * inv) % p]

def solve(fs, p, depth=0):
    """True, False or None (undecided)"""
    if depth > 80: raise RuntimeError('depth')
    cs = [min(val(x, p) for x in f) for f in fs]
    hs = [tuple(x // p**c for x in f) for f, c in zip(fs, cs)]
    for h, c in zip(hs, cs):          # h constant mod p: its square class is the same on the whole class
        if h[0] % p == 0 and h[1] % p == 0 and (c % 2 == 1 or pow(h[2] % p, (p - 1)//2, p) != 1): return False
    R = set()
    for h in hs: R.update(roots_mod(h, p))
    undecided = False
    if all(c % 2 == 0 for c in cs):
        e = (p - 1)//2
        ok = lambda u: all(pow((h[0]*u*u + h[1]*u + h[2]) % p, e, p) == 1 for h in hs)
        if p <= 5000:
            if any(ok(u) for u in range(p)): return True
        else:
            rng = random.Random(p + depth)
            if any(ok(rng.randrange(p)) for _ in range(3000)): return True
            undecided = True
    for r in sorted(R):
        sub = [(f[0]*p*p, (2*f[0]*r + f[1])*p, f[0]*r*r + f[1]*r + f[2]) for f in fs]
        res = solve(sub, p, depth + 1)
        if res: return True
        if res is None: undecided = True
    return None if undecided else False

def status2(v, k):
    v %= 2**k
    if v == 0: return None
    e = 0
    while v % 2 == 0: v //= 2; e += 1
    if k - e < 3: return None
    return 1 if (e % 2 == 0 and v % 8 == 1) else 0

def tree2(fs, a, k):
    st = [status2(c2*a*a + c1*a + c0, k) for (c2, c1, c0) in fs]
    if 0 in st: return False
    if all(s == 1 for s in st): return True
    if k >= 60: raise RuntimeError('depth')
    return tree2(fs, a, k + 1) or tree2(fs, a + 2**k, k + 1)

def soluble(m, n, p):
    """True, False or None (undecided; also when a search is cut off at its depth limit)."""
    try: return _soluble(m, n, p)
    except RuntimeError: return None

def _soluble(m, n, p):
    fs = forms(m, n); rev = [(c0, c1, c2) for (c2, c1, c0) in fs]
    if p == 2:
        return tree2(fs, 0, 1) or tree2(fs, 1, 1) or tree2(rev, 0, 1)
    a = solve(fs, p)
    if a: return True
    b = solve([(f[0]*p*p, f[1]*p, f[2]) for f in rev], p)
    if b: return True
    return None if (a is None or b is None) else False

def bad_primes(m, n):
    N = 2*m*n*(m*m - n*n)*(m*m + n*n)*(m**8 + 14*m**4*n**4 + n**8)
    return sorted(int(q) for q in pari.factor(N)[0])

def report(m, n, BOUND=20000):
    small = [int(q) for q in pari.primes([2, BOUND])]
    ps = sorted(set(small) | set(bad_primes(m, n)))
    ins = [p for p in ps if soluble(m, n, p) is False]
    und = [p for p in ps if p not in ins and soluble(m, n, p) is None]
    return len(ps), ps[-1], ins, und

if __name__ == '__main__':
    for arg in sys.argv[1:]:
        m, n = (int(x) for x in arg.split(','))
        k, top, ins, und = report(m, n)
        print(arg, 'primes tested', k, 'largest', top, 'insoluble at', ins, 'undecided at', und, flush=True)
