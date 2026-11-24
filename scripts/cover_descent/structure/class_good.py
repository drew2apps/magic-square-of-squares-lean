#!/usr/bin/env python3
"""The class at the primes of good reduction: the table of COMPLETE.md, section 7.  Python and numpy.

For a prime p not dividing 2mn(m^2-n^2)(m^2+n^2)(m^8+14m^4n^4+n^8) the curve of the four forms has good reduction, every
point over F_p lifts to a p-adic point, and the class of a point is a unit: for p = 1 (mod 4) it is a residue or a
non-residue (chi(D), D = prod (w_j + x_j); at a point where one norm is 0 that factor is x_j, and where y_j = 0 it is
2x_j), for p = 3 (mod 4) there is one class.  The set of classes that occur depends only on m/n mod p.  This program
lists, for every prime up to the bound, the ratios m/n mod p for which one class is missing or there is no point.

Above 4356 nothing is missing, by the bound of Weil for a curve of genus 33 (for p = 3 mod 4: genus 17, above 1156).

usage: class_good.py                 the table (p = 1 mod 4 up to 4400, p = 3 mod 4 up to 1200) and two checks; 4 minutes
       class_good.py --check B       recompute the rows below B and compare with class_good_table.txt; seconds
"""
import sys, os, re, random
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
def primes(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]
def tables(p):
    rt = np.full(p, -1, dtype=np.int64); a = np.arange(p, dtype=np.int64)
    rt[(a * a) % p] = a
    chi = np.where(rt >= 0, 1, -1); chi[0] = 0
    return rt, chi
def image(p, m, n, rt, chi, U, V):
    """(a residue class occurs, a non-residue class occurs) among the points (u : v) of the arrays U, V"""
    ok = np.ones(len(U), dtype=bool); D = np.ones(len(U), dtype=np.int64)
    for x, y in [(m * V, n * U), (n * V, m * U), ((m + n) * (U - V), (m - n) * (U + V)), ((m - n) * (U - V), (m + n) * (U + V))]:
        x = x % p; y = y % p; N = (x * x + y * y) % p
        w = rt[N]; ok &= (w >= 0)
        D = (D * (np.where(N == 0, x, np.where(y == 0, 2 * x, w + x)) % p)) % p
    c = chi[D][ok]
    return (bool((c == 1).any()), bool((c == -1).any()))
def good(p, m, n):
    W = (pow(m, 8, p) + 14 * pow(m, 4, p) * pow(n, 4, p) + pow(n, 8, p)) % p
    return (m * n * (m * m - n * n) * (m * m + n * n) * W) % p != 0
def table(p):
    """for every m/n mod p of good reduction (n = 1, 2 <= m <= p - 2): which classes occur"""
    rt, chi = tables(p)
    U = np.concatenate([np.arange(p, dtype=np.int64), [1]]); V = np.concatenate([np.ones(p, dtype=np.int64), [0]])
    out = {}
    for m in range(2, p - 1):
        if not good(p, m, 1): continue
        im = image(p, m, 1, rt, chi, U[:96], V[:96]) if p > 200 else (False, False)     # a shortcut only: both seen already
        if not (im[0] and im[1]): im = image(p, m, 1, rt, chi, U, V)
        out[m] = im
    return out
def row(p):
    t = table(p)
    if p % 4 == 3: e = sorted(m for m in t if not (t[m][0] or t[m][1])); a = b = []
    else: e = sorted(m for m in t if t[m] == (False, False)); a = sorted(m for m in t if t[m] == (True, False)); b = sorted(m for m in t if t[m] == (False, True))
    return (e, a, b) if (e or a or b) else None
def fmt(p, r): return '%d E %s R %s N %s\n' % (p, ','.join(map(str, r[0])), ','.join(map(str, r[1])), ','.join(map(str, r[2])))
def load():
    T = {}
    for line in open(os.path.join(HERE, 'class_good_table.txt')):
        if line.startswith('#'): continue
        a = line.rstrip('\n').split(' '); f = lambda i: [int(x) for x in a[i].split(',') if x]
        T[int(a[0])] = (f(2), f(4), f(6))
    return T
HEAD = '# p, then the ratios m/n mod p (n = 1, 2 <= m <= p - 2, good reduction) for which: E there is no point; R every point has a residue as its class; N every point has a non-residue\n'
if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == '--check':
        B = int(sys.argv[2]); T = load(); bad = 0; n = 0
        for p in primes(B):
            if p < 5: continue
            r = row(p); n += 1
            if (r is None) != (p not in T) or (r is not None and tuple(map(list, r)) != tuple(T[p])): bad += 1; print('DIFFERENT at', p)
        print('primes from 5 below %d: %d; rows of class_good_table.txt that differ: %d' % (B, n, bad)); sys.exit(1 if bad else 0)
    B1, B3 = 4400, 1200; rows = {}
    for p in primes(B1):
        if p >= 5 and (p % 4 == 1 or p <= B3):
            r = row(p)
            if r: rows[p] = r
    open(os.path.join(HERE, 'class_good_table.txt'), 'w').write(HEAD + ''.join(fmt(p, rows[p]) for p in sorted(rows)))
    r1 = [p for p in rows if p % 4 == 1]; r3 = [p for p in rows if p % 4 == 3]
    print('p = 1 (mod 4), up to %d: %d primes with a row, the last one %d; ratios mod p: no point %d, only residues %d, only non-residues %d' % (B1, len(r1), max(r1), sum(len(rows[p][0]) for p in r1), sum(len(rows[p][1]) for p in r1), sum(len(rows[p][2]) for p in r1)))
    print('p = 3 (mod 4), up to %d: %d primes with a row, the last one %d; ratios mod p without a point: %d' % (B3, len(r3), max(r3), sum(len(rows[p][0]) for p in r3)))
    print('the primes: %s' % ' '.join(map(str, sorted(rows))))
    # check 1: the local test of lean_cert.py modulo p^2 (points near the zeros of the forms included)
    sys.path.insert(0, os.path.join(HERE, '..')); import lean_cert as L
    tot = bad = 0
    for p in [13, 17, 29, 37, 41]:
        t = table(p); nr = next(a for a in range(2, p) if pow(a, (p - 1)//2, p) == p - 1); M = p * p
        for m in t:
            for mm, nn in ([(m, 1), (m + p, 1 + 2 * p)] if p < 30 else [(m, 1)]):
                surv = L.survivors(mm, nn, M, p); tot += 1
                bad += ((not L.killed(1, surv, M), not L.killed(nr, surv, M)) != t[m])
    print('check 1, against the test modulo p^2 of lean_cert.py, p = 13, 17, 29, 37, 41: cases %d, different %d' % (tot, bad))
    # check 2: exact p-adic classes (the function cls of class_local.py) at integer points, for integer ratios
    src = open(os.path.join(HERE, 'class_local.py')).read()
    i = src.index('def chi('); j = src.index('def cls('); k = j + re.search(r'\n(?=\S)', src[j:]).start()
    ns = {}; exec(src[i:k], ns); cls = ns['cls']; random.seed(7); tot = bad = 0
    for p in sorted(r1):
        full = table(p); free = [t for t in full if full[t] == (True, True)]
        picks = [(kd, t) for i, kd in enumerate('ERN') for t in random.sample(rows[p][i], min(2, len(rows[p][i])))] + [('B', t) for t in random.sample(free, min(2, len(free)))]
        for kind, tau in picks:
            m = tau + p * random.randrange(1, 50); n = 1 + p * random.randrange(1, 50)
            forms = lambda u: [m*m + (n*u)**2, n*n + (m*u)**2, ((m + n)*(u - 1))**2 + ((m - n)*(u + 1))**2, ((m - n)*(u - 1))**2 + ((m + n)*(u + 1))**2]
            us = [(u, 1) for u in range(p)] + [(1, 0)]
            for u0 in range(p):
                if any(f % p == 0 for f in forms(u0)):
                    for kk in range(p):
                        if any(f % (p*p) == 0 for f in forms(u0 + kk*p)): us += [(u0 + kk*p + p*p*s, 1) for s in range(1, p)]
            seen = set(cls(m, n, u, v, p) for u, v in us) - {None}; tot += 1
            bad += (seen != {'E': set(), 'R': {'1'}, 'N': {'e'}, 'B': {'1', 'e'}}[kind])
    print('check 2, exact classes at integer points (all residues of u, and points with a norm divisible by p^2), for the %d primes = 1 (mod 4) of the table: cases %d, different %d' % (len(r1), tot, bad))
