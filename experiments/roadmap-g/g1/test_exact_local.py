"""Positive controls for exact_local.py (EXPERIMENT, written by an AI model; not read by a person).
  1. At small odd p, the exact set is contained in tuple_descent.Local's set at depth K (which can only be too large),
     and equals it once K is large enough (the depth used is printed).
  2. Random integer points (u, v), many of them close to the p-adic roots of the components: every point at which all
     four forms have the class of a rational number has its label tuple in the exact set.
  3. The descent with the exact test reproduces the repository's result on tuple_selmer_997.txt's first ratios
     (nothing left) and keeps the known points of two degenerate ratios.
usage: test_exact_local.py [--write]     (without --write, checks test_exact_local.txt)"""
import sys, os, random
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import exact_local as X
TD = X.TD

def point_labels(L, forms, u, v):
    p = L.p; lab = []
    for al, be in forms:
        z = (al[0]*u + be[0]*v, al[1]*u + be[1]*v)
        dead, det, l = L.cls(z, 10**6) if L.kind == 3 else exact_comp(L, z)
        if dead: return None
        lab.append(l)
    return tuple(lab)
def exact_comp(L, z):
    p = L.p; out = []
    for s in (1, -1):
        w = (z[0] + s*L.iota*z[1]) % L.M; v = X.vp(w, p); assert v is not None and v < L.N - 2
        out.append((v % 2, X.leg(w // p**v, p)))
    return out[0] != out[1], True, out[0]

def main():
    lines = []; random.seed(11)
    ratios = [(23, 7), (41, 9), (97, 56), (120, 119), (627, 364), (1189, 820), (4999, 1802), (9001, 1575), (35, 12), (21, 20), (47, 23), (161, 73)]
    for m, n in ratios:
        forms, sup, known = TD.magic(m, n)
        bad = sorted(q for q in set(x for P in sup for x in P) if 2 < q < 60)
        for p in sorted(set([3, 5, 7, 11, 13, 17, 29, 37, 41] + bad)):
            Lx, ex = X.allowed_at(forms, p)
            K = 2; prev = None
            while True:
                rep = TD.Local(forms, p, K).allowed()
                assert ex <= rep, ('exact not inside the repository set', m, n, p, K)
                if rep == ex or K >= (7 if p < 10 else 4): break
                K += 1
            lines.append('(%d,%d) p=%d: exact %d tuples; repository at depth %d: %d; %s' % (m, n, p, len(ex), K, len(rep), 'equal' if rep == ex else 'exact strictly smaller'))
    # 2. random points
    tested = alive = 0
    for m, n in ratios[:6]:
        forms, sup, known = TD.magic(m, n)
        for p in [5, 13, 17, 29, 3, 7, 11, 19, 1000003, 999983]:
            L, ex = X.allowed_at(forms, p)
            roots = []
            if L.kind == 1:
                for al, be in forms:
                    for s in (1, -1):
                        a, b = L.comp(be, s), L.comp(al, s)
                        if b % p: roots.append((-a * pow(b, -1, L.M)) % L.M)
            for i in range(4000):
                if roots and i % 2: t = (random.choice(roots) + p**random.randint(1, 6) * random.randrange(1, p)) % L.M; u, v = t, 1
                else: u, v = random.randrange(-10**12, 10**12), random.randrange(1, 10**12)
                if random.random() < 0.3: u, v = v, u * p**random.randint(1, 3)
                g = gcd(u, v); u, v = u // g, v // g
                lab = point_labels(L, forms, u, v); tested += 1
                if lab is None: continue
                alive += 1
                assert lab in ex, ('a point outside the exact set', m, n, p, u, v, lab)
    lines.append('random points: %d tested, %d with every form of rational class, all inside the exact set' % (tested, alive))
    # 3. descent
    for m, n in [(1189, 820), (4999, 1802), (35, 12), (21, 20), (120, 119), (47, 23), (161, 73)] + [tuple(int(x) for x in l.split()[0].strip('():').split(',')) for l in open(os.path.join(X.REPO, 'scripts', 'cover_descent', 'tuple_selmer_997.txt')) if l.startswith('(')][:6]:
        forms, sup, known = TD.magic(m, n)
        primes = sorted(set(q for P in sup for q in P) | set(TD.primerange(2, 151)))
        S, used = X.descend_exact(forms, sup, known, primes)
        Sr, usedr, sizes = TD.descend(forms, sup, known, 150)
        lines.append('(%d,%d): exact descent leaves %s; repository descent leaves %s' % (m, n, S if S else 'nothing', Sr if Sr else 'nothing'))
        assert set(S) <= set(Sr) and set(known) <= set(S)
    return lines

if __name__ == '__main__':
    out = '\n'.join(main()) + '\n'; f = os.path.join(HERE, 'test_exact_local.txt')
    if '--write' in sys.argv: open(f, 'w').write(out); print(out, end='')
    else:
        print(out, end='')
        if os.path.exists(f): assert open(f).read() == out, 'differs from test_exact_local.txt'; print('# matches test_exact_local.txt')
