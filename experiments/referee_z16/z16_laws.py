"""Reviewer z16 (an AI model; not read by a person): EXPLORATION_5 §3, laws (a) and (b), recomputed.

(a) From EXPLORATION_4's tables (data/out_*.jsonl: the 6,904 residual 4-tuples at the 937 ratios) and the repository's
    lists of the 9,736 stage-1 ratios (roadmap-g/g2/g2_997_all.txt, g2_15000_all.txt).  Shapes: S1 = (mn d0 d8,
    (m^2-n^2) dL dU) ~ (2mn, 1), S2 ~ (1, 2(m^2-n^2)) up to squares.
(b) Membership of (1, 2mn, 2mn) in Sel^2 of E_08: y^2 = x(x - A^2)(x - C^2), (A, B, C) = (m^2-n^2, 2mn, m^2+n^2),
    with MY OWN local images (random p-adic points, exact rational arithmetic; the image is the span of the classes
    found, checked to have the size 4 (odd p) or 8 (p = 2) that |E(Q_p)/2E(Q_p)| has for full 2-torsion), on a fresh
    sample (seed 2016, 3,000 ratios m < 4000, C not a square, as law_b.py) and on the session's rule
       opposite parity: in Sel^2  <=>  no prime = 5 mod 8 divides m^2 - n^2;   both odd: never.
    Also the proof sketch's 2-adic gap: the ratios with mn = 2 mod 4 (hypothesis H2 of note2 fails there).
Deterministic.  usage: python3 -I z16_laws.py [quick]"""
import sys, os, json, random, collections, math
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
from sympy import factorint

def sqf(x):
    return math.prod(q for q, e in factorint(abs(x)).items() if e % 2)

# ---------------- (a)
def law_a():
    rows = []
    for f in ('out_997.jsonl', 'out_15000_0.jsonl', 'out_15000_1.jsonl'):
        rows += [json.loads(l) for l in open(os.path.join(HERE, 'data', f))]
    S1, S2 = {}, {}
    for r in rows:
        m, n = r['m'], r['n']
        for t, k in r['tuples']:
            a = sqf(m*n*t[0]*t[1]); b = sqf((m*m - n*n)*t[2]*t[3])
            if (a, b) == (sqf(2*m*n), 1): S1.setdefault((m, n), []).append(t)
            if (a, b) == (1, sqf(2*(m*m - n*n))): S2.setdefault((m, n), []).append(t)
    allr = []
    for f in ('g2_997_all.txt', 'g2_15000_all.txt'):
        for l in open(os.path.join(REPO, 'experiments', 'roadmap-g', 'g2', f)):
            if l.strip() and l[0].isdigit(): allr.append(tuple(map(int, l.split()[0].split(','))))
    print('(a) class-keeping four-form ratios (the 9,736 of g2_*_all.txt, called stage-1 by the session): %d (distinct %d); ratios with residual tuples: %d; S1-shape at %d ratios, S2-shape at %d, both %d'
          % (len(allr), len(set(allr)), len(rows), len(S1), len(S2), len(set(S1) & set(S2))))
    def oddexp_ok(x): return all(q % 8 in (1, 7) for q, e in factorint(x).items() if q > 2 and e % 2)
    def even_bad(x): return sorted({q % 8 for q, e in factorint(x).items() if q > 2 and not e % 2 and q % 8 not in (1, 7)})
    for name, S, X in (('S1', S1, lambda m, n: (m - n)*(m + n)), ('S2', S2, lambda m, n: m*n)):
        good = sum(oddexp_ok(X(m, n)) for m, n in S)
        suff = sum(1 for m, n in set(allr) if (m, n) not in S and oddexp_ok(X(m, n)))
        ev = collections.Counter(tuple(even_bad(X(m, n))) for m, n in S)
        par = collections.Counter('opposite' if (m + n) % 2 else 'both odd' for m, n in S)
        print('    %s: every odd prime of odd exponent in %s is +-1 mod 8 at %d of %d; condition holds without such a tuple at %d class-keeping'
              'stage-1 ratios; residues mod 8 of the primes not +-1 mod 8 that occur (to even power): %s; parity %s'
              % (name, 'm^2-n^2' if name == 'S1' else 'mn', good, len(S), suff, dict(ev), dict(par)))
    return S1, S2

# ---------------- (b)
def vcls(x, p):
    """Class of a nonzero rational x in Q_p^*/squares: (v mod 2, unit class) with unit class = Legendre bit (p odd)
    or unit mod 8 (p = 2)."""
    x = Fr(x); a, b = x.numerator, x.denominator; v = 0
    while a % p == 0: a //= p; v += 1
    while b % p == 0: b //= p; v -= 1
    u = a * b                                          # a/b = a*b / b^2
    if p == 2: return (v % 2, u % 8)
    return (v % 2, 0 if pow(u % p, (p - 1)//2, p) == 1 else 1)
def is_sq(x, p):
    c = vcls(x, p); return c == (0, 1) if p == 2 else c == (0, 0)
def cmul(c1, c2, p):
    if p == 2: return ((c1[0] + c2[0]) % 2, (c1[1] * c2[1]) % 8)
    return ((c1[0] + c2[0]) % 2, c1[1] ^ c2[1])

def local_image(e, p, rng, want):
    """Span of the classes (x - e1, x - e2, x - e3) of points of E(Q_p), from the 2-torsion and random points."""
    e1, e2, e3 = e
    found = set()
    for i in range(3):                                 # points of order 2: the i-th entry is the product of the others
        ei = e[i]; c = [None]*3
        for j in range(3):
            if j != i: c[j] = vcls(ei - e[j], p)
        c[i] = cmul(c[(i + 1) % 3], c[(i + 2) % 3], p); found.add(tuple(c))
    span = {tuple(((0, 1) if p == 2 else (0, 0)) for _ in range(3))}
    def add(c):
        nonlocal span
        if c in span: return
        span = span | {tuple(cmul(a, b, p) for a, b in zip(c, s)) for s in span}
    for c in found: add(c)
    tries = 0
    while len(span) < want and tries < 200000:
        tries += 1
        k = rng.randrange(-8, 9); base = rng.choice((0, e1, e2, e3))
        w = rng.randrange(1, p**6) * rng.choice((1, -1))
        if w % p == 0: continue
        x = base + (Fr(p)**k) * w
        vals = [x - e1, x - e2, x - e3]
        if any(v == 0 for v in vals): continue
        if not is_sq(vals[0]*vals[1]*vals[2], p): continue
        add(tuple(vcls(v, p) for v in vals))
    assert len(span) == want, ('local image incomplete', e, p, len(span))
    return span

def in_sel(m, n, rng):
    A, B, C = m*m - n*n, 2*m*n, m*m + n*n
    e = (0, A*A, C*C)
    cl = (1, sqf(B), sqf(B))
    primes = sorted({2} | set(factorint(A)) | set(factorint(B)) | set(factorint(C)))
    for p in primes:
        img = local_image(e, p, rng, 8 if p == 2 else 4)
        if tuple(vcls(x, p) for x in cl) not in img: return False, p
    return True, None                                  # real place: (1, beta, beta) with beta > 0 is the class of x > C^2

def law_b(N, seed):
    rng = random.Random(seed); rows = []
    while len(rows) < N:
        m = rng.randrange(2, 4000); n = rng.randrange(1, m)
        if math.gcd(m, n) != 1: continue
        C = m*m + n*n
        if math.isqrt(C)**2 == C: continue
        ins, at = in_sel(m, n, rng)
        five = any(q % 8 == 5 for q in factorint(m*m - n*n))
        par = 'opposite' if (m + n) % 2 else 'both odd'
        h2 = 'mn=2 mod 4' if (m*n) % 4 == 2 else ('4|mn' if par == 'opposite' else '-')
        rows.append((m, n, ins, five, par, h2, at))
    c = collections.Counter((r[4], r[3], r[2]) for r in rows)
    print('(b) %d random coprime ratios (seed %d), m < 4000, C not a square: (parity, a prime = 5 mod 8 divides m^2-n^2, in Sel^2): %s'
          % (N, seed, sorted(c.items())))
    rule = sum(r[2] == (r[4] == 'opposite' and not r[3]) for r in rows)
    print('    rule "opposite parity and no prime = 5 mod 8 divides m^2 - n^2": agrees on %d of %d' % (rule, len(rows)))
    c2 = collections.Counter((r[5], r[2]) for r in rows if r[4] == 'opposite' and not r[3])
    print('    opposite parity, no prime = 5 mod 8: by 2-adic type (H2 fails when mn = 2 mod 4): %s' % sorted(c2.items()))
    c3 = collections.Counter((r[4], r[6] if r[6] is None else ('2' if r[6] == 2 else str(r[6] % 8) + ' mod 8')) for r in rows if not r[2])
    print('    where membership fails (first prime found): %s' % sorted(c3.items(), key=str))

if __name__ == '__main__':
    quick = 'quick' in sys.argv
    law_a()
    law_b(300 if quick else 3000, 2016)
