"""Reviewer z16 (an AI model; not read by a person): soundness tests of the five-form local test, with my own class code.

A rational solution has classes (d0, d8, dL, dU, d4) with d_j gamma_j(u, v) a square in Q(i); at every prime p its
label tuple is the label tuple of the p-adic point (u, v).  So the property an exclusion needs is: at every prime p,
every primitive (u, v) in Z_p^2 at which all five gamma_j have the class of a rational number gives a label tuple in
exact_local.allowed_at(forms5, p).  Tested here three ways, none using tuple_descent.Local.cls:

  A. planted rational points: the degenerate solutions t = +-1 (m^2+n^2 a square) and t = 0, infinity (2(m^2+n^2) a
     square), classes by Lemma S1 done in integers (my_class); every one must pass at every bad prime and at the good
     primes below 60, and survive the join of z16_five.residual (end to end, support of d4 included);
  B. p-adic points built directly: integers (u, v) chosen at random, near every root of every component of every form
     (p = 1 mod 4: t = t* + p^j r in both charts), with prescribed 2-adic shapes at p = 2; a point is used when my own
     local class (my_label: valuations and Legendre symbols, or at 2 squares of units modulo 8) says all five forms
     have a rational class.  Its label tuple must be in the allowed set;
  C. my_label against tuple_descent.Local.cls at every point of B (agreement of the two class computations).
Seeds fixed; output deterministic.  usage: python3 -I z16_witness.py [quick]"""
import sys, os, random, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import z16_five as Z
X, TD, fact = Z.X, Z.TD, Z.fact

def sqf(x):
    return math.prod(q for q, e in fact(x).items() if e % 2)

def my_class(z):
    """Lemma S1 in integers: z = (x, y) with x^2 + y^2 = w^2; z = r * square with r = 2(x + w) (or 2(x - w) if that is
    0); the class is the squarefree part of |r| (-1 = i^2)."""
    x, y = z; N = x*x + y*y; w = math.isqrt(N); assert w*w == N and N > 0
    r = 2*(x + w) if x + w != 0 else 2*(x - w)
    return sqf(abs(r))

def vp(x, p):
    v = 0
    while x % p == 0: x //= p; v += 1
    return v, x

def chi(a, p):                                  # 0 = square (the repository's convention), 1 = non-square
    return 0 if pow(a % p, (p - 1)//2, p) == 1 else 1

SQU8 = {((a*a - b*b) % 8, (2*a*b) % 8) for a in range(8) for b in range(8) if (a + b) % 2}
def my_label(z, p, iota=None, M=None):
    """None if some gamma has no rational class at p (or precision is short); else the rational class as
    p = 2: d in {1, 2, 3, 6}; p odd: (v mod 2, chi(unit)) for p = 1 mod 4, (v mod 2, 0) for p = 3 mod 4."""
    x, y = z
    if p == 2:
        e = 0
        while (x + y) % 2 == 0: x, y = (x + y)//2, (y - x)//2; e += 1         # divide by 1 + i
        hits = []
        for d in (1, 2, 3, 6):
            a, b, ee = d*x, d*y, e
            while (a + b) % 2 == 0: a, b = (a + b)//2, (b - a)//2; ee += 1
            if ee % 2 == 0 and (a % 8, b % 8) in SQU8: hits.append(d)            # (1+i)^2 = 2i is a square times i... see below
        # (1+i)^(2k) = (2i)^k and 2i = (1+i)^2, so an even power of 1+i is a square: only the unit part matters
        return hits[0] if len(hits) == 1 else (None if not hits else ('MANY', tuple(hits)))
    if p % 4 == 3:
        v = 0
        while x % p == 0 and y % p == 0: x //= p; y //= p; v += 1
        if chi(x*x + y*y, p): return None
        return (v % 2, 0)
    out = []
    for s in (1, -1):
        c = (x + s*iota*y) % M
        if c == 0: return None
        v, u = vp(c, p)
        if v > 40: return None
        out.append((v % 2, chi(u, p)))
    return out[0] if out[0] == out[1] else None

def lab_to_repo(L, p, lab):
    return L.dlabel(lab) if p == 2 else lab

def sqrt_m1(p, N):
    g = 2
    while pow(g, (p - 1)//2, p) != p - 1: g += 1
    r = pow(g, (p - 1)//4, p); M = p**N
    for _ in range(N.bit_length() + 2): r = (r - (r*r + 1)*pow(2*r, -1, M)) % M
    assert (r*r + 1) % M == 0
    return r, M

def forms5(m, n):
    forms, sup, known = TD.magic(m, n)
    return forms + [Z.gamma4(m, n)]

def val(f, u, v): return (f[0][0]*u + f[1][0]*v, f[0][1]*u + f[1][1]*v)

def test_B(m, n, p, rng, n_random=1500, n_near=40, depth=8):
    F = forms5(m, n); L, al = X.allowed_at(F, p)
    iota = M = None
    if p % 4 == 1: iota, M = sqrt_m1(p, 60)
    pts = []
    for _ in range(n_random):
        u, v = rng.randrange(-10**9, 10**9), rng.randrange(1, 10**9)
        pts.append((u, v))
    if p % 4 == 1:
        for A, B in F:
            for s in (1, -1):
                a = (A[0] + s*iota*A[1]) % M; b = (B[0] + s*iota*B[1]) % M
                k = min(vp(a, p)[0] if a else 50, vp(b, p)[0] if b else 50)    # a common factor p^k (gamma4 at p | m^2+n^2)
                a, b, M2 = a // p**k, b // p**k, M // p**k
                # chart (t, 1): a t + b = 0; chart (1, t): a + b t = 0 (needs a unit coefficient)
                for chart, (c0, c1) in ((0, (a, b)), (1, (b, a))):
                    if c0 % p == 0: continue
                    root = (-c1 * pow(c0, -1, M2)) % M2
                    for j in range(1, depth + 1):
                        for _ in range(n_near):
                            t = (root + p**j * rng.randrange(1, 10**6)) % (p**(j + 30))
                            pts.append((t, 1) if chart == 0 else (1, t))
    if p == 2:
        for _ in range(n_random):
            a, b = rng.randrange(0, 6), rng.randrange(0, 6)
            u = 2**a * rng.randrange(1, 10**6) * rng.choice((1, -1)); v = 2**b * rng.randrange(1, 10**6)
            g = math.gcd(u, v); pts.append((u//g, v//g))
    tested = fails = disagree = 0; first = None
    for u, v in pts:
        if math.gcd(u, v) != 1: continue
        labs = [my_label(val(f, u, v), p, iota, M) for f in F]
        if any(l is None for l in labs): continue
        assert not any(isinstance(l, tuple) and l and l[0] == 'MANY' for l in labs), labs
        tup = tuple(lab_to_repo(L, p, l) for l in labs)
        tested += 1
        if tup not in al:
            fails += 1; first = first or (u, v, tup)
        # C: agreement with the repository's class code (exact values, huge disc radius)
        rl = [L.cls(val(f, u, v), 10**6) for f in F] if p != 2 else [TD.Local(F, 2, 12).cls(val(f, u, v), 10**6) for f in F]
        if any(r[0] or not r[1] for r in rl) or tuple(r[2] for r in rl) != tup: disagree += 1
    return tested, fails, disagree, first

def test_A(m, n, extra_good=60):
    C = m*m + n*n; F = forms5(m, n)
    pts = []
    if math.isqrt(C)**2 == C: pts += [(1, 1), (-1, 1)]
    if math.isqrt(2*C)**2 == 2*C: pts += [(0, 1), (1, 0)]
    tups = []
    for u, v in pts:
        vals = [val(f, u, v) for f in F]
        for z in vals: assert math.isqrt(z[0]**2 + z[1]**2)**2 == z[0]**2 + z[1]**2     # a rational point
        tups.append(tuple(my_class(z) for z in vals))
    primes = sorted(set(Z.bad_primes(m, n)) | {q for q in range(3, extra_good) if all(q % r for r in range(2, q))})
    nprime = nfail = 0
    for p in primes:
        L, al = X.allowed_at(F, p)
        for t in tups:
            nprime += 1; nfail += tuple(L.dlabel(d) for d in t) not in al
    T5 = Z.residual(m, n)
    lost = [t for t in tups if t not in T5]
    return len(pts), nprime, nfail, len(T5), lost

def main():
    quick = 'quick' in sys.argv
    rng = random.Random(16)
    # A: every coprime degenerate ratio m <= 200 (quick: <= 60), and 40 random degenerate ratios with m up to 10^5
    degs = [(m, n) for m in range(2, 61 if quick else 201) for n in range(1, m) if math.gcd(m, n) == 1 and
            (math.isqrt(m*m + n*n)**2 == m*m + n*n or math.isqrt(2*(m*m + n*n))**2 == 2*(m*m + n*n))]
    r2 = random.Random(7)
    big = set()
    while len(big) < (8 if quick else 40):
        a, b = r2.randrange(2, 300), r2.randrange(1, 300)
        if a <= b or math.gcd(a, b) != 1 or (a - b) % 2 == 0: continue
        m, n = a*a - b*b, 2*a*b                    # m^2 + n^2 a square
        if r2.random() < 0.5: m, n = m + n, abs(m - n)   # its partner: m^2 + n^2 twice a square
        m, n = max(m, n), min(m, n); g = math.gcd(m, n); m, n = m//g, n//g
        if m > n > 0: big.add((m, n))
    tot = [0, 0, 0, 0]; lostany = []
    for m, n in degs + sorted(big):
        npts, nprime, nfail, n5, lost = test_A(m, n)
        tot[0] += 1; tot[1] += npts; tot[2] += nprime; tot[3] += nfail
        if lost: lostany.append((m, n, lost))
    print('A. planted degenerate solutions: %d ratios (%d with m <= %d, %d larger), %d points; %d (point, prime) checks, %d '
          'refused; refused by the join: %d ratios' % (tot[0], len(degs), 60 if quick else 200, len(big), tot[1], tot[2], tot[3], len(lostany)))
    for x in lostany: print('   LOST', x)
    # B, C
    ratios = [(14903, 9696), (24599, 5207), (3320, 2303), (627, 364), (465, 416), (90384, 58537), (115513, 18217), (12, 5), (7, 1)]
    if quick: ratios = ratios[:3]
    T = [0, 0, 0]
    for m, n in ratios:
        ps = sorted(set(q for q in Z.bad_primes(m, n) if q < 10**7) | {3, 5, 7, 11, 13, 17, 29, 37, 41})
        line = []
        for p in ps:
            t, f, d, first = test_B(m, n, p, rng, n_random=300 if quick else 1500)
            T[0] += t; T[1] += f; T[2] += d; line.append('%d:%d/%d/%d' % (p, t, f, d))
            if first: print('   FAIL', (m, n), p, first)
        print('B/C (%d,%d) prime:points/refused/disagree  %s' % (m, n, ' '.join(line)))
    print('B. p-adic points with five rational classes: %d tested, %d refused by the allowed set' % (T[0], T[1]))
    print('C. my class code against tuple_descent.Local.cls at those points: %d disagreements' % T[2])

if __name__ == '__main__':
    main()
