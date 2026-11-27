"""Review x6 (run 6), written by an AI model; not read by a person.  The Weil step of exact_local: (1) the smallest p at
which weil_all decides every pattern, for each number of roots (p = 1 mod 4) and of quadratics (p = 3 mod 4); (2) on
593 (ratio, prime) cases, the allowed set with weil_all forced False (listing only) against the default.  Prints only.
usage: python3 -I weil_check.py"""
import sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "roadmap-g", "g1"))
import exact_local as X
TD = X.TD
import cypari2; P = cypari2.Pari()
from math import isqrt
orig = X.weil_all
def run(forms, p, force_brute):
    X.weil_all = (lambda p_, d, e: False) if force_brute else orig
    try: return X.allowed_at(forms, p)[1]
    finally: X.weil_all = orig
# 1. thresholds: smallest p with weil_all true, for g roots (kind 1) and g quadratics (kind 3)
for g in range(1, 9):
    lo = 3
    while not orig(lo, [1]*g, g): lo = lo*2
    a, b = lo//2, lo
    while a < b:
        c = (a + b)//2
        if orig(c, [1]*g, g): b = c
        else: a = c + 1
    print('kind 1, %d distinct roots: Weil decides every pattern from p >= %d' % (g, a))
for g in range(1, 5):
    for ex in (0, 4):
        a = 3
        while not orig(a, [2]*g, ex): a += 1
        print('kind 3, %d quadratics, %d excluded: from p >= %d' % (g, ex, a))
# 2. forced listing against the default, many (ratio, prime)
n_cases = 0; diffs = 0; t0 = time.time()
ratios = [(627, 364), (291, 224), (1109503, 597671), (993260025, 766111057), (10125, 4459), (4999, 1802), (1189, 820)]
for m, n in ratios:
    forms, sup, known = TD.magic(m, n)
    bad = sorted(set(int(q) for x in (m*n*(m*m-n*n), m*m+n*n, m**8+14*m**4*n**4+n**8) for q in P.factor(x)[0]))
    primes = [q for q in bad if 2 < q <= 5*10**6] + [int(q) for q in P.primes([3, 400])] + [int(P.nextprime(10**6 + 17*m % 1000)), int(P.nextprime(4*10**6 + n % 1000)), 1000003, 4999999 if P.isprime(4999999) else int(P.precprime(5*10**6))]
    for p in sorted(set(primes)):
        a, b = run(forms, p, False), run(forms, p, True)
        n_cases += 1
        if a != b:
            diffs += 1; print('DIFFERENCE', m, n, p, len(a), len(b), 'listing a subset' if b <= a else 'listing NOT a subset')
print('# %d (ratio, prime) cases with the Weil step and with listing forced: %d differences' % (n_cases, diffs))
