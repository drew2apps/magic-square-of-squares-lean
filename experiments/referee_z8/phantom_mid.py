"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
H1 tested good primes to 20,000 and leaned on H3 (wrong) above 1153.  This tests the known and extra tuples of the 24
phantom ratios of H1 at every good prime p = 1 mod 4 with 17 <= p < 1200 and 20000 < p < 67583 with my_local.allowed
(label vector of a tuple: (chi(d_j / p))).  Known tuples must be allowed (asserted).  67,583 is the least p with
(p - 511)^2 > 66564 p: from there Hasse-Weil (quoted) on the genus-129 cover that carries a label vector gives a point
outside the 8 * 64 points over the branch points, so no good prime p >= 67,583 removes anything.  (Good p = 3 mod 4 keep
everything for these ratios: they have rational points.)  Rests on my_local.py.
usage: phantom_mid.py [--write]   (phantom_mid.txt; about 6 minutes)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit, isprime, good, load_phantoms
import my_local as M
p = 512
while not ((p - 511) ** 2 > 66564 * p): p += 1
out = ['z8: H1\'s 24 ratios at good primes p = 1 mod 4 (my local test). Written by an AI model; not read by a person.',
       'least p with (p - 511)^2 > 66564 p: %d' % p]
R = load_phantoms()
def lab(T, p): return sum(((0 if pow(d % p, (p - 1) // 2, p) == 1 else 1) << j) for j, d in enumerate(T))
for lo, hi in ((16, 1200), (20000, 67583)):
    ps = [p for p in range(lo + 1, hi) if p % 4 == 1 and isprime(p)]
    cuts = []; nchk = 0
    for p in ps:
        ch = M.chitab(p)
        for (m, n), (known, extra) in sorted(R.items()):
            if not good(m, n, p): continue
            a = M.allowed(m * pow(n, -1, p) % p, 1, p, ch); nchk += 1
            for T in known: assert a >> lab(T, p) & 1, ('a known tuple is not allowed', m, n, p, T)
            for T in extra:
                if not a >> lab(T, p) & 1: cuts.append(((m, n), p, T))
    out.append('%d < p < %d: %d primes, %d (ratio, good prime) pairs, extra tuples cut: %d %s' % (lo, hi, len(ps), nchk, len(cuts), cuts[:5]))
emit(out, 'phantom_mid.txt')
