#!/usr/bin/env python3
"""Review x4 (run 6): a second, independent sample of PARI's ellap (as used by lvalue_stream.py --ap) against the repository's exact
order certificate lball_exact53.ap_one on the curve w*f0*f8 of (6557,324): random primes (seed 6, not the bundle's
seed 1), the 300 largest primes <= Mmax and 300 primes just above 1000.  Evidence about the condition, not a removal.  Written by an AI model; not read by a person.  python3 apsample.py [K]
(K random primes, default 30000; about 100 s).  Prints only."""
import os, sys, random, time
R = os.environ.get('MSS_REPO', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
K = int(sys.argv[1]) if len(sys.argv) > 1 else 30000
os.environ['MSS_PAIR'] = '6557,324'; os.environ['MSS_CURVE'] = 'w*f0*f8'
sys.path.insert(0, R + '/experiments/census-wide/lvalue')
import lvalue_big as LB
import cypari2
LX, L53 = LB.LX, LB.L53
pari = cypari2.Pari(); pari.allocatemem(10**9)
N, a, LD = LB.the_curve(); M = LX.max_terms(N)
Em = pari.ellinit(pari.ellminimalmodel(pari.ellinit(list(a))))
assert int(pari.ellglobalred(Em)[0]) == N
bad = {p for p, T in LD['rows'].items() if T['red'] != 'good'}
rng = random.Random(6); P = set()
while len(P) < K:
    p = int(pari.nextprime(rng.randrange(1000, M)))
    if p <= M and p not in bad: P.add(p)
q = M
for _ in range(300): q = int(pari.precprime(q)); P.add(q); q -= 1
q = 1000
for _ in range(300): q = int(pari.nextprime(q + 1)); P.add(q)
P = sorted(p for p in P if p not in bad)
st = L53.new_stats(); dis = []; t = time.time()
for p in P:
    x = L53.ap_one(a, LD, p, st); y = int(pari.ellap(Em, p))
    if x != y: dis.append((p, x, y))
print('N = %d, Mmax = %d; %d primes (%d random, seed 6; 300 largest <= Mmax; 300 above 1000): ap_one and PARI ellap disagree at %d; stats %s'
      % (N, M, len(P), K, len(dis), dict(sorted(st.items()))))
sys.stderr.write('[%.0f s]\n' % (time.time() - t))
print('disagreements:', dis[:10]); assert not dis
