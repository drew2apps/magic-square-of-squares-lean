"""Review x6 (run 6), written by an AI model; not read by a person.  The control (1109503,597671): arithmetic, the 6 tuples, explicit local points for each at 2 and at every odd bad
prime, and the good primes = 1 mod 4 from 20,000 to 600,000 (above 595,717 exact_local's Weil bound puts every
pattern on the generic residues of a good prime, so every tuple passes there; good primes = 3 mod 4 give every tuple the same labels)."""
import sys, os, time, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import mylocal as ML
X, TD, P = ML.X, ML.TD, ML.P
from math import isqrt, gcd
m, n = 1109503, 597671
S2 = m*m + n*n; s = isqrt(S2 // 2)
print('m^2 + n^2 = %d = 2 * %d^2: %s; %d = %s; gcd %d; m, n odd: %s' % (S2, s, 2*s*s == S2, s, P.factor(s), gcd(m, n), m % 2 == 1 and n % 2 == 1))
forms, sup, known = TD.magic(m, n)
bad = sorted(set(int(q) for x in (2*m*n*(m*m-n*n), S2, m**8+14*m**4*n**4+n**8) for q in P.factor(x)[0]))
print('bad primes:', bad)
print('known (t = 0, t = infinity):', known)
primes = sorted(set(bad) | set(int(q) for q in P.primes([2, 20000])))
t0 = time.time(); T, used = X.descend_exact(forms, sup, known, primes)
print('descent, every bad prime and good primes <= 20000: %d tuples: %s (%.0f s)' % (len(T), T, time.time() - t0))
rnd = random.Random(7)
for tup in T:
    miss = [q for q in bad if ML.witness(forms, (m, n), q, tup, rnd) is None]
    w2 = ML.witness(forms, (m, n), 2, tup, rnd)
    print('  %s: explicit point at 2: %s; at the odd bad primes: %s' % (tup, w2, 'all found' if not miss else 'MISSING at %s' % miss))
t0 = time.time(); fails = {}; cnt = 0
for q in P.primes([20001, 600000]):
    q = int(q)
    if q % 4 != 1 or q in bad: continue
    L, al = X.allowed_at(forms, q); cnt += 1
    for tup in T:
        if tuple(L.dlabel(d) for d in tup) not in al: fails.setdefault(tup, []).append(q)
print('good primes = 1 mod 4 in (20000, 600000]: %d tested (%.0f s); tuples failing somewhere: %s' % (cnt, time.time() - t0, fails or 'none'))
