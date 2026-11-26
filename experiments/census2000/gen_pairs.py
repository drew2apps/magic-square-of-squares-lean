"""Pairs 1 <= n < m, coprime, M0 < m <= M1, that survive Result 1b of G1_FULL_SYSTEM.md
(the 484 classes of modp_classes.txt and the 5 mod 8 condition). Experiment, written by an AI model (Claude).
usage: gen_pairs.py M0 M1 out.json"""
import sys, json, math, collections
R='/tmp/claude-0/-home-user/9ec7a8d3-42dc-53c8-8518-1a663ac809f1/scratchpad/repo/scripts/g1_full/'
cls = collections.defaultdict(set)
for line in open(R+'modp_classes.txt'):
    p, r = line.split()[:2]; cls[int(p)].add(int(r))
def odd5(N):
    x = N; d = 2
    while d * d <= x:
        if x % d == 0:
            k = 0
            while x % d == 0: x //= d; k += 1
            if d % 8 == 5 and k % 2: return True
        d += 1 if d == 2 else 2
    return x > 1 and x % 8 == 5
def inclass(m, n):
    return any(n % p and (m * pow(n, -1, p)) % p in rs for p, rs in cls.items())
def issq(x): r = math.isqrt(x); return r*r == x
M0, M1 = int(sys.argv[1]), int(sys.argv[2])
out = []
for m in range(max(2, M0+1), M1+1):
    for n in range(1, m):
        if math.gcd(m, n) == 1 and not odd5(m*m+n*n) and not inclass(m, n): out.append((m, n))
S = m2 = lambda m,n: m*m+n*n
typ = collections.Counter('square' if issq(m*m+n*n) else 'twice' if (m*m+n*n)%2==0 and issq((m*m+n*n)//2) else 'other' for m,n in out)
print(len(out), dict(typ))
json.dump({'open': out}, open(sys.argv[3], 'w'))
