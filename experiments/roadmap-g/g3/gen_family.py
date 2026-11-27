"""G3 (experiment, written by an AI model; not read by a person): the ratios m : n of the two degenerate types
(m^2 + n^2 a square, or twice a square), coprime, M0 < m <= M1, that pass the same congruences as
../../census2000/gen_pairs.py (the 484 classes and the 5 mod 8 condition; for these types the 5 mod 8 condition always holds).
Enumerates Pythagorean pairs directly instead of all pairs.  usage: gen_family.py M0 M1 out.json"""
import sys, os, json, math, collections
sys.argv_saved = list(sys.argv)
R = os.path.join(os.environ.get('MSS_REPO', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')), 'scripts', 'g1_full', '')
cls = collections.defaultdict(set)
for line in open(R + 'modp_classes.txt'):
    p, r = line.split()[:2]; cls[int(p)].add(int(r))
def inclass(m, n): return any(n % p and (m * pow(n, -1, p)) % p in rs for p, rs in cls.items())
def issq(x): r = math.isqrt(x); return r*r == x
M0, M1 = int(sys.argv[1]), int(sys.argv[2])
S = set()
a = 2
while a*a - 1 <= M1 * 2:
    for b in range(1, a):
        if (a - b) % 2 and math.gcd(a, b) == 1:
            x, y = a*a - b*b, 2*a*b; m, n = max(x, y), min(x, y)
            if M0 < m <= M1: S.add((m, n))                                   # square type
            s, d = (x + y), abs(x - y)                                     # twice-square type: (m+n)^2+(m-n)^2 = 2(m^2+n^2)
            m2, n2 = max(s, d), min(s, d)
            if M0 < m2 <= M1 and math.gcd(m2, n2) == 1: S.add((m2, n2))
    a += 1
out = sorted(p for p in S if not inclass(*p))
for m, n in out: assert issq(m*m + n*n) or ((m*m + n*n) % 2 == 0 and issq((m*m + n*n)//2))
typ = collections.Counter('square' if issq(m*m + n*n) else 'twice' for m, n in out)
print(len(S), 'degenerate coprime ratios;', len(out), 'pass the congruences', dict(typ))
json.dump({'open': out}, open(sys.argv[3], 'w'))
