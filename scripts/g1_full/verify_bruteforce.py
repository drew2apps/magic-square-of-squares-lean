# Independent brute force: all primitive (a, b) modulo M, plain Python, no projective reduction.
import json, math, random
d = json.load(open('local5_200.json'))
def brute(m, n, M, p):
    sq = {x*x % M for x in range(M)}
    for a in range(M):
        for b in range(M):
            if a % p == 0 and b % p == 0: continue
            g4 = 2*(m*m+n*n)*(a*a+b*b)
            if all(v % M in sq for v in (m*m*b*b+n*n*a*a, n*n*b*b+m*m*a*a, g4, g4-8*m*n*a*b, g4+8*m*n*a*b)): return False
    return True
random.seed(1); bad = 0; tested = 0
items = list(d['closed'].items()); random.shuffle(items)
for key, rs in items[:400]:
    m, n = map(int, key.split(','))
    for r in rs[:2]:
        if r.startswith('p='): M = p = int(r[2:])
        elif r.startswith('2^'): p = 2; M = 2**int(r[2:])
        elif '^' in r: p, k = map(int, r.split('^')); M = p**k
        else: continue
        if M > 700: continue
        tested += 1
        if not brute(m, n, M, p): bad += 1; print('MISMATCH', m, n, r)
# negative control: pairs with a trivial point must never be obstructed
neg = 0
for m, n in d['open'][:20]:
    for (M, p) in [(32,2),(7,7),(17,17),(23,23),(31,31),(41,41),(47,47),(49,7),(25,5),(9,3)]:
        if brute(m, n, M, p): neg += 1; print('FALSE OBSTRUCTION', m, n, M)
print('closures re-checked by brute force:', tested, 'mismatches:', bad, '| false obstructions on trivial-point pairs:', neg)
