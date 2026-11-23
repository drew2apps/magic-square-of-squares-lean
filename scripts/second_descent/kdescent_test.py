"""kdescent.py on curves defined over Q with all points of order 2 rational: y^2 = x(x + a)(x + b), a, b random integers.
There rank E(K) = rank E(Q) + rank E^(-1)(Q), and PARI gives bounds for both.  The Selmer bound over K must never be below
the rank; it can exceed it by an even number (the 2-part of the Tate-Shafarevich group over K).
usage: kdescent_test.py seed N"""
import random, sys
import cypari2
import kdescent as kd
pari = cypari2.Pari()
def bounds(a, b, tw):
    E = pari.ellinit([0, tw*(a + b), 0, a*b, 0])        # y^2 = x(x + tw a)(x + tw b)
    lo, hi = 0, 10**9
    for F in pari.ellisomat(pari.ellminimalmodel(E), 0, 1)[0]:
        r = pari.ellrank(pari.ellinit(F), 1); lo = max(lo, int(r[0])); hi = min(hi, int(r[1]))
    return lo, hi
rng = random.Random(int(sys.argv[1])); N = int(sys.argv[2]); bad = 0; eq = 0; sharp = 0; fail = 0; dist = {}
for _ in range(N):
    a = rng.randint(-60, 60); b = rng.randint(-60, 60)
    if a == 0 or b == 0 or a == b: continue
    kb = kd.selmer_bound([(0, 0), (-a, 0), (-b, 0)])
    if kb is None: fail += 1; continue
    l1, h1 = bounds(a, b, 1); l2, h2 = bounds(a, b, -1)
    if kb < l1 + l2: bad += 1; print('BELOW THE RANK', a, b, kb, (l1, h1), (l2, h2))
    if l1 == h1 and l2 == h2:
        sharp += 1; eq += kb == l1 + l2; dist[kb - l1 - l2] = dist.get(kb - l1 - l2, 0) + 1
print('curves y^2 = x(x+a)(x+b): local image not found: %d; Selmer bound over K below the PARI lower bound for rank E(Q) + rank E^(-1)(Q): %d' % (fail, bad))
print('curves with both ranks known: %d; bound over K minus rank over K: %s' % (sharp, sorted(dist.items())))
