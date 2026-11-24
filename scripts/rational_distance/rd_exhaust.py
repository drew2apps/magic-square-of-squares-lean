#!/usr/bin/env python3
"""The first step of rd.py on a few lines, again by exhaustion (EXPERIMENT): tuple_bruteforce.py of ../cover_descent,
which shares no code with tuple_descent.py, on the four forms of the line y = r/s, at the primes up to 61 and the
primes of the supports above 61 (199 on the line 26/225).
For each line: the tuples kept by the exhaustion, the tuples left by the descent (primes up to 150), and whether every
tuple kept has four equal classes.  usage: rd_exhaust.py r,s [r,s ...]"""
import sys, os
from itertools import product
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'cover_descent'))
import tuple_bruteforce as bf
import rd
for arg in sys.argv[1:]:
    r, s = (int(x) for x in arg.split(','))
    forms = [((s, 0), (0, r)), ((s, 0), (-s, r)), ((s, 0), (0, r - s)), ((s, 0), (-s, r - s))]
    sup = [bf.factor(2*r*s), bf.factor(2*r*s), bf.factor(2*s*abs(r - s)), bf.factor(2*s*abs(r - s))]
    cand = [bf.products(P) for P in sup]; alive = None
    for p in sorted(set(bf.primes_to(61)) | set(q for P in sup for q in P)):
        N = 1
        while p**(N + 1) <= 5000: N += 1
        boxes = bf.kept(forms, cand, p, p**N)
        if alive is None:
            alive = set()
            for box in boxes: alive |= set(product(*box))
        else: alive = set(t for t in alive if any(all(t[j] in box[j] for j in range(4)) for box in boxes))
    T, mixed = rd.line(r, s)
    print('y = %d/%d: kept by the exhaustion: %d %s; left by the descent: %d; the descent set is inside the exhaustion set: %s; every tuple kept by the exhaustion has four equal classes: %s'
          % (r, s, len(alive), sorted(alive) if len(alive) <= 6 else '', len(T), set(T) <= alive, all(len(set(t)) == 1 for t in alive)), flush=True)
