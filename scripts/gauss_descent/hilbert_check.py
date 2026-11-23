#!/usr/bin/env python3
"""The class labels of gauss_descent.py against PARI's Hilbert symbols, on random Gaussian integers.

Two elements x, y of K = Q(i) are in the same class at the places above p exactly when x y is a square in every completion
K_P, P | p, and that holds exactly when the Hilbert symbol (x y, b)_P is 1 for all b in K_P*.  It is enough to take b from a
set that spans K_P*/K_P*^2; that the list used (B below, p, and the Gaussian primes above p) spans is checked by counting the distinct symbol vectors
(16 at 1+i, 4 at an inert prime, 4 at each of the two places of a split prime).
The program's labels (Local.cls: valuation parity and quadratic character; at 1+i the unit modulo pi^5) are written from
scratch; the symbols are PARI's nfhilbert.  They must agree: same label <=> all symbols 1.

usage: hilbert_check.py N [seed]      N random pairs for each prime"""
import sys, random
import cypari2
import gauss_descent as gd
pari = cypari2.Pari()
nf = pari('nfinit(y^2 + 1)')
def el(z): return pari('Mod(%d + %d*y, y^2 + 1)' % (z[0], z[1]))
B = [(0, 1), (1, 1), (3, 0), (5, 0), (1, 2), (2, 1), (7, 0), (2, 3), (3, 2), (1, 4), (11, 0), (13, 0), (4, 1), (2, 5)]
Bp = [el(b) for b in B]

def vec(z, P, Bl): return tuple(int(pari.nfhilbert(nf, el(z), b, P)) for b in Bl)

def rnd(p):
    while True:
        z = (random.randint(-60, 60), random.randint(-60, 60))
        if z == (0, 0): continue
        for _ in range(random.choice([0, 0, 0, 1, 1, 2, 3])):          # vary the valuations
            z = gd.gmul(z, random.choice([(p, 0), (1, 1)] if p == 2 else [(p, 0)] + ([gd.gprime(p), (gd.gprime(p)[0], -gd.gprime(p)[1])] if p % 4 == 1 else [])))
        return z

def main(N, seed):
    random.seed(seed); total = 0
    for p in [2, 3, 5, 7, 11, 13, 17, 29, 37]:
        L = gd.Local(p); Ps = list(pari.idealprimedec(nf, p)); seen = [set() for _ in Ps]; same = 0; und = 0
        Bl = Bp + [el((p, 0))] + ([el(gd.gprime(p)), el((gd.gprime(p)[0], -gd.gprime(p)[1]))] if p % 4 == 1 else [])     # uniformisers at every place above p
        for _ in range(N):
            x, y = rnd(p), rnd(p); cx, cy = L.cls(x), L.cls(y)
            if None in cx or None in cy: und += 1; continue
            xy = gd.gmul(x, y); vs = [vec(xy, P, Bl) for P in Ps]
            for s, P in zip(seen, Ps): s.add(vec(x, P, Bl))
            lab = cx == cy; sym = all(all(c == 1 for c in v) for v in vs)
            assert lab == sym, (p, x, y, cx, cy, vs)
            same += lab; total += 1
        want = 16 if p == 2 else 4
        assert all(len(s) == want for s in seen), (p, [len(s) for s in seen])
        print('p = %2d (%s): %d pairs, label equal in %d, and exactly then all symbols are 1; %d classes seen at each of %d place(s)' % (p, L.kind, N - und, same, want, len(Ps)), flush=True)
    print('%d pairs in all: the labels and the Hilbert symbols agree' % total)

if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]) if len(sys.argv) > 2 else 1)
