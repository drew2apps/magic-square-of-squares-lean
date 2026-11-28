#!/usr/bin/env python3
"""2-Selmer group of E' = E/<T3>: y^2 = x (x - e2')(x - e3') over K=Q(i) (same code as selmer.py)."""
import sys, json
sys.path.insert(0, '/mnt/user-data/outputs/r16/scripts/second_descent')
from sage.all__sagemath_schemes import *
import kdescent as kd
gd = kd.gd
K = QuadraticField(-1, 'i'); I = K.gen()
eg = [(0, 0), (-63145210128, -36055447954), (-896553529872, 511925434546)]
# Selmer group (copy of kd.selmer_bound, returning the group)
import random as pyrandom
rng = pyrandom.Random(1); primes = set([2])
for a in range(3):
    for b in range(a + 1, 3): primes |= kd.gauss_factor(kd.gsub(eg[a], eg[b]))
primes = sorted(primes); basis = [(0, 1), (1, 1)]
for q in primes:
    if q == 2: continue
    if q % 4 == 3: basis.append((q, 0))
    else:
        g = gd.gprime(q); basis += [g, (g[0], -g[1])]
g = len(basis); cons = []
for q in primes:
    for pl in kd.places_above(q):
        img = kd.local_image(eg, pl, rng); assert img is not None
        d = pl.dim; M = [pl.vec(b) for b in basis]
        for phi in kd.nullspace(img, 2*d):
            row = 0
            for k in range(g):
                a = sum((phi >> t & 1) & M[k][t] for t in range(d)) & 1
                c = sum((phi >> (d + t) & 1) & M[k][t] for t in range(d)) & 1
                row |= a << k | c << (g + k)
            cons.append(row)
sel = kd.nullspace(cons, 2*g)
print('primes', primes); print('basis', basis); print('dim Sel', len(sel))
def elt(v):
    b1 = K(1); b2 = K(1)
    for k in range(g):
        z = basis[k][0] + I*basis[k][1]
        if v >> k & 1: b1 *= z
        if v >> (g + k) & 1: b2 *= z
    return b1, b2
out = []
for mask in range(1 << len(sel)):
    v = 0
    for k in range(len(sel)):
        if mask >> k & 1: v ^= sel[k]
    b1, b2 = elt(v); b3 = elt((v ^ (v >> g)) & ((1 << g) - 1))[0]   # b3 = b1*b2 mod squares
    out.append([[int(b1[0]), int(b1[1])], [int(b2[0]), int(b2[1])], [int(b3[0]), int(b3[1])]])
    print(mask, b1.factor(), '|', b2.factor(), '|', b3.factor())
json.dump(dict(e=eg, selmer=out), open('selmer_Eprime.json', 'w'))
