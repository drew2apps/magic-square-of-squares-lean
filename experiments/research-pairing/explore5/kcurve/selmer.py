#!/usr/bin/env python3
"""Model of E over K=Q(i) for (m,n)=(627,364), quartic y^2 = L1*L3*M1*M2, and its 2-Selmer group (uses r16 kdescent.py)."""
import sys, json, random
sys.path.insert(0, '/mnt/user-data/outputs/r16/scripts/second_descent')
from sage.all__sagemath_schemes import *
import kdescent as kd
gd = kd.gd
K = QuadraticField(-1, 'i'); I = K.gen(); RK = PolynomialRing(K, 's'); s = RK.gen()
m, n = 627, 364; p = m + I*n; pb = m - I*n
Fs = [m - n*s, n + m*s, p*s - pb, pb*s + p]; Q = prod(Fs)
roots = [-f[0]/f[1] for f in Fs]; r = roots[3]
A = Q.derivative()(r)
xs = [A/(rk - r) for rk in roots[:3]]
den = lcm([x[j].denominator() for x in xs for j in (0, 1)])
e = [x*den*den for x in xs]
# check: Q(r + A*den^2/X) * X^4/(A*den)^2 == prod(X - e_k)
X = RK.gen()
lhs = sum(Q.derivative(k)(r)/factorial(k)*(A*den*den)**k*X**(4 - k) for k in range(1, 5))/(A*den)**2
assert lhs == prod(X - ek for ek in e), "model check failed"
d2 = e[1] - e[0]; d3 = e[2] - e[0]
OK = K.ring_of_integers(); gg = K.ideal(d2, d3).gens_reduced()[0]
c = prod(pp.gens_reduced()[0]**(ee//2) for pp, ee in K.ideal(gg).factor())
shift = e[0]; scale = c*c
e = [K(0), d2/scale, d3/scale]
print('shift', shift, 'scale c =', c, ' x_orig = (X*c^2 + shift)/den^2,  s = r + A/x_orig')
eg = [(int(z[0]), int(z[1])) for z in e]
print('Q(s) =', Q); print('r =', r, ' A =', A, ' den =', den); print('e =', eg)
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
json.dump(dict(m=m, n=n, e=eg, den=int(den), shift=[int(shift[0]), int(shift[1])], c=str(c), A=str(A), r=str(r), selmer=out), open('selmer.json', 'w'))
