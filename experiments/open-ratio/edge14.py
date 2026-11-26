"""Edge 1-4 in detail: Sel^phi(E) (X on E', E' = curve 4) and the images of the torsion of E' and of E."""
import sys; sys.path.insert(0, '.')
from phi_descent import *
a, b = EDGES['1-4']
ap = gsc(-2, a); bp = gadd(gd.gmul(a, a), gsc(-4, b))
primes, basis, sE, sEp = edge(a, b, 1)
def elem(v):
    z = (1, 0)
    for k, g in enumerate(basis):
        if v >> k & 1: z = gd.gmul(z, g)
    return z
print("a =", a, " b =", b); print("a' =", ap, " b' =", bp)
# sqrt(b) in Z[i]?
pb = pari(f"nfroots(nfinit(y^2+1), x^2 - ({b[0]} + {b[1]}*y))"); print("sqrt(b) in K:", pb)
print("Sel^phi(E) basis (X on E'):", [elem(v) for v in sEp])
print("Sel^phihat(E') basis (x on E):", [elem(v) for v in sE])
# torsion of E' and E via PARI
K = pari("bnfinit(y^2+1,1)")
for name, (A, B) in (("E'", (ap, bp)), ("E", (a, b))):
    Ee = pari.ellinit([0, f"{A[0]}+{A[1]}*y", 0, f"{B[0]}+{B[1]}*y", 0], K)
    T = pari.elltors(Ee); print(name, "torsion", T[0], T[1], "gens:", pari.lift(T[2]))
    pts = []
    # all torsion points: enumerate group
    gens = T[2]; ords = [int(o) for o in T[1]]
    import itertools
    for cs in itertools.product(*[range(o) for o in ords]):
        P = pari("[0]")
        for c, g in zip(cs, gens): P = pari.elladd(Ee, P, pari.ellmul(Ee, g, c))
        if P == pari("[0]"): continue
        X = pari.lift(P[0])
        pts.append((cs, X))
    for cs, X in pts:
        print("   ", cs, "x =", X)
