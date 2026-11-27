"""Review x6 (run 6), written by an AI model; not read by a person.  Local labels at explicit points, computed here
without exact_local or tuple_descent.Local: own square root of -1 mod p^60 (PARI and Newton), own valuations and symbols;
points at random in both charts and t* + p^e w (e = 1..11) around every p-adic root of every component (p = 1 mod 4) or
the real part of the Gaussian root (p = 3 mod 4); at 2, (t, 1), (1, 2t) for t < 2^12 and random points, with the
repository's key2.  witness(...) returns a point whose labels are the tuple's (up to a unit), or None."""
import random, sys, os
from math import gcd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "roadmap-g", "g1"))
import cypari2
P = cypari2.Pari()
import exact_local as X
TD = X.TD
def iota_mod(p, N):
    r = int(P.lift(P.sqrt(P.Mod(-1, p))))
    M = p**N
    for _ in range(12 + N.bit_length()): r = (r - (r*r + 1)*pow(2*r, -1, M)) % M
    assert (r*r + 1) % M == 0
    return r, M

def vp(x, p):
    if x == 0: return None
    v = 0
    while x % p == 0: x //= p; v += 1
    return v

def leg(a, p):
    a %= p
    return 0 if pow(a, (p - 1)//2, p) == 1 else 1

def label_point(forms, p, u, v, io=None):
    """labels of the four forms at the integer point (u, v), my own code; None if some class is not rational."""
    out = []
    for al, be in forms:
        z0, z1 = al[0]*u + be[0]*v, al[1]*u + be[1]*v
        if p % 4 == 3:
            w = min(x for x in (vp(z0, p), vp(z1, p)) if x is not None)
            a, b = z0 // p**w, z1 // p**w
            if leg(a*a + b*b, p): return None
            out.append((w % 2, 0))
        else:
            r, M = io
            labs = []
            for s in (1, -1):
                c = (z0 + s*r*z1) % M; w = vp(c, p)
                assert w is not None and w < 40, 'precision'
                labs.append((w % 2, leg(c // p**w, p)))
            if labs[0] != labs[1]: return None
            out.append(labs[0])
    return tuple(out)

CAP = 3000

def points(forms, p, rnd, n_rand=1500):
    N = 60; M = p**N; io = iota_mod(p, N) if p % 4 == 1 else None
    pts = []
    if p == 2:
        pts = [(t, 1) for t in range(1, 2**12)] + [(1, 2*t) for t in range(2**12)] + [(rnd.randrange(-2**40, 2**40), rnd.randrange(1, 2**40)) for _ in range(4000)]
    else:
        for _ in range(n_rand):
            pts.append((rnd.randrange(-p**3, p**3), rnd.randrange(1, p**3)))
            pts.append((rnd.randrange(1, p**3), p**rnd.randint(1, 4)*rnd.randrange(1, p**2)))
        cen = [(0, 0), (1, 0)]
        for al, be in forms:
            if p % 4 == 1:
                r = io[0]
                for s in (1, -1):
                    a = (al[0] + s*r*al[1]) % M; b = (be[0] + s*r*be[1]) % M
                    va, vb = vp(a, p), vp(b, p)
                    if va is not None and (vb is None or va <= vb): cen.append((0, (-(b // p**va) * pow(a // p**va, -1, M)) % M))
                    if vb is not None and (va is None or vb <= va): cen.append((1, (-(a // p**vb) * pow(b // p**vb, -1, M)) % M))
            else:
                na = al[0]**2 + al[1]**2; num = (-(be[0]*al[0] + be[1]*al[1]), -(be[1]*al[0] - be[0]*al[1])); w = vp(na, p) or 0
                if all(vp(x, p) is None or vp(x, p) >= w for x in num): cen.append((0, (num[0] // p**w) * pow(na // p**w, -1, M) % M))
        for c, t0 in cen:
            for e in range(1, 12):
                for _ in range(30):
                    t = (t0 + p**e * rnd.randrange(1, p**2)) % M
                    if c == 0: pts.append((t if t < M//2 else t - M, 1))
                    elif t % p == 0: pts.append((1, t))
    found = {}
    for u, v in pts:
        g = gcd(u, v); u, v = u // g, v // g
        gs = [(al[0]*u + be[0]*v, al[1]*u + be[1]*v) for al, be in forms]
        if any(z == (0, 0) for z in gs): continue
        if p == 2:
            ks = tuple(TD.key2(z)[1] for z in gs)
            if not all(k in TD.D2 for k in ks): continue
            L = TD.Local(forms, 2, 1)
            for kk in (ks, tuple(L.times3(k) for k in ks)): found.setdefault(kk, (u, v))
        else:
            lab = label_point(forms, p, u, v, io)
            if lab is None: continue
            found.setdefault(lab, (u, v))
            if p % 4 == 1: found.setdefault(tuple((l[0], 1 - l[1]) for l in lab), ('unit', u, v))
    return found

def dlab(d, p):
    if p == 2: return TD.key2((d, 0))[1]
    v = vp(d, p); return (v % 2, leg(d // p**v, p) if p % 4 == 1 else 0)

def divisors(sup):
    out = []
    for S in sup:
        ds = [1]
        for q in S: ds += [x*q for x in ds]
        out.append(ds)
    return out

WIT = {}
def witness(forms, key, p, tup, rnd):
    if (key, p) not in WIT: WIT[(key, p)] = points(forms, p, rnd)
    return WIT[(key, p)].get(tuple(dlab(d, p) for d in tup))

