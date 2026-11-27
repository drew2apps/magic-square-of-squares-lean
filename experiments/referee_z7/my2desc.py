#!/usr/bin/env python3
"""Review z7 of roadmap H (H7) and of scripts/second_descent/rank0_ct.py: an independent complete 2-descent.
Written by an AI model (review agent z7, 27 November 2026); not read by a person.  Experiment, not counted.

Written from the textbook (Silverman, AEC X.1.4; Cremona, Algorithms for modular elliptic curves, 3.6) before the
reviewer read rank_second.py or the H7 scripts.  For y^2 = (x - e1)(x - e2)(x - e3), e_i distinct integers:
    Sel^2(E) = { d = (d1, d2) in Q(S,2)^2 : res_v(d) in W_v for every v in S u {inf} },
S = {-1} u primes of 2 (e1-e2)(e1-e3)(e2-e3), W_v = image of the local Kummer map P -> (x - e1, x - e2) (with the
usual values at the 2-torsion points).  W_v is the span of the images of the 2-torsion and of sampled points of
E(Q_v) (x rational with f(x) a square in Q_v, tested exactly), taken until it has the order |E(Q_v)/2E(Q_v)| =
|E(Q_v)[2]| |2|_v^{-1}: 4 (p odd), 8 (p = 2), 2 (v = inf; all roots real).  Exact arithmetic (Python integers and
Fractions); sympy factors the discriminant.  No PARI.
Rests on (quoted): the definition of the 2-Selmer group and the order of E(Q_v)/2E(Q_v).

    python3 my2desc.py        self-test on y^2 = x^3 - n^2 x (prints my2desc.txt)
"""
import sys, random
from fractions import Fraction as Fr
from sympy import factorint

def vp(a, p):
    a = int(a); v = 0
    while a % p == 0: a //= p; v += 1
    return v, a

def loc_class(q, p):
    """class of a non-zero rational q in Q_p*/Q_p*^2 as a tuple of bits; p = 0 is the real place"""
    q = Fr(q); assert q != 0
    if p == 0: return (1 if q < 0 else 0,)
    vn, un = vp(q.numerator, p); vd, ud = vp(q.denominator, p)
    v = vn - vd; u = un * ud
    if p == 2:
        u %= 8
        return (v % 2, 1 if u in (3, 7) else 0, 1 if u in (5, 7) else 0)
    leg = pow(u % p, (p - 1) // 2, p)
    return (v % 2, 0 if leg == 1 else 1)

def is_sq_local(q, p):
    q = Fr(q)
    if q == 0: return True
    return all(b == 0 for b in loc_class(q, p))

def span(vecs):
    """fully reduced F2 basis of bit vectors (ints), as (pivot, vector) pairs"""
    basis = []
    for v in vecs:
        for piv, b in basis:
            if v >> piv & 1: v ^= b
        if v:
            piv = v.bit_length() - 1
            basis = [(pp, (bb ^ v) if bb >> piv & 1 else bb) for pp, bb in basis]
            basis.append((piv, v))
    return basis

def tobits(t):
    r = 0
    for i, b in enumerate(t): r |= (b & 1) << i
    return r

class Curve:
    def __init__(self, e1, e2, e3):
        self.e = [int(e1), int(e2), int(e3)]
        assert len(set(self.e)) == 3
        e1, e2, e3 = self.e
        N = 2 * abs((e1 - e2) * (e1 - e3) * (e2 - e3))
        self.primes = sorted(int(p) for p in factorint(N).keys())
        self.S = [-1] + self.primes
    def f(self, x):
        e1, e2, e3 = self.e
        return (x - e1) * (x - e2) * (x - e3)
    def kummer(self, x):
        e1, e2, e3 = self.e
        if x == e1: return (Fr((e1 - e2) * (e1 - e3)), Fr(e1 - e2))
        if x == e2: return (Fr(e2 - e1), Fr((e2 - e1) * (e2 - e3)))
        return (Fr(x - e1), Fr(x - e2))
    def locvec(self, pair, p):
        return tobits(loc_class(pair[0], p) + loc_class(pair[1], p))
    def local_image(self, p, rng, tries=20000):
        e1, e2, e3 = self.e
        tdim = 1 if p == 0 else (3 if p == 2 else 2)
        vecs = [self.locvec(self.kummer(e1), p), self.locvec(self.kummer(e2), p)]
        B = span(vecs); es = sorted(self.e)
        for it in range(tries):
            if len(B) >= tdim: break
            if p == 0:
                if rng.randint(0, 1) == 0: x = Fr(es[0]) + Fr(rng.randint(1, 999), 1000) * (es[1] - es[0])
                else: x = Fr(es[2]) + rng.randint(1, 10 ** 6)
            else:
                mode = rng.randint(0, 3)
                if mode == 0: x = Fr(rng.randint(-p ** 8, p ** 8))
                elif mode == 1: x = Fr(rng.choice(self.e)) + Fr(p) ** rng.randint(0, 40) * rng.randint(1, p ** 6)
                elif mode == 2: x = Fr(rng.randint(1, p ** 6) * rng.choice([1, -1]), p ** (2 * rng.randint(1, 20)))
                else: x = Fr(rng.choice(self.e)) + Fr(rng.randint(1, p ** 6) * rng.choice([1, -1]), p ** rng.randint(1, 30))
            fx = self.f(x)
            if fx == 0 or not is_sq_local(fx, p): continue
            vecs.append(self.locvec(self.kummer(x), p)); B = span(vecs)
        if len(B) != tdim: raise RuntimeError('local image at %d: dim %d, expected %d' % (p, len(B), tdim))
        return [b for _, b in B]
    def selmer(self, seed=1):
        rng = random.Random(seed)
        places = [0] + self.primes
        n = len(self.S)
        gens = [(Fr(self.S[i]), Fr(1)) for i in range(n)] + [(Fr(1), Fr(self.S[i])) for i in range(n)]
        quot = []
        for p in places:
            Wb = span(self.local_image(p, rng))
            quot.append((p, Wb, 2 * (3 if p == 2 else (1 if p == 0 else 2))))
        def image(g):
            out = 0; shift = 0
            for p, Wb, width in quot:
                v = self.locvec(g, p)
                for piv, b in Wb:
                    if v >> piv & 1: v ^= b
                out |= v << shift; shift += width
            return out
        basis = []; kernel = []
        for i, g in enumerate(gens):
            img, comb = image(g), 1 << i
            for piv, (bi, bc) in basis:
                if img >> piv & 1: img ^= bi; comb ^= bc
            if img: basis.append((img.bit_length() - 1, (img, comb)))
            else: kernel.append(comb)
        self.sel = []
        for _, c in span(kernel):
            d1 = Fr(1); d2 = Fr(1)
            for i in range(len(gens)):
                if c >> i & 1: d1 *= gens[i][0]; d2 *= gens[i][1]
            self.sel.append((d1, d2))
        return len(self.sel)
    def classvec(self, pair):
        """bit vector of the class in Q(S,2)^2 of a pair of rationals (assumed in Q(S,2)^2)"""
        out = 0
        for k, q in enumerate(pair):
            q = Fr(q)
            for i, s in enumerate(self.S):
                b = (1 if q < 0 else 0) if s == -1 else (vp(q.numerator, s)[0] - vp(q.denominator, s)[0]) % 2
                out |= b << (k * len(self.S) + i)
        return out
    def sel_contains(self, pair):
        B = span([self.classvec(s) for s in self.sel]); v = self.classvec(pair)
        for piv, b in B:
            if v >> piv & 1: v ^= b
        return v == 0

if __name__ == '__main__':
    print('# my2desc.py self-test (review z7; written by an AI model): dim Sel^2 of y^2 = x^3 - n^2 x by the complete 2-descent')
    print('# reference: PARI 2.17.4 ellrank [r,R,s] recorded at review time (not computed here); R + s + 2 must equal dim Sel^2')
    REF = {1: (0, 0, 0), 5: (1, 1, 0), 6: (1, 1, 0), 7: (1, 1, 0), 13: (1, 1, 0), 14: (1, 1, 0), 15: (1, 1, 0), 17: (0, 0, 2),
           21: (1, 1, 0), 22: (1, 1, 0), 23: (1, 1, 0), 29: (1, 1, 0), 30: (1, 1, 0), 31: (1, 1, 0), 34: (2, 2, 0), 41: (2, 2, 0),
           65: (2, 2, 0), 73: (0, 0, 2), 82: (0, 0, 2), 89: (0, 0, 2), 97: (0, 0, 2), 105: (0, 0, 0), 113: (0, 2, 0),
           137: (2, 2, 0), 145: (2, 2, 0), 193: (0, 0, 2), 205: (1, 1, 2), 210: (2, 2, 0), 219: (2, 2, 0), 226: (2, 2, 0),
           241: (0, 0, 2)}
    for n, (r, R, s) in REF.items():
        d = Curve(-n, 0, n).selmer()
        assert d == R + s + 2, n
        print('n = %d: dim Sel^2 = %d; PARI [%d, %d, %d]: R + s + 2 = %d, agrees' % (n, d, r, R, s, R + s + 2))
    print('# all agree')
