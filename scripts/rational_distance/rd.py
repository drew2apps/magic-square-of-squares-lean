#!/usr/bin/env python3
"""The rational distance problem for the unit square, line by line (EXPERIMENT; RATIONAL_DISTANCE.md).

A point (x, y) with rational coordinates has rational distances to the four corners 0, 1, i, 1 + i exactly when the
four numbers z, z - 1, z - i, z - 1 - i (z = x + i y) have square norms.  On the line y = r/s, with x = u/v in lowest
terms, these are, up to the common factor s v, the linear forms
    s u + i r v,   s (u - v) + i r v,   s u + i (r - s) v,   s (u - v) + i (r - s) v
over Z[i].  A form of square norm is a squarefree positive integer (its class) times a square of Q(i); the class is
supported on the primes of 2 r s (first two forms) and of 2 s (r - s) (last two).  First step: the tuples of four classes
that pass every local test (tuple_descent.py of ../cover_descent, the program used for the magic square).

Second step, by hand (RATIONAL_DISTANCE.md): if the four classes are equal, z - a_j = D w_j^2 with one D, and
(w_1 : w_2 : w_3 : w_4) is a point over Q(i) of the curve
    C:  w_3^2 = (1 - i) w_1^2 + i w_2^2,   w_4^2 = -i w_1^2 + (1 + i) w_2^2,
of genus 1, whose Jacobian is y^2 = x^3 - x.  That curve has rank 0 over Q(i), so C has exactly the eight points
(1 : +-1 : +-1 : +-1), where w_1^2 = w_2^2, which is impossible.  So a line on which the first step leaves only tuples of
four equal classes has no point with four rational distances.

usage: rd.py line r,s                 the tuples left on the line y = r/s
       rd.py scan S0 S1 [K]           every line y = r/s in lowest terms with S0 <= s <= S1 and 0 < r/s < 1/2
                                      (with K: 1 < r/s <= K instead); one line of output per line, and a summary
       rd.py curve                    points of C modulo primes p = 1 mod 4, against y^2 = x^3 - x
       rd.py control [C]              points with THREE rational distances (search with denominator below C, default 400):
                                      the first step with three forms has to keep the tuple of each of them"""
import sys, os
from math import gcd, isqrt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'cover_descent'))
import tuple_descent as td
from sympy import factorint

def square(r, s):
    forms = [((s, 0), (0, r)), ((s, 0), (-s, r)), ((s, 0), (0, r - s)), ((s, 0), (-s, r - s))]
    pr = lambda x: sorted(set([2]) | set(factorint(abs(x))))
    sup = [pr(r*s), pr(r*s), pr(s*(r - s)), pr(s*(r - s))]
    known = sorted(set([tuple(td.cls_of((al[0]*u + be[0]*v, al[1]*u + be[1]*v)) for al, be in forms) for u, v in ((1, 0),)]))
    return forms, sup, known

def line(r, s):
    forms, sup, known = square(r, s)
    T, used, sizes = td.descend(forms, sup, known, 150)
    mixed = [t for t in T if len(set(t)) > 1]
    return T, mixed

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'line':
        r, s = (int(x) for x in sys.argv[2].split(',')); T, mixed = line(r, s)
        print('y = %d/%d: tuples left: %s' % (r, s, T)); print('all four classes equal in every tuple: %s' % (not mixed))
    elif mode == 'scan':
        S0, S1 = int(sys.argv[2]), int(sys.argv[3]); K = int(sys.argv[4]) if len(sys.argv) > 4 else None
        n = 0; bad = []; hist = {}
        for s in range(S0, S1 + 1):
            for r in (range(1, s) if K is None else range(s + 1, K*s + 1)):
                if gcd(r, s) != 1 or (K is None and 2*r >= s): continue
                T, mixed = line(r, s); n += 1; hist[len(T)] = hist.get(len(T), 0) + 1
                if mixed: bad.append('%d/%d' % (r, s))
                print('%d/%d %d %s' % (r, s, len(T), 'equal' if not mixed else 'MIXED ' + str(mixed[:4])), flush=True)
        print('lines y = r/s with %d <= s <= %d, %s: %d; tuples left (number: lines): %s; lines with a tuple whose four classes are not all equal: %d %s'
              % (S0, S1, '0 < r/s < 1/2' if K is None else '1 < r/s <= %d' % K, n, sorted(hist.items()), len(bad), bad))
    elif mode == 'curve':
        for p in (5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97):
            for iota in [x for x in range(p) if (x*x + 1) % p == 0]:
                sq = {}
                for w in range(p): sq.setdefault(w*w % p, []).append(w)
                n = 0
                for w1, w2 in [(1, t) for t in range(p)] + [(0, 1)]:
                    n += len(sq.get(((1 - iota)*w1*w1 + iota*w2*w2) % p, []))*len(sq.get((-iota*w1*w1 + (1 + iota)*w2*w2) % p, []))
                e = 1 + sum(1 + (1 if pow((x**3 - x) % p, (p - 1)//2, p) == 1 else (0 if (x**3 - x) % p == 0 else -1)) for x in range(p))
                print('p = %d, i -> %d: points of C: %d; points of y^2 = x^3 - x: %d; equal: %s' % (p, iota, n, e, n == e))
    elif mode == 'control':
        C = int(sys.argv[2]) if len(sys.argv) > 2 else 400
        issq = lambda x: x >= 0 and isqrt(x)**2 == x
        pts = [(a, b, c) for c in range(2, C) for b in range(1, 2*c) if b != c for a in range(1, 2*c)
               if a != c and gcd(gcd(a, b), c) == 1 and issq(a*a + b*b) and issq((a - c)**2 + b*b) and issq(a*a + (b - c)**2)]
        kept = eq = 0
        for a, b, c in pts:
            g = gcd(b, c); r, s = b//g, c//g; g2 = gcd(a, c); u, v = a//g2, c//g2
            forms, sup, known = square(r, s)
            tup = tuple(td.cls_of((al[0]*u + be[0]*v, al[1]*u + be[1]*v)) for al, be in forms[:3])
            T, used, sizes = td.descend(forms[:3], sup[:3], [k[:3] for k in known], 150)
            kept += tup in T; eq += len(set(tup)) == 1
        print('points (a/c, b/c), c < %d, with rational distances to (0,0), (1,0), (0,1): %d, the first: %s' % (C, len(pts), pts[:4]))
        print('the tuple of the point is kept by the first step on its line (three forms): %d of %d; points whose three classes are equal: %d' % (kept, len(pts), eq))
