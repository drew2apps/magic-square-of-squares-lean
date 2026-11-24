#!/usr/bin/env python3
"""Perfect triangle (Heron triangle with three rational medians), one line y = r/s at a time (EXPERIMENT; PERFECT_TRIANGLE.md).
A = 0, B = 1, C = z = x + iy.  Sides |z|, |z - 1|, 1; medians |z + 1|/2, |z - 2|/2, |2z - 1|/2; area y/2.
Forms on the line y = r/s, x = u/v (times s v):  z, z - 1, z + 1, z - 2, 2z - 1.  First step: tuple_descent.py of
../cover_descent.  Written by a reviewer instance of the AI model on 24 November 2026; kept with small changes.
usage: pt.py line r,s | pt.py scan S0 S1 [R]     (lines with S0 <= s <= S1 and 0 < r/s <= R, default 2)"""
import sys, os, time
from math import gcd, isqrt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'cover_descent'))
import tuple_descent as td
from sympy import factorint

def triangle(r, s, which=(0, 1, 2, 3, 4)):
    allf = [((s, 0), (0, r)), ((s, 0), (-s, r)), ((s, 0), (s, r)), ((s, 0), (-2*s, r)), ((2*s, 0), (-s, 2*r))]
    forms = [allf[j] for j in which]
    P = sorted(set([2]) | set(factorint(r*s)))
    sup = [P for _ in forms]
    known = sorted(set([tuple(td.cls_of((al[0]*u + be[0]*v, al[1]*u + be[1]*v)) for al, be in forms) for u, v in ((1, 0),)]))
    return forms, sup, known

def line(r, s, which=(0, 1, 2, 3, 4), pmax=150):
    forms, sup, known = triangle(r, s, which)
    T, used, sizes = td.descend(forms, sup, known, pmax)
    return T, known

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'line':
        r, s = (int(x) for x in sys.argv[2].split(',')); T, known = line(r, s)
        print('y = %d/%d: known (infinity): %s; tuples left (%d): %s' % (r, s, known, len(T), T))
    elif mode == 'scan':
        S0, S1 = int(sys.argv[2]), int(sys.argv[3]); R = float(sys.argv[4]) if len(sys.argv) > 4 else 2.0
        n = 0; hist = {}; bad = []; t0 = time.time()
        for s in range(S0, S1 + 1):
            for r in range(1, int(R*s) + 1):
                if gcd(r, s) != 1: continue
                T, known = line(r, s); n += 1
                g0 = known[0]
                # degenerate shape: g*(1,1,1,1,2) up to squares, g squarefree
                deg = all(len(set(t[:4])) == 1 and td.sqf(2*t[0]*t[4]) == 1 for t in T)
                hist[len(T)] = hist.get(len(T), 0) + 1
                if not deg: bad.append('%d/%d' % (r, s))
                print('%d/%d %d %s' % (r, s, len(T), 'degenerate-only' if deg else 'OTHER ' + str([t for t in T if not (len(set(t[:4])) == 1 and td.sqf(2*t[0]*t[4]) == 1)][:4])), flush=True)
        print('lines y = r/s, %d <= s <= %d, 0 < r/s <= %g: %d lines; tuples left (number: lines): %s; lines with a tuple not of the shape g*(1,1,1,1,2): %d %s'
              % (S0, S1, R, n, sorted(hist.items()), len(bad), bad))
