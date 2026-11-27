#!/usr/bin/env python3
"""Review z7 of roadmap H, H7 claim 3: the partner map carries the nine cells, not only the five forms.
Written by an AI model (review agent z7, 27 November 2026); not read by a person.  Experiment, not counted.

With M = m + n, N = m - n (and the halved pair M/2, N/2) and t' = (t + 1)/(t - 1), the nine cells of global.py's
classify_t for (M, N) at t' are the nine cells of (m, n) at t, permuted, times one common factor (4/(t-1)^2, resp.
1/(t-1)^2).  So SOLUTION, degenerate and no are preserved (G1Relations.lean's partner_forms / partner_forms_half give
the five forms in Lean; the four cells (mt +- n)^2, (nt +- m)^2 are checked here).  And (6881,6233) -> (6557,324).
sympy only; no PARI.

    python3 partner_cells.py         prints partner_cells.txt   (a few seconds)
"""
import sympy as sp
from math import gcd
print('# partner_cells.py (review z7; written by an AI model)')
m, n, t = sp.symbols('m n t')
def cells(m, n, t):
    s = (m * m + n * n) / sp.Integer(2)
    F0 = m * m + n * n * t * t; F8 = n * n + m * m * t * t; F4 = s * (1 + t * t); L = F4 - 2 * m * n * t; U = F4 + 2 * m * n * t
    return [(m * t + n) ** 2, (m * t - n) ** 2, (n * t + m) ** 2, (n * t - m) ** 2, F0, F8, F4, L, U]
c = cells(m, n, t); tp = (t + 1) / (t - 1)
for nm, (M, N, sc) in (('(m+n, m-n)', (m + n, m - n, 4)), ('((m+n)/2, (m-n)/2)', ((m + n) / 2, (m - n) / 2, 1))):
    cp = [sp.simplify(x * (t - 1) ** 2 / sc) for x in cells(M, N, tp)]
    perm = [next((j for j, y in enumerate(c) if sp.simplify(x - y) == 0), None) for x in cp]
    assert sorted(perm) == list(range(9))
    print('partner %s at t\' = (t+1)/(t-1): cell i times (t-1)^2/%d = cell %s of (m, n) at t' % (nm, sc, perm))
g = gcd(6881 + 6233, 6881 - 6233)
print('(6881,6233): m, n odd, g = %d, halved partner ((m+n)/2, (m-n)/2) = (%d,%d)' % (g, (6881 + 6233) // g, (6881 - 6233) // g))
print('# end')
