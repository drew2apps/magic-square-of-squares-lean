#!/usr/bin/env python3
"""For every coprime n < m <= N with m^2 + n^2 a square or twice a square, build the degenerate magic square
of squares of BARRIER.md and check it exactly.

Square:  a b c / d e f / g h i, entries A = a^2, ..., I = i^2.  Theorem 1 of PAPER.md:
  b + d = m u, h + f = n u, b - d = n v, f - h = m v, 8 e^2 = (m^2 + n^2)(u^2 + v^2).
Type 1 (m^2 + n^2 = k^2): u = v = 2, roots b = m + n, d = m - n, f = m + n, h = n - m, e = k; the square is
  k^2 (m+n)^2 (m-n)^2 / (m-n)^2 k^2 (m+n)^2 / (m+n)^2 (m-n)^2 k^2.
Type 2 (m^2 + n^2 = 2k^2): u = 0, v = 2, roots b = n, d = -n, f = m, h = -m, e = k; the square is
  m^2 n^2 k^2 / n^2 k^2 m^2 / k^2 m^2 n^2.
Checked: the eight line sums, that the nine entries are squares of the stated roots (the edge entries and the
centre) or squares (the corners), the four equations of Theorem 1 and the centre equation, and that the entries
are not all equal.  For type 2, b + d = h + f = 0, so the ratio is carried by b - d = 2n, f - h = 2m."""
import sys
from math import gcd, isqrt


def sq(x): return x >= 0 and isqrt(x)**2 == x


def magic(M):
    s = sum(M[0])
    L = [r for r in M] + [list(c) for c in zip(*M)] + [[M[i][i] for i in range(3)], [M[i][2 - i] for i in range(3)]]
    return all(sum(l) == s for l in L)


N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
cnt = {1: 0, 2: 0}; fail = 0
for m in range(2, N + 1):
    for n in range(1, m):
        if gcd(m, n) != 1: continue
        S = m*m + n*n
        ty = 1 if sq(S) else 2 if (S % 2 == 0 and sq(S // 2)) else 0
        if not ty: continue
        k = isqrt(S) if ty == 1 else isqrt(S // 2)
        if ty == 1:
            u, v = 2, 2; b, d, f, h = m + n, m - n, m + n, n - m
            M = [[k*k, b*b, (m - n)**2], [d*d, k*k, f*f], [(m + n)**2, h*h, k*k]]
        else:
            u, v = 0, 2; b, d, f, h = n, -n, m, -m
            M = [[m*m, b*b, k*k], [d*d, k*k, f*f], [k*k, h*h, n*n]]
        ok = (magic(M) and all(sq(x) for r in M for x in r)
              and M[0][1] == b*b and M[1][0] == d*d and M[1][2] == f*f and M[2][1] == h*h
              and b + d == m*u and h + f == n*u and b - d == n*v and f - h == m*v
              and 8*k*k == S*(u*u + v*v) and M[1][1] == k*k
              and len({x for r in M for x in r}) == 3)
        if ok: cnt[ty] += 1
        else: fail += 1; print('FAIL', m, n, M)
print('N = %d: degenerate magic squares of squares checked: %d with m^2+n^2 a square, %d with twice a square; failures %d'
      % (N, cnt[1], cnt[2], fail))
