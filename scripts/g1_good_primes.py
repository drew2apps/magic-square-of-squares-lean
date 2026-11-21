#!/usr/bin/env python3
"""G1 family: which ratios (m:n) have no point modulo a prime of good reduction?

The three conditions
    a^2 = m^2 V^2 + n^2 U^2,   b^2 = 2 (m^2 + n^2)(U^2 + V^2),   c^2 = n^2 V^2 + m^2 U^2
define a curve C of genus 5 (a cover of the (U:V)-line branched at six points) whenever
p does not divide 2 m n (m^2 - n^2)(m^2 + n^2). This script counts the points of C over F_p
for every such ratio and every odd prime up to the limit, reports the ratios with no point,
and checks the Weil bound |N - (p + 1)| <= 10 sqrt(p).

usage: g1_good_primes.py [limit]      (default 400)
"""
import math
import sys


def primes(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


LIM = int(sys.argv[1]) if len(sys.argv) > 1 else 400
violations, worst, obstructed = 0, 0.0, {}
for p in primes(LIM):
    if p == 2:
        continue
    chi = [0] * p
    for x in range(1, p):
        chi[x * x % p] = 1
    for x in range(1, p):
        if chi[x] == 0:
            chi[x] = -1
    pts = [(1, v) for v in range(p)] + [(0, 1)]
    n = 1
    for m in range(1, p):                    # the ratios (1:0) and (0:1) have bad reduction
        if (m * m - 1) % p == 0 or (m * m + 1) % p == 0:
            continue                         # bad reduction
        N = 0
        for (u, v) in pts:
            q1 = (m * m * v * v + n * n * u * u) % p
            q2 = (2 * (m * m + n * n) * (u * u + v * v)) % p
            q3 = (n * n * v * v + m * m * u * u) % p
            N += (1 + chi[q1]) * (1 + chi[q2]) * (1 + chi[q3])
        dev = abs(N - (p + 1)) / (10 * math.sqrt(p))
        worst = max(worst, dev)
        if dev > 1:
            violations += 1
        if N == 0:
            obstructed.setdefault(p, []).append(m)
print('odd primes up to', LIM)
print('Weil bound violations:', violations, '| worst |N-(p+1)| / (10 sqrt p) =', round(worst, 3))
for p, ms in obstructed.items():
    print('p =', p, ': no point for m/n =', ms)
print('smallest prime with p + 1 - 10 sqrt(p) > 0:',
      next(p for p in primes(200) if p + 1 - 10 * math.sqrt(p) > 0))
