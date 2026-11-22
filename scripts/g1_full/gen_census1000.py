#!/usr/bin/env python3
"""Writes MagicSquareOfSquares/G1Census1000.lean: the census of the full G1 system for m <= 1000.

Each pair is closed by a ratio class modulo one of the 34 primes (looked up in a bit mask), by one
of the two 2-adic classes, or by a prime 5 mod 8 dividing m^2 + n^2 exactly once. The pairs closed
by none of these are listed as exceptions. The kernel evaluates the census in blocks."""
import math, os
here = os.path.dirname(os.path.abspath(__file__))
N = 1000
table = {}
for line in open(os.path.join(here, 'modp_classes.txt')):
    p, r = map(int, line.split()); table.setdefault(p, set()).add(r)
primes = sorted(table)
sqmask = lambda p: sum(1 << s for s in {x * x % p for x in range(p)})
rmask = lambda p: sum(1 << r for r in table[p])
five = [5, 13, 29, 37, 53, 61, 109, 149, 157, 173, 181, 269, 373, 541, 557, 701, 773, 3461, 12277]
def cond8(a, b): return (a % 4 == 2 and b % 2 == 1) or (b % 4 == 2 and a % 2 == 1)
def cond32(a, b): return ((a % 8 in (1, 7)) and (b % 8 in (3, 5))) or ((a % 8 in (3, 5)) and (b % 8 in (1, 7)))
def closed(m, n):
    for p in primes:
        if n % p and (m * pow(n, p - 2, p)) % p in table[p]: return True
    if cond8(m, n) or cond32(m, n): return True
    s = m * m + n * n
    return any(s % p == 0 and s % (p * p) != 0 and n % p != 0 for p in five)
exc = [(m, n) for m in range(2, N + 1) for n in range(1, m) if math.gcd(m, n) == 1 and not closed(m, n)]
issq = lambda x: math.isqrt(x) ** 2 == x
triv = sum(1 for m, n in exc if issq(m * m + n * n) or ((m * m + n * n) % 2 == 0 and issq((m * m + n * n) // 2)))
npairs = sum(1 for m in range(2, N + 1) for n in range(1, m) if math.gcd(m, n) == 1)
# blocks of about 25,000 pairs (m, n) with n < m
bounds = [0]; acc = 0
for m in range(0, N + 1):
    acc += m
    if acc >= 25000: bounds.append(m + 1); acc = 0
if bounds[-1] != N + 1: bounds.append(N + 1)
blocks = list(zip(bounds[:-1], bounds[1:]))
rows = ',\n'.join(f'  ({p}, {sqmask(p)}, {rmask(p)})' for p in primes)
block_thms = '\n\n'.join(f'theorem block_{k} : blockOk {lo} {hi - lo} = true := by decide +kernel' for k, (lo, hi) in enumerate(blocks))
disj = ' ∨ '.join([f'm < {blocks[0][1]}'] + [f'({lo} ≤ m ∧ m < {hi})' for lo, hi in blocks[1:]])
cases = '\n'.join(f'  · exact pairOk_of_block block_{k} (by omega) (by omega) hn' for k in range(len(blocks)))
T = open(os.path.join(here, 'census1000.template')).read()
out = (T.replace('@ROWS@', rows).replace('@FIVE@', str(five)).replace('@FIVEALTS@', ' | '.join(['rfl'] * len(five)))
        .replace('@EXC@', str(exc)).replace('@NEXC@', str(len(exc))).replace('@TRIV@', str(triv))
        .replace('@NPAIRS@', f'{npairs:,}').replace('@NCLOSED@', f'{npairs - len(exc):,}')
        .replace('@BLOCKS@', block_thms).replace('@DISJ@', disj).replace('@CASES@', cases)
        .replace('@RCASES@', ' | '.join(['h'] * len(blocks))))
open(os.path.join(here, '..', '..', 'MagicSquareOfSquares', 'G1Census1000.lean'), 'w', encoding='utf-8').write(out)
print('pairs', npairs, 'exceptions', len(exc), 'with a trivial point', triv, 'blocks', len(blocks))
