#!/usr/bin/env python3
"""Guide only (PARI's ellrank, floating point inside): the exact Selmer bounds of selmer_square_5000.json against
PARI's rank bounds [lo, hi] on 120 random curves with m <= 1500, R supplied.  Nothing here is used in a proof.

    python3 selmer_square_pari.py [--write]      (writes selmer_square_pari.txt with --write)

Expected: hi <= s on every curve (the 2-Selmer rank is at most the 2-isogeny bound: dim Sel^2(E) <= dim Sel(E) +
dim Sel(E') - 1, since T' is not in phi(E(Q)[2])); and lo = hi = s mod 2 (the parity remark of SELMER_SQUARE.md).
Written by an AI model (Claude); not read by a person.
"""
import os, sys, json, random
from collections import Counter
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import sigma_exact as SE
import cypari2
WRITE = '--write' in sys.argv
pari = cypari2.Pari(); pari.allocatemem(4 * 10**8)
out = []
def say(s): print(s, flush=True); out.append(s)
d = json.load(open(os.path.join(here, 'selmer_square_5000.json')))
rnd = random.Random(5); sample = [r for r in d if r['m'] <= 1500]; rnd.shuffle(sample); sample = sample[:120]
bad = par = 0; c = Counter()
for r in sample:
    m, nn = r['m'], r['n']; M = SE.model(m, nn, 2)
    E = pari.ellinit([0, M['a2'], 0, M['a4'], M['a6']])
    S = m*m + nn*nn; xR = -4*m*nn*S*(m*m - nn*nn); R = [xR, pari.ellordinate(E, xR)[0]]
    rr = pari.ellrank(E, 0, [R]); lo, hi = int(rr[0]), int(rr[1])
    c[(r['bound'], lo, hi)] += 1
    if hi > r['bound'] or lo > hi: bad += 1; say('INCONSISTENT %d %d s=%d PARI [%d,%d]' % (m, nn, r['bound'], lo, hi))
    if lo % 2 == r['bound'] % 2 and hi % 2 == r['bound'] % 2: par += 1
say('selmer_square_pari.py (guide only): %d curves, PARI upper bound above s on %d; lo, hi of the parity of s on %d'
    % (len(sample), bad, par))
for k in sorted(c): say('  s = %d, PARI [%d, %d]: %d curves' % (k[0], k[1], k[2], c[k]))
assert bad == 0
if WRITE: open(os.path.join(here, 'selmer_square_pari.txt'), 'w').write('\n'.join(out) + '\n')
