"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
Claim 1 of H3 says good primes p >= 1154 remove nothing.  This scans every good residue r = m/n mod p at every prime
p = 1 mod 4 with LO <= p <= HI with my_local.allowed and prints, per prime: p, good residues, residues that lose a
vector, least size of an allowed set; then the primes with a loss.  Rests on my_local.py.
usage: scan_losses.py [LO HI] [--write]   (default 1153 5000 -> scan_losses.txt; about 12 minutes on one core)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit, isprime, good, ARGS
import my_local as M
lo, hi = (int(ARGS[0]), int(ARGS[1])) if len(ARGS) >= 2 else (1153, 5000)
out = ['z8: losses of the four-class set at good primes p = 1 mod 4, %d <= p <= %d (every good residue). Written by an AI model; not read by a person.' % (lo, hi),
       'columns: p, good residues, residues losing a vector, least allowed size']
withloss = []
for p in range(lo, hi + 1):
    if p % 4 != 1 or not isprime(p): continue
    ch = M.chitab(p); lost = 0; tot = 0; mn = 16
    for r in range(1, p):
        if not good(r, 1, p): continue
        a = M.allowed(r, 1, p, ch); tot += 1
        if a != (1 << 16) - 1: lost += 1
        mn = min(mn, bin(a).count('1'))
    out.append('%d %d %d %d' % (p, tot, lost, mn))
    if lost: withloss.append(p)
out.append('primes with a loss: %d; the largest: %s' % (len(withloss), max(withloss) if withloss else '-'))
emit(out, 'scan_losses.txt' if (lo, hi) == (1153, 5000) else 'scan_losses_%d_%d.txt' % (lo, hi))
