"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
octic_cut.py's argument at every octic prime of badprimes_left.py (24 ratios), against exact_local's verdict there.
Rests on octic_cut.py and on badprimes_left.txt (committed; rerun badprimes_left.py first if needed).
usage: octic_cut_all.py [--write]   (octic_cut_all.txt; a few s)"""
import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit, load_phantoms, HERE
import octic_cut as O
R = load_phantoms()
out = ['z8: octic_cut.py at every octic prime above 20,000 of the 24 ratios, against exact_local (badprimes_left.txt). Written by an AI model; not read by a person.']
cur = None; agree = 0; tot = 0
for l in open(os.path.join(HERE, 'badprimes_left.txt')):
    mm = re.match(r'\((\d+),(\d+)\): \d+ such primes', l.strip())
    if mm: cur = (int(mm.group(1)), int(mm.group(2))); continue
    mm = re.match(r'(\d+) \((m\^2\+n\^2|octic), p mod 4 = \d\): \d+ allowed label tuples, (.*)', l.strip())
    if not mm or mm.group(2) != 'octic': continue
    p = int(mm.group(1)); el_cut = mm.group(3).startswith('cuts')
    known, extra = R[cur]
    r, valid, kok, cut = O.check(cur[0], cur[1], p, known, extra)
    tot += 1; agree += (len(cut) > 0) == el_cut
    out.append('%s p = %d: %d meetings, %s; known %s; argument cuts %d of %d extra; exact_local %s' % (cur, p, len(r), 'all with v = 1' if valid else 'NOT ALL v = 1', 'ok' if kok else 'VIOLATE', len(cut), len(extra), 'cuts' if el_cut else 'keeps'))
out.append('octic primes: %d; the argument cuts exactly where exact_local cuts: %d of %d' % (tot, agree, tot))
emit(out, 'octic_cut_all.txt')
