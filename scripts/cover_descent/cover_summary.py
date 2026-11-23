#!/usr/bin/env python3
"""Counts from the outputs of cover_batch.py.  usage: cover_summary.py cover_34.txt cover_123.txt"""
import sys, re
for fn in sys.argv[1:]:
    L = [l for l in open(fn) if l.startswith('(')]
    none = [l for l in L if 'no class is possible' in l]; ex = [l for l in L if 'exactly the classes of the known points' in l]; more = [l for l in L if 'MORE' in l]
    sq = [l for l in ex if "'t = 1, -1'" in l]; tw = [l for l in ex if "'t = 0, infinity'" in l]
    print('%s: %d ratios; no class possible: %d; exactly the class of the known points: %d (%d at t = 1, -1 and %d at t = 0, infinity); one or more further classes: %d' % (fn, len(L), len(none), len(ex), len(sq), len(tw), len(more)))
    if more: print('  further classes: ' + '; '.join('%s %s' % (re.match(r'\(\d+,\d+\)', l).group(0), re.search(r'left: (\[[^\]]*\])', l).group(1)) for l in more))
    need = [re.match(r'\(\d+,\d+\)', l).group(0) for l in L if 'good primes up to 150;' not in l]
    print('  ratios that needed good primes above 150: %d: %s' % (len(need), ' '.join(need)))
