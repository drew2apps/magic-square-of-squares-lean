"""Reviewer z16 (an AI model; not read by a person): summary of z16_five.py on the 937 ratios (out937.jsonl),
against EXPLORATION_5 §2's numbers, the session's own f4_all.jsonl (data/), the pairing flags of EXPLORATION_4
(data/out_*.jsonl), and the tier-A condition of G1Primes.lean / IsMSS.edge_norm: every prime = 5 mod 8 divides
m^2 + n^2 to an even power.  usage: python3 -I z16_summary.py"""
import os, json, collections
from sympy import factorint
HERE = os.path.dirname(os.path.abspath(__file__))
def load(f): return [json.loads(l) for l in open(os.path.join(HERE, f)) if l.startswith('{')]
mine = {(r['m'], r['n']): r for r in load('out937.jsonl')}
pair = {(r['m'], r['n']): r for f in ('out_997.jsonl', 'out_15000_0.jsonl', 'out_15000_1.jsonl') for r in load('data/' + f)}
sess = {}
for r in load('data/f4_all.jsonl'):
    sess.setdefault((r['m'], r['n']), r)                          # the first record (the list run) for each ratio
print('ratios: mine %d, EXPLORATION_4 tables %d, same set: %s' % (len(mine), len(pair), set(mine) == set(pair)))
n4 = sum(r['n4'] for r in mine.values()); n4p = sum(len(r['tuples']) for r in pair.values())
same4 = all(mine[k]['n4'] == len(pair[k]['tuples']) for k in mine)
print('four-form tuples: mine %d, EXPLORATION_4 %d, equal at every ratio: %s' % (n4, n4p, same4))
keep = [k for k, r in mine.items() if r['n5']]
ext = sum(len(r['ext']) for r in mine.values()); n5 = sum(r['n5'] for r in mine.values())
print('five forms: ratios that keep a 5-tuple %d of %d; 4-tuples that extend %d of %d; 5-tuples %d' % (len(keep), len(mine), ext, n4, n5))
agree = all(sorted(map(tuple, mine[k]['ext'])) == sorted(map(tuple, sess[k]['ext'])) and mine[k]['n5'] == sess[k]['n5'] for k in mine if k in sess)
print('against the session\'s f4_all.jsonl (%d of the 937 present): same n5 and same extending 4-tuples at every ratio: %s' % (sum(k in sess for k in mine), agree))
# pairing survivors (no curve kills; 'outside' counted as not killed, as f4_all.py does; how many 'outside' there are)
outside = sum(1 for r in pair.values() for t, k in r['tuples'] for v in k.values() if v == 'outside')
surv = {k: [tuple(t) for t, kk in r['tuples'] if not any(v == 1 for v in kk.values())] for k, r in pair.items()}
sr = [k for k, s in surv.items() if s]
print('pairing flags: %d entries "outside Sel^2"; pairing survivors %d tuples at %d ratios: %s' % (outside, sum(map(len, surv.values())), len(sr), sorted(sr)))
se = {k: [t for t in surv[k] if t in set(map(tuple, mine[k]['ext']))] for k in sr}
print('pairing survivors that extend to a 5-tuple: %s' % {k: v for k, v in se.items() if v})
for k in sorted(sr):
    print('   %s: survivors %s; 5-tuples here %d' % (k, surv[k], mine[k]['n5']))
# the condition at primes = 5 mod 8
def bad5(m, n): return sorted(q for q, e in factorint(m*m + n*n).items() if q % 8 == 5 and e % 2)
b5 = [k for k in mine if bad5(*k)]
print('ratios of the 937 with a prime = 5 mod 8 to an odd power in m^2 + n^2 (excluded by G1Primes.lean, tier A): %d' % len(b5))
print('   of these, keep a 5-tuple: %d;   ratios without such a prime: %d, of which keep a 5-tuple: %d'
      % (sum(1 for k in b5 if mine[k]['n5']), len(mine) - len(b5), sum(1 for k in mine if k not in b5 and mine[k]['n5'])))
t4 = sum(mine[k]['n4'] for k in mine if k not in b5); e4 = sum(len(mine[k]['ext']) for k in mine if k not in b5)
print('   four-form tuples at the ratios without such a prime: %d, of which extend: %d' % (t4, e4))
print('   pairing-survivor ratios with such a prime: %s' % sorted(k for k in sr if k in b5))
# the support of d4: d4 in the 5-tuples left is 2^a times the primes of odd exponent in m^2+n^2 (a check of the support)
def oddpart(m, n): return {q for q, e in factorint(m*m + n*n).items() if e % 2 and q > 2}
viol = [(k, t) for k in keep for t in mine[k]['t5'] if not {q for q in factorint(t[4])} <= ({2} | oddpart(*k))]
print('5-tuples whose d4 has a prime of even exponent in m^2 + n^2: %d' % len(viol))
# the three samples beyond m = 15000 (EXPLORATION_5 §2 table: 46 of 269 biased ratios keep a 5-tuple; 0 survivors extend)
samp = {(r['m'], r['n']): r for r in load('out_samples.jsonl')}
for s in ('biased', 'random', 'fewp'):
    rows = [r for r in load('data/samp_%s.jsonl' % s) if 'tuples' in r]
    k = [(r['m'], r['n']) for r in rows]
    same4 = all(samp[x]['n4'] == len(r['tuples']) for x, r in zip(k, rows))
    keep = sum(1 for x in k if samp[x]['n5']); ext = sum(len(samp[x]['ext']) for x in k)
    sv = [(x, tuple(t)) for x, r in zip(k, rows) for t, kk in r['tuples'] if not any(v == 1 for v in kk.values())]
    print('sample %s: %d ratios, four-form tuples equal to the session\'s at every ratio: %s; keep a 5-tuple %d; 4-tuples that extend %d of %d; '
          'pairing survivors %d, of which extend %d' % (s, len(k), same4, keep, ext, sum(len(r['tuples']) for r in rows), len(sv),
                                                         sum(list(t) in samp[x]['ext'] for x, t in sv)))
# the nine open census ratios (and the partner of the last ghost): the five-form tuples left are the degenerate ones
import math, sys
sys.path.insert(0, HERE)
import z16_witness as W
for r in load('out_open.jsonl'):
    m, n = r['m'], r['n']; C = m*m + n*n; F = W.forms5(m, n)
    pts = ([(1, 1), (-1, 1)] if math.isqrt(C)**2 == C else []) + ([(0, 1), (1, 0)] if math.isqrt(2*C)**2 == 2*C else [])
    deg = sorted(tuple(W.my_class(W.val(f, u, v)) for f in F) for u, v in pts)
    if not pts:
        print('(%d,%d), not degenerate: four-form %d, five-form %d' % (m, n, r['n4'], r['n5'])); continue
    print('open (%d,%d): four-form %d, five-form %d; the five-form tuples are exactly the %d degenerate ones: %s'
          % (m, n, r['n4'], r['n5'], len(deg), sorted(map(tuple, r['t5'])) == deg))
