#!/usr/bin/env python3
"""Checks for dm_sigma.py (G1_FULL_SYSTEM.md, "The curves of F0*F8*GL").  Python only, no Sage.

    python3 check_sigma.py                 the committed records: counts, no alert, every listed t degenerate
    python3 check_sigma.py FRESH.jsonl     also compare fresh records (from dm_sigma.py) with the committed ones

The committed records are dm_sigma_306.jsonl (every coprime ratio m <= 1000 with m^2 + n^2 a square or twice a
square, both curves, and the partner pairs outside that range) and dm_sigma_24.jsonl (the 24 ratios that were open
in the census before 24 November 2026 and four partners outside it, with ellrank effort up to 4 and Magma's upper
bounds of magma_sigma2_raw.txt)."""
import sys, os, json
from math import gcd, isqrt
D = os.path.dirname(os.path.abspath(__file__))
def load(f): return [json.loads(l) for l in open(os.path.join(D, f)) if l.strip()]
def partner(m, n): return (m + n, m - n) if (m + n) % 2 else ((m + n)//2, (m - n)//2)
def sq(x): r = isqrt(x); return r*r == x
KEYS = ('pair', 'curve', 'rank', 'status', 'amax', 'four')
bad = 0
recs = {}
for f in ('dm_sigma_306.jsonl', 'dm_sigma_24.jsonl'):
    for r in load(f):
        k = (tuple(r['pair']), r['curve'], f)
        assert k not in recs, k
        recs[k] = r
        if r['status'] == 'ALERT' or r['status'].startswith('error'): print('BAD', f, r['pair'], r['curve'], r['status']); bad += 1
        for t, v in r.get('four', {}).items():
            if v != 'degenerate': print('BAD', f, r['pair'], r['curve'], 't =', t, v); bad += 1
deg = [(m, n) for m in range(2, 1001) for n in range(1, m) if gcd(m, n) == 1 and (sq(m*m + n*n) or ((m*m + n*n) % 2 == 0 and sq((m*m + n*n)//2)))]
closed = lambda p, f: any(recs.get((p, c, f), {}).get('status') == 'closed' for c in ('E+', 'E-'))
for p in deg:
    for c in ('E+', 'E-'): assert (p, c, 'dm_sigma_306.jsonl') in recs, (p, c)
direct = [p for p in deg if closed(p, 'dm_sigma_306.jsonl')]
withp = [p for p in deg if closed(p, 'dm_sigma_306.jsonl') or closed(partner(*p), 'dm_sigma_306.jsonl')]
pyth = [p for p in deg if sq(p[0]**2 + p[1]**2)]
print('%d ratios with m <= 1000 and m^2 + n^2 a square (%d) or twice a square (%d): the four corners are squares only at '
      'the degenerate t for %d by their own curves, %d with the partner pair' % (len(deg), len(pyth), len(deg) - len(pyth), len(direct), len(withp)))
print('twice a square, closed by their own curves: %d' % len([p for p in direct if p not in pyth]))
left = [p for p in deg if p not in withp]
print('left (%d, %d of the first type): %s' % (len(left), len([p for p in left if p in pyth]), ' '.join('(%d,%d)' % p for p in left)))
OPEN24 = [(91,60),(151,31),(195,28),(209,120),(223,167),(329,89),(420,341),(440,279),(572,555),(608,105),(609,200),(627,364),
          (680,111),(713,503),(719,161),(741,580),(760,39),(761,79),(791,569),(799,721),(809,409),(861,620),(945,248),(991,263)]
c24 = [p for p in OPEN24 if closed(p, 'dm_sigma_24.jsonl') or closed(partner(*p), 'dm_sigma_24.jsonl')]
c24b = [p for p in OPEN24 if p in withp]
assert c24 == c24b, (c24, c24b)
print('of the 24 ratios open in the census before 24 November 2026: %d closed: %s' % (len(c24), ' '.join('(%d,%d)' % p for p in c24)))
print('open: %s' % ' '.join('(%d,%d)' % p for p in OPEN24 if p not in c24))
FOUR = [(72,65),(137,7),(425,168),(593,257),(740,539),(828,173),(840,83),(920,129),(923,757)]
c4 = [p for p in FOUR if p in withp]
print('of the 9 ratios closed in the four-point tier: %d closed here as well: %s' % (len(c4), ' '.join('(%d,%d)' % p for p in c4)))
if len(sys.argv) > 1:
    fresh = [json.loads(l) for l in open(sys.argv[1]) if l.startswith('{')]
    n = 0
    for r in fresh:
        k = (tuple(r['pair']), r['curve'], 'dm_sigma_24.jsonl')
        if k not in recs: k = (k[0], k[1], 'dm_sigma_306.jsonl')
        old = recs[k]
        for key in KEYS:
            if json.dumps(old.get(key), sort_keys=True) != json.dumps(r.get(key), sort_keys=True):
                print('DIFFERS', k, key, old.get(key), r.get(key)); bad += 1
        n += 1
    print('fresh records compared: %d' % n)
print('BAD %d' % bad if bad else 'OK')
sys.exit(1 if bad else 0)
