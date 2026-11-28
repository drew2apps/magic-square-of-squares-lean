"""Reviewer z16 (an AI model; not read by a person): EXPLORATION_5 §3's ghost data, checked against two tier-A filters
that the four-form pipeline does not apply: the 484 residue classes (scripts/g1_full/modp_classes.txt, five forms
modulo 34 primes; Lean obstructsFull) and "every prime = 5 mod 8 divides m^2 + n^2 to an even power" (G1Primes.lean,
IsMSS.edge_norm).  Inputs: data/out_*.jsonl (EXPLORATION_4, the 937), data/samp_*.jsonl (EXPLORATION_5's three
samples beyond m = 15000), out_samples.jsonl and out937.jsonl (z16_five.py).  usage: python3 -I z16_ghost.py"""
import os, json, collections
from sympy import factorint
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
CL = collections.defaultdict(set)
for l in open(os.path.join(REPO, 'scripts', 'g1_full', 'modp_classes.txt')):
    p, r = map(int, l.split()); CL[p].add(r)
def in484(m, n):
    return sorted(p for p in CL if n % p and (m * pow(n, -1, p)) % p in CL[p])
def bad5(m, n): return sorted(q for q, e in factorint(m*m + n*n).items() if q % 8 == 5 and e % 2)
def sqf(x):
    r = 1
    for q, e in factorint(abs(x)).items():
        if e % 2: r *= q
    return r
def shape(m, n, t):
    a = sqf(m*n*t[0]*t[1]); b = sqf((m*m - n*n)*t[2]*t[3])
    if (a, b) == (sqf(2*m*n), 1): return 'S1'
    if (a, b) == (1, sqf(2*(m*m - n*n))): return 'S2'
    return 'other'
def load(f): return [json.loads(l) for l in open(os.path.join(HERE, f)) if l.startswith('{')]
pops = {'937 (m <= 15000)': load('data/out_997.jsonl') + load('data/out_15000_0.jsonl') + load('data/out_15000_1.jsonl')}
for s in ('biased', 'random', 'fewp'): pops['sample ' + s] = [r for r in load('data/samp_%s.jsonl' % s) if 'tuples' in r]
five = {}
for f in ('out937.jsonl', 'out_samples.jsonl'):
    if os.path.exists(os.path.join(HERE, f)):
        for r in load(f): five[(r['m'], r['n'])] = r
allsurv = []
for name, rows in pops.items():
    c = collections.Counter(); surv = []
    for r in rows:
        m, n = r['m'], r['n']; k = (m, n)
        a, b = bool(in484(m, n)), bool(bad5(m, n))
        c['ratios'] += 1; c['in the 484 classes'] += a; c['prime = 5 mod 8 to odd power in C'] += b; c['either'] += a or b
        if k in five:
            c['five-form run'] += 1; c['keep a 5-tuple'] += bool(five[k]['n5'])
            c['keep a 5-tuple and pass both filters'] += bool(five[k]['n5']) and not (a or b)
        for t, kk in r['tuples']:
            c['tuples'] += 1
            if not any(v == 1 for v in kk.values()):
                surv.append((k, tuple(t), shape(m, n, t), a, b))
    c['pairing survivors'] = len(surv); c['ratios with a survivor'] = len({s[0] for s in surv})
    print('%-18s %s' % (name, dict(c)))
    allsurv += surv
print('survivors (ratio, tuple, shape, in the 484 classes, prime = 5 mod 8 to odd power in C, extends to a 5-tuple):')
for k, t, sh, a, b in allsurv:
    ext = (list(t) in five[k]['ext']) if k in five else None
    print('   %s %s %s 484:%s 5mod8:%s ext:%s  C = %s' % (k, t, sh, a, b, ext, dict(factorint(k[0]**2 + k[1]**2))))
print('partner (24599,5207) (not in the data files): 484 classes %s, primes = 5 mod 8 of odd exponent %s, C = %s'
      % (in484(24599, 5207), bad5(24599, 5207), dict(factorint(24599**2 + 5207**2))))
print('survivor ratios: %d (+ the partner); shapes %s' % (len({s[0] for s in allsurv}), dict(collections.Counter(s[2] for s in allsurv))))
# the 937 restricted to the ratios that pass both tier-A filters (the population that matters for a magic square)
rows = pops['937 (m <= 15000)']
keep = [r for r in rows if not in484(r['m'], r['n']) and not bad5(r['m'], r['n'])]
k5 = [r for r in keep if five[(r['m'], r['n'])]['n5']]
print('937 passing the 484 classes and the 5 mod 8 condition: %d ratios, %d four-form tuples, %d pairing survivors; '
      'with F4: %d ratios keep a 5-tuple, %d 4-tuples extend, %d 5-tuples'
      % (len(keep), sum(len(r['tuples']) for r in keep), sum(1 for r in keep for t, kk in r['tuples'] if not any(v == 1 for v in kk.values())),
         len(k5), sum(len(five[(r['m'], r['n'])]['ext']) for r in k5), sum(five[(r['m'], r['n'])]['n5'] for r in k5)))
print('   by the 484 classes alone: %d of 937 excluded; by the 5 mod 8 condition alone: %d; by the 484 classes, among the 211 that keep a 5-tuple: %d'
      % (sum(1 for r in rows if in484(r['m'], r['n'])), sum(1 for r in rows if bad5(r['m'], r['n'])),
         sum(1 for r in rows if five[(r['m'], r['n'])]['n5'] and in484(r['m'], r['n']))))
