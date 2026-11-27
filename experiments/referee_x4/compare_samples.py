#!/usr/bin/env python3
"""Review x4 (run 6): compares the item-2 reruns in samples_11k/ (made by samples_11k.sh with the repository's
programs, PARI 2.17.4) with the committed records of the 9000 < m <= 11000 run.  Prints only; ends with assertions.
Written by an AI model; not read by a person.   python3 compare_samples.py"""
import os, json, collections
HERE = os.path.dirname(os.path.abspath(__file__)); S = os.path.join(HERE, 'samples_11k')
CW = os.path.abspath(os.path.join(HERE, '..', 'census-wide'))
def jl(p):
    for l in open(p):
        if l.startswith('{'): yield json.loads(l)
def strip(o):
    o = dict(o); o.pop('sec', None)
    if 'partner' in o: o['partner'] = strip(o['partner'])
    return o
ok = True
old = {tuple(o['pair']): o for k in (0, 1) for o in jl(os.path.join(CW, 's11_%d.jsonl' % k))}
new = [o for k in (0, 1) for o in jl(os.path.join(S, 'sd_%d.jsonl' % k))]
same = sum(strip(o) == strip(old[tuple(o['pair'])]) for o in new); ok &= same == len(new) == 10
print('second descent (wide_pari.py): %d of %d records identical apart from "sec"; statuses %s' % (same, len(new), dict(collections.Counter(o['status'] for o in new))))
a = {}
for k in (0, 1): a.update(json.load(open(os.path.join(CW, 'g11k_out%d.json' % k))))
b = json.load(open(os.path.join(S, 'g_out.json')))
same = sum(b[k] == a[k] for k in b); ok &= same == len(b) == 10
print('rank 0 (global.py): %d of %d records identical; %d close' % (same, len(b), sum(any(c['closes'] for c in v) for v in b.values())))
old = {(tuple(o['of']), tuple(o['pair']), o['curve']): o for k in (0, 1) for o in jl(os.path.join(CW, 'dm11k%d.jsonl' % k))}
new = {(tuple(o['of']), tuple(o['pair']), o['curve']): o for o in jl(os.path.join(S, 'dm_out.jsonl'))}
st = all(new[k]['status'] == old[k]['status'] for k in new)
cl = [k for k in new if new[k]['status'] == 'closed']
clsame = all(strip(new[k]) == strip(old[k]) for k in cl)
lower = sum(1 for k in new if new[k]['status'] != 'closed' and new[k]['rank'][1] == old[k]['rank'][1] and new[k]['rank'][0] != old[k]['rank'][0])
upper = sum(1 for k in new if new[k]['rank'][1] != old[k]['rank'][1])
ok &= st and clsame and upper == 0 and len(new) == 80
print('Demjanenko-Manin (dm_pari.py): %d records; every status identical: %s; the %d closed records identical apart from "sec": %s; '
      'upper rank bounds identical: %s; %d "rank not 1" records differ only in PARI\'s lower bound' % (len(new), st, len(cl), clsame, upper == 0, lower))
def last(f): return [json.loads(l) for l in open(f) if l.startswith('{')][-1]
for d in ('220962238', '-220962238'):
    x = strip(last(os.path.join(S, 'route_10512_%s.txt' % d))); y = strip(last(os.path.join(CW, 'route', 'log_10512_%s.txt' % d)))
    ok &= x == y and x['closed'] is False and x['rank_upper_isogeny_class'] == 3
    print('route (route_pari.py) (10512,145), d = %s: identical apart from "sec": %s; closed %s; rank_E1 %s; isogeny-class upper bound %s; left %s'
          % (d, x == y, x['closed'], x['rank_E1'], x['rank_upper_isogeny_class'], x['left']))
    r1 = [l.split()[3:12] for l in open(os.path.join(S, 'rh_10512_%s.txt' % d))]; r0 = [l.split()[3:12] for l in open(os.path.join(CW, 'route', 'rh_10512_%s.txt' % d)) if not l.startswith('PARI')]
    ok &= r1 == r0
    print('rank_hard.py effort 8, d = %s: %s (committed: %s)' % (d, [' '.join(r) for r in r1], 'same' if r1 == r0 else 'DIFFERENT'))
cd = os.path.abspath(os.path.join(HERE, '..', '..', 'scripts', 'cover_descent'))
L1 = open(os.path.join(S, 'list_10000_11000.txt')).read().splitlines()
L0 = [l for l in open(os.path.join(cd, 'structure', 'class_selmer_15000_list.txt')).read().splitlines() if 10000 < int(l.split(',')[0]) <= 11000]
T1 = open(os.path.join(S, 'tuple_10000_11000.txt')).read().splitlines()
T0 = [l for l in open(os.path.join(cd, 'tuple_selmer_15000.txt')).read().splitlines() if 10000 < int(l[1:].split(',')[0]) <= 11000]
ok &= L1 == L0 and T1 == T0 and all(l.endswith('nothing is left') for l in T1)
print('class step 10000 < m <= 11000 (sympy): list identical to the committed one: %s (%d ratios); four-class descent identical: %s; nothing left for %d of %d'
      % (L1 == L0, len(L1), T1 == T0, sum(l.endswith('nothing is left') for l in T1), len(T1)))
print('all samples reproduce: %s' % ok); assert ok
