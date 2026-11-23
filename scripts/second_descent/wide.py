#!/usr/bin/env python3
"""Both descents on the ratios with a degenerate solution beyond the census.  EXPERIMENT.
For M0 < m <= M1: every coprime pair that avoids the 484 residue classes and has m^2 + n^2 a square or twice a square.
For each: the first descent (must leave exactly 2 tuples, those of the two degenerate solutions, which must differ), then the
second descent; if that does not close the ratio, the same for its partner.  One JSON line per ratio.
usage: wide.py M0 M1 classes_file        (classes_file: ../g1_full/modp_classes.txt)
       wide.py --summary FILE [FILE ...]"""
import sys, json, time
def summary(files):
    R = [json.loads(l) for f in files for l in open(f)]
    st = {}
    for r in R: st[r['status']] = st.get(r['status'], 0) + 1
    print('%d ratios with a degenerate solution that avoid the 484 classes, m from %d to %d' % (len(R), min(r['pair'][0] for r in R), max(r['pair'][0] for r in R)))
    print('first descent: exactly 2 tuples left, those of the two degenerate solutions, for %d of %d; the same for all %d partners that were run' % (
        sum(r.get('first') == 2 and r.get('differ') for r in R), len(R), sum('partner' in r and r['partner'].get('first') == 2 and r['partner'].get('differ') for r in R)))
    print('second descent:', ', '.join('%s %d' % kv for kv in sorted(st.items())))
    ts = set(t + ': ' + v for r in R for t, v in list(r.get('rational_t', {}).items()) + list(r.get('partner', {}).get('rational_t', {}).items()))
    print('rational t returned by the equations:', ', '.join(sorted(ts)))
    print('not closed (%d):' % st.get('open', 0), ' '.join('%d,%d' % tuple(r['pair']) for r in R if r['status'] == 'open'))
if len(sys.argv) > 1 and sys.argv[1] == '--summary': summary(sys.argv[2:]); sys.exit(0)
from math import gcd, isqrt
import numpy as np
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gauss_descent'))
import gauss_descent as gd, scan
from sage.all__sagemath_schemes import *
import second_descent as sd

K = sd.K; I = sd.I
def trivial(m, n):
    S = m*m + n*n
    return isqrt(S)**2 == S or (S % 2 == 0 and isqrt(S//2)**2 == S//2)

def first(m, n):
    """Number of tuples left by the first descent, and whether the two degenerate solutions have different tuples."""
    c, last = scan.descent_fast(m, n, pmax=200)
    S = m*m + n*n
    sols = [(1, 1), (1, -1)] if isqrt(S)**2 == S else [(0, 1), (1, 0)]
    forms = [(I*n, K(m)), (I*m, K(n)), (K(1), I), (m + I*n, -n - I*m), (m + I*n, n + I*m)]
    vals = [[A*u + B*v for A, B in forms] for u, v in sols]
    differ = not all((a/b).is_square() for a, b in zip(*vals))
    return c, differ

def one(m, n):
    c, differ = first(m, n); out = dict(pair=[m, n], first=c, differ=differ)
    if c != 2 or not differ: out['status'] = 'first descent: not the expected 2 tuples'; return out
    r = sd.close(m, n)
    out['status'] = r['status']; out['by'] = r.get('by')
    if r['status'] == 'closed':
        cc = r['curves'][-1]; out['rank_bounds'] = cc['rank_bounds']; out['rational_t'] = cc['rational_t']
    return out

if __name__ == '__main__':
    M0, M1 = int(sys.argv[1]), int(sys.argv[2]); cl = {}
    for l in open(sys.argv[3]):
        if l.startswith('#') or not l.strip(): continue
        a = l.split(); cl.setdefault(int(a[0]), set()).update(int(x) for x in a[1:])
    for m in range(M0 + 1, M1 + 1):
        for n in range(1, m):
            if gcd(m, n) != 1 or not trivial(m, n): continue
            if any(n % p and (m*pow(n, -1, p)) % p in rs for p, rs in cl.items()): continue
            t0 = time.time()
            try:
                o = one(m, n)
                if o['status'] == 'open':
                    a, b = m + n, m - n
                    if a % 2 == 0: a, b = a//2, b//2
                    po = one(a, b); o['partner'] = po; 
                    if po['status'] == 'closed': o['status'] = 'closed by partner'
            except Exception as ex:
                o = dict(pair=[m, n], status='error %s %s' % (type(ex).__name__, str(ex)[:80]))
            print(json.dumps(o, default=str), flush=True)
