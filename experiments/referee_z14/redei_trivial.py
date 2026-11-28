"""z14: how many of mixtest's soluble Redei pairs are trivial (d = beta)?  The sampled pairs are recovered by replaying
mixtest.run's random draws (same seed, same calls in the same order; checked on the counts it recorded).  Written by an AI model (reviewer z14); not read
by a person.   usage: python3 -I redei_trivial.py | grep -v '^PARI' | diff - redei_trivial.txt"""
import os, sys, json, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv = sys.argv[:1]
sys.path.insert(0, os.path.join(HERE, '..', 'conductor-t10'))
import t10_pairct as T
T.PARI.allocatemem(10 ** 8, 4 * 10 ** 8)
from fractions import Fraction as Fr
RI = T.RI
def sqf(x): x = Fr(x); return RI.sqfree_part(x.numerator * x.denominator)
R = {}
for f in sorted(glob.glob(os.path.join(HERE, 'mix_*.jsonl'))):
    for l in open(f): r = json.loads(l); R[(r['m'], r['n'])] = r
import random
NA_ALL, NREDEI = 24, 150
def factor(n):
    n = abs(int(n))
    if n < 2: return {}
    f = T.PARI.factor(n); return {int(q): int(e) for q, e in zip(f[0], f[1])}
def replay(m, n, k, Vmask):
    rng = random.Random(1000003 * m + n)
    for _ in range(3 if k > 1 else 0): rng.sample(range(k), 2)
    if 2 ** k - 1 <= NA_ALL: masks = list(range(1, 2 ** k))
    else:
        masks = set(1 << i for i in range(k))
        while len(masks) < NA_ALL: masks.add(rng.randrange(1, 2 ** k))
        masks = sorted(masks)
    for _ in range(12): rng.choice(masks); rng.randrange(k)
    nb = [x for x in masks if bin(x).count('1') > 1]
    if nb: rng.sample(nb, min(2, len(nb)))
    gAB = [-1] + sorted(factor(2 * (m * m - n * n) * 2 * m * n))
    rng.sample(range(2 ** len(gAB)), min(24, 2 ** len(gAB)))
    pairs = [(x, y) for x in Vmask for y in Vmask if x and y]
    if len(pairs) > NREDEI: pairs = rng.sample(pairs, NREDEI)
    return pairs, len(masks)
exact = bound = 0; sol_tot = 0; nmask_ok = 0
for key in sorted(R):
    r = R[key]; sel = r['sel']; k = len(sel)
    def elt(mask):
        c1 = c2 = 1
        for i in range(k):
            if mask >> i & 1: c1 *= sel[i][0]; c2 *= sel[i][1]
        return (sqf(c1), sqf(c2))
    Vmask = [mask for mask in range(2 ** k) if mask == 0 or elt(mask)[0] == 1]
    pairs, na = replay(key[0], key[1], k, Vmask)
    nmask_ok += na == r['na']
    sol_tot += r['redei'].get('soluble', 0)
    for x, y in pairs:
        if x == y:
            d = elt(x)[1]      # x^2 - d y^2 = d z^2 soluble iff (d, -1)_v = 1 at every v
            ps = sorted(set(factor(2 * d)) | {2})
            exact += all(T.RS.hilbert(d, -1, p) == 1 for p in ps) and T.RS.hilbert(d, -1, 'inf') == 1
    bound += len(pairs) != (r['redei'].get('soluble', 0) + r['redei'].get('insoluble', 0))
print('replay: number of a-classes as recorded on %d of %d curves; number of Redei pairs as recorded on %d of %d curves' % (nmask_ok, len(R), len(R) - bound, len(R)))
print('soluble Redei pairs %d; with d = beta %d; with d != beta %d' % (sol_tot, exact, sol_tot - exact))
