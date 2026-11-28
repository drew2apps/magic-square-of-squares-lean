"""z15: section 3's "reason": at each odd bad prime p, is the allowed set of exact_local.allowed_at (label tuples in
(F_2^2)^4 or (F_2)^4; dlabel is a homomorphism there) a coset of a subspace?  If every local set were a coset (and the
set at 2 too), the residual set -- an intersection of preimages under homomorphisms -- would be a coset or empty.
Counts over the 937 ratios; p = 2 (labels of tuple_descent's D2) is not tested.  Written by an AI model (reviewer
z15), not read by a person."""
import os, sys, json, collections
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
ns = {'__file__': os.path.join(HERE, 'z15_residual_all.py'), '__name__': 'z15res'}
exec(open(os.path.join(HERE, 'z15_residual_all.py')).read().split('def read_list')[0], ns)
X = ns['X']; TD = ns['TD']; fact = ns['fact']
def bits(labs):
    v = 0
    for lab in labs:
        for b in lab: v = (v << 1) | b
    return v
def is_coset(S):
    S = list(S); v0 = S[0]; basis = []
    for x in S:
        x ^= v0
        for b in basis:
            if x ^ b < x: x ^= b
        if x: basis.append(x); basis.sort(reverse=True)
    return 2 ** len(basis) == len(S)
res = [json.loads(l) for l in open(os.path.join(HERE, 'z15_residual_all.jsonl'))]
coset_of = {}
import cypari2
P = cypari2.Pari()
def vec(t):
    v = set()
    for i, d in enumerate(t):
        f = P.factor(d)
        for q, e in zip(f[0], f[1]):
            if int(e) % 2: v.add((i, int(q)))
    return frozenset(v)
loc_stats = collections.Counter(); ratio_stats = collections.Counter()
for r in res:
    m, n = r['m'], r['n']; V = [vec(t) for t in r['tuples']]
    # residual coset test as in z15_shapes.py
    idx = {}; rows = []
    for v in V:
        x = 0
        for e in v ^ V[0]:
            if e not in idx: idx[e] = len(idx)
            x |= 1 << idx[e]
        rows.append(x)
    basis = []
    for x in rows:
        for b in basis:
            if x ^ b < x: x ^= b
        if x: basis.append(x); basis.sort(reverse=True)
    rc = 2 ** len(basis) == len(V)
    forms, sup, known = TD.magic(m, n)
    bad = set()
    for x in (2 * m * n * (m * m - n * n), m * m + n * n, m ** 8 + 14 * m ** 4 * n ** 4 + n ** 8): bad |= set(fact(x))
    noncoset_odd = []
    for p in sorted(bad):
        if p == 2: continue
        L, al = X.allowed_at(forms, p)
        if len(al) == len(L.labels()) ** len(forms): continue
        c = is_coset(set(bits(a) for a in al))
        loc_stats[c] += 1
        if not c: noncoset_odd.append(p)
    ratio_stats[('residual coset' if rc else 'residual NOT coset', 'an odd local set is not a coset' if noncoset_odd else 'all odd local sets are cosets')] += 1
print('z15_localcoset.py: written by an AI model (reviewer z15), not read by a person.')
print('odd bad primes with a non-trivial allowed set, over the 937 ratios: coset %d, not a coset %d' % (loc_stats[True], loc_stats[False]))
for k, v in sorted(ratio_stats.items()): print('   ', k, v)
