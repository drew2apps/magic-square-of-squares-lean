"""H3 (part 3): the shape of the allowed sets of h3_weil.json (EXPERIMENT, written by an AI model; not read by a person).
For the good primes p = 1 mod 4 with 17 <= p <= 1153 (91 of them; the session wrote 96) and every good residue r = m/n mod p: how many sets of each size,
whether each is an affine subspace of F_2^4 (a coset of a subgroup), and, for each size s, the largest prime with a set
of at most s vectors.  Every set is closed under v -> v + (1,1,1,1) (the class of a non-square unit), so sizes are even.
usage: h3_sizes.py [--write]   (h3_sizes.txt; reads h3_weil.json)"""
import sys, os, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
T = {int(p): {int(r): a for r, a in d.items()} for p, d in json.load(open(os.path.join(HERE, 'h3_weil.json'))).items()}
def vecs(a): return [b for b in range(16) if a >> b & 1]
def affine(a):
    v = vecs(a)
    if not v: return True
    sub = set(x ^ v[0] for x in v)
    return all((x ^ y) in sub for x in sub for y in sub)
size = collections.Counter(); aff = collections.Counter(); mx = {}
for p in sorted(T):
    for r, a in T[p].items():
        c = len(vecs(a)); size[c] += 1; aff[affine(a)] += 1
        assert all((b ^ 15) in vecs(a) for b in vecs(a)), 'a set not closed under the unit class'
        for s in range(0, 17, 2):
            if c <= s: mx[s] = max(mx.get(s, 0), p)
out = ['H3 (part 3): the allowed sets at the good primes 17 <= p <= 1153, p = 1 mod 4. Written by an AI model; not read by a person. Experiment, not counted.',
       'every set is closed under adding (1,1,1,1): yes',
       'sets by size: ' + ', '.join('%d: %d' % (s, size[s]) for s in sorted(size)),
       'affine subspaces of F_2^4: %d of %d (so the good-prime test is not linear over F_2)' % (aff[True], sum(aff.values())),
       'largest prime with a set of at most s vectors: ' + ', '.join('%d: %d' % (s, mx[s]) for s in sorted(mx) if s < 16)]
s = '\n'.join(out) + '\n'; f = os.path.join(HERE, 'h3_sizes.txt')
if '--write' in sys.argv: open(f, 'w').write(s)
else:
    sys.stdout.write(s)
    if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
