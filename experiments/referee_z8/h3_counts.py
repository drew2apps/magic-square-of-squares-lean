"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
Claim 3 of the review: the numbers of REPORT_H3.md recounted from H3's committed h3_weil.json and h3_tables.json, with
my own code (and sympy's prime list for the number of primes).
usage: h3_counts.py [--write]   (h3_counts.txt; under 1 s)"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit, H3
from sympy import primerange
T = {int(p): {int(r): a for r, a in d.items()} for p, d in json.load(open(os.path.join(H3, 'h3_weil.json'))).items()}
J = {int(p): {int(r): a for r, a in d.items()} for p, d in json.load(open(os.path.join(H3, 'h3_tables.json'))).items()}
FULL = (1 << 16) - 1; ps = sorted(T)
def vecs(a): return [b for b in range(16) if a >> b & 1]
def affine(a):
    v = vecs(a)
    if not v: return True
    sub = set(x ^ v[0] for x in v)
    return all((x ^ y) in sub for x in sub for y in sub)
want = [p for p in primerange(17, 1154) if p % 4 == 1]
out = ['z8: REPORT_H3 numbers recounted from h3_weil.json / h3_tables.json. Written by an AI model; not read by a person.',
       'primes in h3_weil.json: %d (%d to %d); primes p = 1 mod 4 with 17 <= p <= 1153 (sympy): %d; equal lists: %s; with 5 and 13: %d' % (len(ps), ps[0], ps[-1], len(want), ps == want, len([p for p in primerange(2, 1154) if p % 4 == 1])),
       'good residues: %d; losing a vector: %d' % (sum(len(T[p]) for p in ps), sum(1 for p in ps for a in T[p].values() if a != FULL)),
       'primes with no loss: %s' % [p for p in ps if all(a == FULL for a in T[p].values())],
       'sets that are not affine subspaces: %d' % sum(1 for p in ps for a in T[p].values() if not affine(a)),
       'empty sets: %d, largest p with one: %d' % (sum(1 for p in ps for a in T[p].values() if a == 0), max(p for p in ps for a in T[p].values() if a == 0)),
       'largest p with a set of at most s vectors: ' + ', '.join('%d: %d' % (s, max(p for p in ps for a in T[p].values() if bin(a).count('1') <= s)) for s in (2, 4, 6, 8, 10, 12, 14)),
       'most vectors lost at a prime above 1069: %d' % max(16 - bin(a).count('1') for p in ps if p > 1069 for a in T[p].values()),
       'at 1153: %d good residues, %d lose a vector' % (len(T[1153]), sum(1 for a in T[1153].values() if a != FULL)),
       'mean allowed size: ' + ', '.join('%d: %.2f' % (p, sum(bin(a).count('1') for a in T[p].values()) / len(T[p])) for p in (17, 193, 197, 509, 1153)),
       'h3_tables.json: primes %d to %d (%d); good residues %d; empty %d; full sets at %s; no good residue at %s' % (min(J), max(J), len(J), sum(len(J[p]) for p in J), sum(1 for p in J for a in J[p].values() if a == 0), {p: sum(1 for a in J[p].values() if a == FULL) for p in sorted(J) if any(a == FULL for a in J[p].values())}, [p for p in sorted(J) if not J[p]]),
       'good residues at 197 < p < 1154: %d' % sum(len(T[p]) for p in ps if p > 197)]
emit(out, 'h3_counts.txt')
