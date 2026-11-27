"""H1, part 5: the partner of each phantom (EXPERIMENT, written by an AI model; not read by a person; nothing here changes
the census).
The partner of m : n is (m + n) : (m - n), halved when both are even, in lowest terms (CENSUS9000.md); it swaps the
square and twice-square types, and the census closes a ratio when its partner closes.  For each phantom of the scan
(h1_scan_k*.txt) and G1's (1109503,597671), the partner goes through the same test as the scan (G1's one(): the class
criterion, then the four-class descent with the exact local test at every odd prime and good primes to 150, 1500, 5000,
20000 while more than the known tuples is left).  Since H3, good primes to 1153 are all that can cut, so "good primes to
20000" is the full local test.
Output h1_partners.txt: per phantom, its partner and whether the partner is a phantom, with its extra tuples.
usage: h1_partners.py [--write]
G1's folder is found through G1_DIR, or else ../../roadmap-g/g1."""
import sys, os, re, glob
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__))
G1 = os.environ.get('G1_DIR', os.path.join(HERE, '..', '..', 'roadmap-g', 'g1'))
if not os.path.isdir(G1): G1 = os.path.join(HERE, '..', 'g1')
sys.path.insert(0, G1)
argv = sys.argv; sys.argv = [argv[0]]
import g1_stress as G
sys.argv = argv

def partner(m, n):
    a, b = m + n, m - n
    g = gcd(a, b); return a // g, b // g

def phantoms():
    out = [(1109503, 597671)]
    for f in sorted(glob.glob(os.path.join(HERE, 'h1_scan_k*.txt'))):
        for l in open(f):
            mm = re.match(r'\((\d+),(\d+)\) PHANTOM', l)
            if mm: out.append((int(mm.group(1)), int(mm.group(2))))
    return sorted(set(out))

def main():
    ph = phantoms(); lines = ['H1: the partners of the phantoms. Written by an AI model; not read by a person. Experiment, not counted.']
    both = 0
    for m, n in ph:
        p = partner(m, n)
        st, d = G.one(*p)
        if d is None: lines.append('(%d,%d) -> partner (%d,%d): FAILED: %s' % (m, n, p[0], p[1], st)); continue
        S, known, pmax = d; extra = sorted(set(S) - set(known))
        if extra: both += 1
        lines.append('(%d,%d) -> partner (%d,%d)%s: %s; good primes to %d; known %s; extra %s' % (m, n, p[0], p[1],
                     ' (also in the list)' if p in ph else '', 'PHANTOM' if extra else 'no extra tuple', pmax, sorted(known), extra if extra else '[]'))
    lines.append('phantoms: %d; partners that are phantoms: %d' % (len(ph), both))
    return lines

if __name__ == '__main__':
    s = '\n'.join(main()) + '\n'; f = os.path.join(HERE, 'h1_partners.txt')
    if '--write' in argv: open(f, 'w').write(s)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
