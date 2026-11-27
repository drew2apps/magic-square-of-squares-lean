"""G4 input (EXPERIMENT, written by an AI model; not read by a person): every coprime square-type ratio m : n
(m^2 + n^2 a square, m > n) with m <= M, sorted, as {"open": [[m, n], ...]}.  No congruence filter.
usage: gen_square.py M   (checks square_<M>.json; --write rewrites it)"""
import sys, os, json, math
M = int([a for a in sys.argv[1:] if not a.startswith('--')][0])
S = set(); a = 2
while a*a - (a-1)**2 <= M:
    for b in range(1, a):
        if (a - b) % 2 and math.gcd(a, b) == 1:
            x, y = a*a - b*b, 2*a*b
            if max(x, y) <= M: S.add((max(x, y), min(x, y)))
    a += 1
s = json.dumps({'open': [list(p) for p in sorted(S)]})
f = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'square_%d.json' % M)
if '--write' in sys.argv: open(f, 'w').write(s)
print('%d ratios' % len(S))
if '--write' not in sys.argv: assert json.load(open(f)) == json.loads(s); print('# matches square_%d.json' % M)
