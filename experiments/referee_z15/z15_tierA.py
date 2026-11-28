"""z15: which of the 937 ratios with residual tuples (and which of the 5 survivor ratios) lie in one of the 484 residue
classes of scripts/g1_full/modp_classes.txt (m = r n mod p, n not 0 mod p, as scripts/search/ratio_pairs.c tests F3), or
have a prime p = 5 mod 8 dividing m^2 + n^2 to an odd power and not dividing n (the hypothesis of G1Primes.lean).
Membership is computed here; the Lean status of the two exclusions is the repository's (not re-checked).
Written by an AI model (reviewer z15), not read by a person."""
import os, json
import cypari2
P = cypari2.Pari(); HERE = os.path.dirname(os.path.abspath(__file__))
CL = [tuple(map(int, l.split())) for l in open(os.path.join(HERE, '..', '..', 'scripts', 'g1_full', 'modp_classes.txt')) if l.strip()]
assert len(CL) == 484
def in484(m, n): return [(p, r) for p, r in CL if n % p and m % p == r * n % p]
def p5(m, n):
    f = P.factor(m * m + n * n)
    return [int(q) for q, e in zip(f[0], f[1]) if int(q) % 8 == 5 and int(e) % 2 and n % int(q)]
L = [(r['m'], r['n']) for r in map(json.loads, open(os.path.join(HERE, 'z15_residual_all.jsonl')))]
a = [r for r in L if in484(*r)]; b = [r for r in L if p5(*r)]; u = set(a) | set(b)
print('z15_tierA.py: written by an AI model (reviewer z15), not read by a person.')
print('937 ratios: in a 484 class %d; p = 5 mod 8 to an odd power in m^2+n^2 %d; either %d; neither %d' % (len(a), len(b), len(u), len(L) - len(u)))
for r in [(3320, 2303), (5623, 1017), (9113, 3200), (12313, 5913), (14903, 9696)]:
    print('  survivor ratio', r, ': 484 classes', in484(*r), '; p = 5 mod 8:', p5(*r))
for r in [(627, 364), (991, 263)]:
    print('  degenerate control', r, ': 484 classes', in484(*r), '; p = 5 mod 8:', p5(*r))
# the ratios not excluded by either: how many has the repository's t10 decided (pilot or my reruns)?
import re
pilot = set()
for f in ('t10_997_a.txt', 't10_997_b.txt', 't10_15000_s100_a.txt', 't10_15000_s100_b.txt'):
    for l in open(os.path.join(HERE, '..', 'conductor-t10', f)):
        mm = re.match(r'ratio \((\d+),(\d+)\)', l)
        if mm: pilot.add((int(mm.group(1)), int(mm.group(2))))
mine = set()
for f in os.listdir(HERE):
    if f.startswith('t10_') and f.endswith('.jsonl'):
        for l in open(os.path.join(HERE, f)): r = json.loads(l); mine.add((r['m'], r['n']))
rest = [r for r in L if r not in u]
print('of the %d not excluded: in the pilot %d, in my t10 reruns (all t10_*.jsonl now) %d, in neither %d' % (len(rest), sum(r in pilot for r in rest), sum(r in mine for r in rest), sum(r not in pilot and r not in mine for r in rest)))
import random
# notA_order.txt (input of the t10 run t10_notA.jsonl) was drawn once, when the reruns were t10_five, t10_ghost (stopped
# after 46 ratios) and the first ~61 of sample80: the 49 then undecided by t10, in the order of seed 71.  Kept as is.
if not os.path.exists(os.path.join(HERE, 'notA_order.txt')):
    left = sorted(r for r in rest if r not in pilot and r not in mine)
    random.Random(71).shuffle(left)
    open(os.path.join(HERE, 'notA_order.txt'), 'w').write(''.join('%d,%d\n' % r for r in left))
