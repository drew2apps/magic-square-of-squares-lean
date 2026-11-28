"""z15: root numbers of E_08, E_LU, E_mix (models checked isomorphic to t10's curves in z15_models.txt) with PARI's
ellrootno, on (1) the session's populations, reproduced draw for draw (rootno.py's random.Random(1) sequence), to check
its counts; (2) complete populations instead of samples: all 997 and all 8,739 class-keeping ratios of the G2 lists, the
937 with residual tuples (101 + 836), the 5 survivor ratios; (3) the 2-Selmer parity check: (-1)^(dim Sel2 - 2) = w
(Monsky 1996; Dokchitser-Dokchitser: the parity of the 2-Selmer rank minus dim E[2](Q) is the root number for E/Q, with
full rational 2-torsion as here) against the session's own dims in out_*.jsonl and t10's dims in my reruns.
Statistics: chi-square of the four patterns against 1/4 each (3 d.f.), z of w_mix = +1 against 1/2, and the phi
coefficient of (w_08, w_LU) with z = phi sqrt(N).  Written by an AI model (reviewer z15), not read by a person."""
import os, json, re, random, collections, math
from math import gcd, isqrt
import cypari2
P = cypari2.Pari(); P.allocatemem(10 ** 9)
HERE = os.path.dirname(os.path.abspath(__file__)); G2 = os.path.join(HERE, '..', 'roadmap-g', 'g2')
def cubic(e): return [0, -(e[0] + e[1] + e[2]), 0, e[0] * e[1] + e[0] * e[2] + e[1] * e[2], -e[0] * e[1] * e[2]]
CACHE = {}
def w(m, n):
    if (m, n) not in CACHE:
        A, B, C = m * m - n * n, 2 * m * n, m * m + n * n
        CACHE[(m, n)] = tuple(int(P.ellrootno(P.ellinit(cubic(e)))) for e in
                              ((0, A * A, C * C), (0, B * B, C * C), (0, 2 * (C * C - A * B), 2 * (C * C + A * B))))
    return CACHE[(m, n)]
def tab(pairs, label):
    N = len(pairs); c = collections.Counter(); cm = collections.Counter()
    for m, n in pairs:
        a, b, x = w(m, n); c[(a, b)] += 1; cm[x] += 1
    chi = sum((c[k] - N / 4) ** 2 / (N / 4) for k in [(1, 1), (1, -1), (-1, 1), (-1, -1)])
    # phi coefficient of the two signs
    pa = sum(v for k, v in c.items() if k[0] == 1) / N; pb = sum(v for k, v in c.items() if k[1] == 1) / N
    p11 = c[(1, 1)] / N
    den = math.sqrt(pa * (1 - pa) * pb * (1 - pb)) if 0 < pa < 1 and 0 < pb < 1 else float('nan')
    phi = (p11 - pa * pb) / den if den == den and den > 0 else float('nan')
    zmix = (cm[1] - N / 2) / (math.sqrt(N) / 2)
    print('%s: N = %d | (w08,wLU) (+,+) %d (+,-) %d (-,+) %d (-,-) %d | chi2(3) = %.2f | phi = %+.3f (z = %+.2f) | w_mix +1: %d, -1: %d (z = %+.2f)'
          % (label, N, c[(1, 1)], c[(1, -1)], c[(-1, 1)], c[(-1, -1)], chi, phi, phi * math.sqrt(N), cm[1], cm[-1], zmix))
print('z15_rootno.py: written by an AI model (reviewer z15), not read by a person.')
print('chi2 with 3 d.f.: 7.81 at p = 0.05, 11.34 at p = 0.01.')
print()
print('(1) the session\'s populations, drawn as rootno.py draws them:')
def load(fn): return [tuple(int(x) for x in l.split()[0].split(',')) for l in open(fn) if re.match(r'\d+,\d+', l)]
def is_sq(x): r = isqrt(x); return r * r == x
K1 = load(os.path.join(G2, 'g2_997_all.txt')); K2 = load(os.path.join(G2, 'g2_15000_all.txt'))
tab(K1, 'class-keeping non-degenerate, m<=5000')
rng = random.Random(1); tab(rng.sample(K2, 2500), 'class-keeping non-degenerate, 5000<m<=15000 (sample)')
need = [tuple(int(x) for x in l.split()[0].split(',')) for l in open(os.path.join(G2, 'g2_997_all.txt')) if 'left' in l]
tab(need, '  of which bad primes leave tuples (m<=5000)')
def rand(H, N, cond):
    out = []
    while len(out) < N:
        m = rng.randint(2, H); n = rng.randint(1, m - 1); S = m * m + n * n
        if gcd(m, n) == 1 and not is_sq(S) and not (S % 2 == 0 and is_sq(S // 2)) and cond(m, n): out.append((m, n))
    return out
tab(rand(15000, 2500, lambda m, n: True), 'random non-degenerate, m<=15000')
tab(rand(15000, 2500, lambda m, n: (m * n) % 4 == 0 or (m - n) % 8 == 0 or (m + n) % 8 == 0), 'random non-degenerate with a 2-adic point')
print()
print('(2) complete populations:')
L937 = [tuple(json.loads(l)[k] for k in ('m', 'n')) for l in open(os.path.join(HERE, 'z15_residual_all.jsonl'))]
tab(K2, 'all class-keeping, 5000<m<=15000')
tab(K1 + K2, 'all class-keeping, m<=15000')
tab(L937, 'the 937 with residual tuples')
tab([r for r in L937 if r[0] > 5000], '  of which 5000<m<=15000 (836)')
tab([r for r in K1 + K2 if r not in set(L937)], 'class-keeping without residual tuples (8,799)')
S5 = [(3320, 2303), (5623, 1017), (9113, 3200), (12313, 5913), (14903, 9696)]
tab(S5, 'the 5 survivor ratios')
print('   per survivor ratio (w08, wLU, wmix):', {r: w(*r) for r in S5})
print()
print('(3) 2-Selmer parity, (-1)^(dim Sel2 - 2) against w:')
rows = [json.loads(l) for f in ('out_997.jsonl', 'out_15000_0.jsonl', 'out_15000_1.jsonl') for l in open(os.path.join(HERE, 'session_data', f))]
bad = collections.Counter(); tot = collections.Counter()
for r in rows:
    W = w(r['m'], r['n'])
    for cv, i in (('08', 0), ('LU', 1), ('0L', 2), ('0U', 2), ('8L', 2), ('8U', 2)):
        tot[cv] += 1; bad[cv] += ((-1) ** (r['dims'][cv] - 2) != W[i])
print('   session dims (out_*.jsonl), 937 ratios: mismatches per curve', {c: '%d of %d' % (bad[c], tot[c]) for c in tot})
bad = collections.Counter(); tot = collections.Counter()
for f in sorted(os.listdir(HERE)):
    if not (f.startswith('t10_') and f.endswith('.jsonl')): continue
    for l in open(os.path.join(HERE, f)):
        r = json.loads(l)
        if 'curves' not in r: continue
        W = w(r['m'], r['n'])
        for cv, d in r['curves'].items():
            i = 0 if cv == '08' else 1 if cv == 'LU' else 2
            tot[cv] += 1; bad[cv] += ((-1) ** (d['sel'] - 2) != W[i])
print('   t10 dims (my reruns, t10_*.jsonl): mismatches per curve', {c: '%d of %d' % (bad[c], tot[c]) for c in sorted(tot)})
