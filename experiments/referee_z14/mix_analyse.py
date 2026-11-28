"""z14: summary of mixtest.py's records (mix_*.jsonl), theoremB_z14.py's (thB_*.jsonl) and rank0_z14.py's (rank0_z14.jsonl).
Written by an AI model (reviewer z14); not read by a person.   usage: python3 -I mix_analyse.py | diff - mix_analyse.txt"""
import os, json, glob, collections
HERE = os.path.dirname(os.path.abspath(__file__))
def load(pat):
    R = []
    for f in sorted(glob.glob(os.path.join(HERE, pat))):
        for l in open(f): R.append(json.loads(l))
    return R
kinds = {}
for l in open(os.path.join(HERE, 'ratios.txt')):
    if l.strip() and l[0] != '#':
        w = l.split(); m, n = map(int, w[0].split(',')); kinds[(m, n)] = (w[1], w[2])
print('# z14 summary. Written by an AI model (reviewer z14); not read by a person.')
R = sorted(load('mix_*.jsonl'), key=lambda r: (r['m'], r['n']))
seen = set(); RR = []
for r in R:
    if (r['m'], r['n']) not in seen: seen.add((r['m'], r['n'])); RR.append(r)
R = RR
print('== mixtest: %d curves (ratios outside the pilot, each its own curve up to isomorphism); kinds %s' % (len(R), dict(collections.Counter(kinds[(r['m'], r['n'])] for r in R))))
def tot(field):
    c = collections.Counter()
    for r in R:
        for k, v in r[field].items(): c[k] += v
    return dict(sorted(c.items()))
print('model = t10 Pair (0L, 0U, 8L, 8U) on all curves:', all(all(r['model_ok']) for r in R))
print('dim Sel^2 histogram:', dict(sorted(collections.Counter(r['dim'] for r in R).items())))
print('DFull basis matrix: failures', sum(r['M_fail'] for r in R), '; diagonal all 1:', all(r['M_diag'] for r in R),
      '; symmetry checks', collections.Counter(x for r in R for x in r['M_symcheck']))
print('(a) cassels.py vs DFull (bilinear), pairs (a, b):', tot('cassels'))
print('    a per curve: %d..%d (all of Sel^2 on %d curves)' % (min(r['na'] for r in R), max(r['na'] for r in R), sum(1 for r in R if r['na'] == 2 ** r['dim'] - 1)))
print('    first disagreements/failures:', [(r['m'], r['n'], r['cassels_bad'][:2]) for r in R if r['cassels_bad']][:3])
print('    cassels.py seed 1 vs seed 77 (equal, equal to DFull):', tot('cassels_seed77'))
print('    DFull direct on non-basis a vs bilinear extension:', tot('dfull_bilinear'))
print('(b) Redei form vs DFull, pairs (d, beta) of V = {d: (1,d,d) in Sel^2}:', tot('redei'))
print('    DFull direct on 3 Redei pairs per curve vs bilinear:', tot('redei_direct'))
print('    pairs with d = beta among the sampled Redei pairs: at most %d (all pairs taken on the %d curves with (2^dim V - 1)^2 <= 150)'
      % (sum(2 ** r['dimV'] - 1 for r in R), sum(1 for r in R if (2 ** r['dimV'] - 1) ** 2 <= 150)))
print('    curves by dim V:', dict(sorted(collections.Counter(r['dimV'] for r in R).items())), '; soluble pairs with pairing -1 on', sum(1 for r in R if r['redei'].get('pairing -1')), 'curves')
print('    curves with at least one soluble pair: %d; distinct soluble pairs are distinct by construction (sampled without repetition per curve)' % sum(1 for r in R if r['redei'].get('soluble')))
print('    first disagreements:', [(r['m'], r['n'], r['redei_bad'][:2]) for r in R if r['redei_bad']][:3])
print('(c) V inside N = norms from Q(sqrt W):', collections.Counter(r['V_in_N'] for r in R), '; V supported on -1 and the primes of 2AB:', collections.Counter(r['V_supported_on_2AB'] for r in R))
print('    (dim V, dim N on <-1, primes of 2AB>):', dict(sorted(collections.Counter((r['dimV'], r['dimN_AB']) for r in R).items())))
print('    V = N (on <-1, primes of 2AB>) on', sum(1 for r in R if r['dimV'] == r['dimN_AB']), 'of', len(R), 'curves')
print('    note conic R+ s^2 - R- D^2 = 2AB d w^2: (d in N, qfsolve found a point):', tot('conic_vs_norm'))
print('    odd primes of W with odd exponent, residues mod 24 seen:', sorted(set(x for r in R for x in r['W_odd_primes_mod8'])), '; W with 2 to an odd power:', sum(r['W_has_2_odd'] for r in R))
T = sorted(load('thB_*.jsonl'), key=lambda r: (r['m'], r['n'], r['curve']))
if T:
    print('== Theorem B (EXPLORATION_5 §1) vs DFull: %d curves on %d ratios' % (len(T), len(set((r['m'], r['n']) for r in T))))
    by = collections.defaultdict(collections.Counter)
    for r in T:
        for k, v in r['stats'].items(): by[r['curve']][k] += v
    for cv in sorted(by): print('  %-12s' % cv, dict(sorted(by[cv].items())))
    print('  first disagreements:', [(r['m'], r['n'], r['curve'], r['bad']) for r in T if r['bad']][:3])
    nd = sum(min(2, r['stats'].get('agree', 0) + r['stats'].get('DISAGREE', 0)) for r in T)
    print('  DFull direct calls %d, equal to the bilinear value on %d; tested %d, agree %d, value -1 on %d; skipped (beta not a norm) %d; local factor changed between two points: %d'
          % (nd, sum(r['stats'].get('DFull direct = bilinear', 0) for r in T), sum(r['stats'].get('agree', 0) + r['stats'].get('DISAGREE', 0) for r in T),
             sum(r['stats'].get('agree', 0) for r in T), sum(r['stats'].get('pairing -1', 0) for r in T),
             sum(r['stats'].get('skipped: beta not a norm from Q(sqrt d)', 0) for r in T), sum(r['stats'].get('local factor depends on the point', 0) for r in T)))
    print('  ratios by kind:', dict(collections.Counter(kinds[(m, n)] for m, n in set((r['m'], r['n']) for r in T))))
K = load('rank0_z14.jsonl')
if K:
    print('== rank0 (EXPLORATION_3 §4) recomputed with full_selmer2 + DFull')
    for st in ('nd5', 'sq5', 'nd6'):
        S = [r for r in K if r['set'] == st]
        if not S: continue
        c = collections.Counter()
        for r in S:
            bA, bB = r['EA']['bound'], r['Etw']['bound']
            c['bound E_A = 0'] += bA == 0; c['bound twist = 0'] += bB == 0; c['min = 0'] += min(bA, bB) == 0
            c['parity %s' % str((bA % 2, bB % 2))] += 1
            c['pairing failures'] += r['EA']['fails'] + r['Etw']['fails']; c['diag not 1'] += (not r['EA']['diag_ok']) + (not r['Etw']['diag_ok'])
        print('  %s: %d ratios (%d distinct): %s' % (st, len(S), len(set((r['m'], r['n']) for r in S)), dict(sorted(c.items()))))
    E3 = {(r['set'], r['m'], r['n']): r for r in map(json.loads, open(os.path.join(HERE, 'e3_rank0_perratio.jsonl')))}
    c = collections.Counter(); diffs = []
    for z in K:
        e = E3[(z['set'], z['m'], z['n'])]
        same = (e['dA'], e['bA'], e['dB'], e['bB']) == (z['EA']['dim'], z['EA']['bound'], z['Etw']['dim'], z['Etw']['bound'])
        c[(z['set'], same)] += 1
        if not same: diffs.append(((z['set'], z['m'], z['n']), 'e3', (e['dA'], e['bA'], e['dB'], e['bB']), 'z14', (z['EA']['dim'], z['EA']['bound'], z['Etw']['dim'], z['Etw']['bound'])))
    print('  (dim, bound) on both curves equal to explore3 bound() (e3_rank0_perratio.jsonl):', dict(sorted(c.items())))
    print('  differences:', diffs[:5])
    for st in ('nd5', 'sq5', 'nd6'):
        S = [r for r in E3.values() if r['set'] == st]
        c = collections.Counter()
        for r in S:
            c['bound E_A = 0'] += r['bA'] == 0; c['bound twist = 0'] += r['bB'] == 0; c['min = 0'] += min(r['bA'], r['bB']) == 0
            c['parity %s' % str((r['bA'] % 2, r['bB'] % 2))] += 1
        print('  explore3 bound() on %s (%d ratios): %s' % (st, len(S), dict(sorted(c.items()))))
