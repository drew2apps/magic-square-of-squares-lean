"""Reviewer z13: the counts of the independent test (indep_a.jsonl, indep_b.jsonl from z13_indep.py), and four
side checks of the proofs that use only the stored Sel^2 bases and z13_formula:
  (e1) at p = 1 mod 4 dividing AC, where d and -d are both squares, the two branches (i), (ii) of T_p agree;
  (e2) the sign of sqrt(d) does not matter: (-1, b_r)_p = 1 at every p | AC (and 2 | A) for b in Sel^2;
  (e3) on E_08 at a ratio of square type, (1, 2mn, 2mn) is in S' and in the radical of Theta;
  (e4) the representative of d matters: with d replaced by q^2 d (q a good prime) formula (2) changes on some b;
  (e5) on E_08 at the other ratios: is (1, 2mn, 2mn) in S' (a Legendre condition) and in the radical of Theta.
Written by an AI model (reviewer z13), not read by a person.  Deterministic; prints a diffable summary."""
import os, sys, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z13_formula as F

recs = []
for fn in ('indep_a.jsonl', 'indep_b.jsonl'):
    p = os.path.join(HERE, fn)
    if os.path.exists(p): recs += [json.loads(l) for l in open(p)]
recs.sort(key=lambda d: (d['cat'], d['m'], d['n'], d['cv']))
err = [d for d in recs if 'error' in d]
ok = [d for d in recs if 'error' not in d]
print('records', len(recs), 'errors', len(err))
for d in err: print('  ERROR', d['m'], d['n'], d['cv'], d['error'][:120])

def span(basis):
    els = [(1, 1)]
    for b in basis:
        els = els + [(F.sqf(x[0] * b[0]), F.sqf(x[1] * b[1])) for x in els]
    return els

tot = collections.Counter(); bycat = collections.defaultdict(collections.Counter)
curves = set(); exS = []; exF = []; exT = []
for d in ok:
    c = bycat[d['cat']]
    key = (d['A'], d['C'])
    c['curves'] += 1; tot['curves'] += 1; curves.add(key)
    c['H2 holds'] += d['H2']; c['A odd'] += d['A'] % 2
    for k in ('agree', 'disagree', 'undefined', 'agree_direct', 'disagree_direct'):
        c['pairs ' + k] += d['cmp'][k]; tot['pairs ' + k] += d['cmp'][k]
    c['theta agree'] += d['theta']['agree']; c['theta disagree'] += d['theta']['disagree']
    tot['theta agree'] += d['theta']['agree']; tot['theta disagree'] += d['theta']['disagree']
    c["S' = {b_r = 1}"] += d['Sprime_eq']; tot["S' = {b_r = 1}"] += d['Sprime_eq']
    c['a evaluated by rank_second'] += d['n_a_direct']; c['a in total'] += d['n_a']
    c['bilinearity checks'] += d['lin_n']; c['bilinearity failures'] += d['lin_bad']
    tot['bilinearity checks'] += d['lin_n']; tot['bilinearity failures'] += d['lin_bad']
    c['rank_second failures'] += len(d['rs_fail'])
    if not d['Sprime_eq'] and len(exS) < 3: exS.append(d)
    if d['cmp']['disagree'] and len(exF) < 3: exF.append(d)
    if d['theta']['disagree'] and len(exT) < 3: exT.append(d)
    if d['A'] % 2 == 1:
        v = min(F.vp(d['C'] - d['A'], 2), 99), min(F.vp(d['C'] + d['A'], 2), 99)
        c['A odd, max v2(C-+A) = %d' % max(v)] += 1
print('curves', tot['curves'], '; distinct curves (A, C)', len(curves))
print('dims of Sel^2:', dict(sorted(collections.Counter(d['dim'] for d in ok).items())))
print("dims of S':", dict(sorted(collections.Counter(d['dimS1'] for d in ok).items())))
for k in sorted(tot): print('total', k, tot[k])
# within the note's hypotheses, on distinct curves ((m, n) -> (m + n, m - n) exchanges E_08 and E_LU, so a curve can
# occur twice), and the smallest disagreement outside H2
seenAC = set(); dist = collections.Counter()
for d in sorted(ok, key=lambda d: (d['C'], d['A'])):
    if not d['H2'] or (d['A'], d['C']) in seenAC: continue
    seenAC.add((d['A'], d['C'])); dist['curves'] += 1
    for k in ('agree', 'disagree', 'undefined', 'agree_direct', 'disagree_direct'): dist['pairs ' + k] += d['cmp'][k]
    dist['theta agree'] += d['theta']['agree']; dist['theta disagree'] += d['theta']['disagree']; dist["S' = {b_r = 1}"] += d['Sprime_eq']
    dist['A odd'] += d['A'] % 2
    dist['A odd, v2(C -+ A) = 5'] += d['A'] % 2 == 1 and max(F.vp(d['C'] - d['A'], 2), F.vp(d['C'] + d['A'], 2)) == 5
    dist['3 | B'] += d['B'] % 3 == 0
print('within H2, distinct curves:', dict(sorted(dist.items())))
bad = sorted((d for d in ok if d['cmp']['disagree']), key=lambda d: (d['C'], d['A']))
if bad:
    d = bad[0]
    print('smallest disagreement (all outside H2):', d['m'], d['n'], d['cv'], 'A, B, C =', d['A'], d['B'], d['C'], 'v2(C - A), v2(C + A) =', F.vp(d['C'] - d['A'], 2), F.vp(d['C'] + d['A'], 2), 'Sel2 basis', d['sel'], d['ex'][:1])
print('curves with a disagreement:', len(bad), 'of which within H2:', sum(d['H2'] for d in bad))
for cat in sorted(bycat):
    print(cat, dict(sorted(bycat[cat].items())))
for d in exS: print("S' differs:", d['m'], d['n'], d['cv'], 'H2', d['H2'], 'Sel2 part', d['S1'], "S'", d['Sprime'])
for d in exF: print('formula (2) disagrees:', d['m'], d['n'], d['cv'], 'H2', d['H2'], d['cmp'], d['ex'][:2])
for d in exT: print('Theta disagrees:', d['m'], d['n'], d['cv'], 'H2', d['H2'], d['theta'], d['theta_ex'][:2])

# side checks on the curves within the note's hypotheses
e1 = collections.Counter(); e2 = collections.Counter(); e3 = collections.Counter(); e4 = []; e5 = collections.Counter()
for d in ok:
    if not d['H2']: continue
    A, B, C = d['A'], d['B'], d['C']
    allsel = span([tuple(x) for x in d['sel']])
    S1 = [c[1] for c in allsel if c[0] == 1]
    places = set(F.factor(A * C)) | ({2} if A % 2 == 0 else set())
    for dl in S1:
        for p in sorted(places):
            sq, msq = F.is_sq_Qp(dl, p), F.is_sq_Qp(-dl, p)
            for b in allsel:
                b_r, b_m = b; b_s = F.sqf(b_r * b_m)
                if sq or msq:
                    e2['(-1, b_r)_p = 1' if F.hilb(-1, b_r, p) == 1 else '(-1, b_r)_p = -1'] += 1
                if sq and msq:
                    t1 = F.hilb(2 * F.padic_sqrt(dl, p), b_r, p)
                    t2 = F.hilb(-A, b_s, p) * F.hilb(-C, b_m, p) * F.hilb((A - C) * F.padic_sqrt(-dl, p), b_r, p)
                    e1['agree' if t1 == t2 else 'DIFFER'] += 1
    if d['cat'] == 'square' and d['cv'] == '08':
        m, n = d['m'], d['n']
        c = F.sqf(2 * m * n)
        inS = c in S1
        rad = inS and all(F.Theta(A, B, C, c, be) == 1 for be in S1)
        e3['(1, 2mn, 2mn) in S\' %s, in the radical of Theta %s, 2mn = +-square %s' % (inS, rad, abs(c) == 1)] += 1
    if d['cv'] == '08' and d['cat'] not in ('square', 'twicesquare'):
        m, n = d['m'], d['n']
        c = F.sqf(2 * m * n)
        inS = c in S1
        rad = inS and all(F.Theta(A, B, C, c, be) == 1 for be in S1)
        e5['m, n of %s parity: (1, 2mn, 2mn) in S\' %s, in the radical %s' % ('opposite' if (m + n) % 2 else 'odd', inS, rad)] += 1
    if len(e4) < 1 and d['cv'] == '08':
        # a b with b_r of odd valuation at some p | AC, and a good prime q that is a non-residue mod that p
        for b in allsel:
            pr = [p for p in places if p != 2 and F.vp(abs(b[0]), p) % 2 == 1]
            if not pr: continue
            p = pr[0]
            q = next(q for q in range(3, 10 ** 5, 2) if F.pari.isprime(q) and (2 * A * B * C) % q and pow(q, (p - 1) // 2, p) == p - 1)
            for dl in S1:
                if not (F.is_sq_Qp(dl, p) or F.is_sq_Qp(-dl, p)): continue
                f1 = F.formula2(A, B, C, dl, b[0], b[1])
                # the same formula with the representative q^2 d (no reduction to the squarefree part)
                val = 1
                for pp in F.factor(B):
                    if pp == 2:
                        if (C + A) % 8 == 0: val *= F.hilb(dl, F.sqf(b[0] * b[1]), 2)
                    elif (C + A) % pp == 0: val *= F.hilb(dl, F.sqf(b[0] * b[1]), pp)
                for pp in sorted(places):
                    D = q * q * dl
                    if F.is_sq_Qp(D, pp): val *= F.hilb(2 * F.padic_sqrt(D, pp), b[0], pp)
                    else: val *= F.hilb(-A, F.sqf(b[0] * b[1]), pp) * F.hilb(-C, b[1], pp) * F.hilb((A - C) * F.padic_sqrt(-D, pp), b[0], pp)
                if val != f1:
                    e4.append((d['m'], d['n'], d['cv'], 'd', dl, 'b', b, 'p', p, 'q', q, 'squarefree d', f1, 'q^2 d', val)); break
            if e4: break
print('(e1) branches (i), (ii) where both apply:', dict(e1))
print('(e2) sign of sqrt(d):', dict(e2))
print('(e3) square type, E_08:', dict(e3))
print('(e4) representative q^2 d instead of d:', e4[:1])
print('(e5) E_08, not of square type:', dict(sorted(e5.items())))
