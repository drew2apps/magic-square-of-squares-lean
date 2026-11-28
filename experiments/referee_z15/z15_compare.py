"""z15: tuple-by-tuple, curve-by-curve comparison of the session's out_*.jsonl (cassels.py/sel.py pairing) with
(A) the run-9 pilot's own output files (experiments/conductor-t10/t10_*.txt, parsed here with my own parser) and
(B) my reruns of the repository's t10_pairct (t10_*.jsonl in this directory: five survivors, sample80, ghost ratios).
Decision per (tuple, curve): killed (t10 'killed' or 'outside Sel2'; session 1 or 'outside') or not (t10 'image' or
'kernel'; session 0).  Also residual sets, Sel^2 dimensions per curve, least killing good primes.
Written by an AI model (reviewer z15), not read by a person.  usage: python3 z15_compare.py > z15_compare.txt"""
import os, re, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
CV = ['08', 'LU', '0L', '0U', '8L', '8U']
ses = {}
for f in ('out_997.jsonl', 'out_15000_0.jsonl', 'out_15000_1.jsonl'):
    for l in open(os.path.join(HERE, 'session_data', f)):
        r = json.loads(l); ses[(r['m'], r['n'])] = r
def sdec(v): return 'K' if v in (1, 'outside') else 'N' if v == 0 else '?'
def tdec(v): return 'K' if v in ('killed', 'outside Sel2', 'outside') else 'N' if v in ('image', 'kernel') else '?'
print('z15_compare.py: written by an AI model (reviewer z15), not read by a person.')
# (A) the pilot files
pil = collections.defaultdict(dict); ratio = None
for f in ('t10_997_a.txt', 't10_997_b.txt', 't10_15000_s100_a.txt', 't10_15000_s100_b.txt'):
    for l in open(os.path.join(HERE, '..', 'conductor-t10', f)):
        mm = re.match(r'ratio \((\d+),(\d+)\)', l)
        if mm: ratio = (int(mm.group(1)), int(mm.group(2))); continue
        mm = re.match(r'    \(([0-9, ]+)\)  good prime (\S+)  (.*)  => ', l)
        if mm:
            t = tuple(int(x) for x in mm.group(1).split(','))
            pil[ratio][t] = {x.split(':')[0]: x.split(':', 1)[1] for x in mm.group(3).split()}
agree = collections.Counter(); setdiff = []; dis = []
for r, d in sorted(pil.items()):
    st = {tuple(t): v for t, v in ses[r]['tuples']}
    if set(st) != set(d): setdiff.append(r)
    for t in d:
        for c in CV:
            a, b = sdec(st[t][c]) if t in st else 'missing', tdec(d[t][c])
            agree[(a, b)] += 1
            if a != b: dis.append((r, t, c, a, b))
print('(A) pilot files: %d ratios, %d tuples, %d (tuple, curve) decisions; (session, pilot) counts %s; residual sets differ at %d ratios; disagreements %d %s'
      % (len(pil), sum(len(d) for d in pil.values()), sum(agree.values()), dict(agree), len(setdiff), len(dis), dis[:5]))
# (B) my reruns
runs = {}; dup = []
for f in sorted(os.listdir(HERE)):
    if f.startswith('t10_') and f.endswith('.jsonl'):
        for l in open(os.path.join(HERE, f)):
            r = json.loads(l); k = (r['m'], r['n'])
            if k in runs:     # a ratio run twice (two of my lists overlap): the two runs must agree exactly
                a, b = runs[k][1], r
                dup.append((k, runs[k][0], f, a.get('tuples') == b.get('tuples') and a.get('curves') == b.get('curves')))
                continue
            runs[k] = (f, r)
agree = collections.Counter(); dis = []; setdiff = []; dimdiff = []; errs = []; cpu = []; inpilot = 0
undec = []; byfile = collections.Counter(); distinct = collections.Counter()
for k, (f, r) in sorted(runs.items()):
    if 'error' in r: errs.append((k, r['error'])); continue
    inpilot += k in pil
    byfile[f] += 1
    st = {tuple(t): v for t, v in ses[k]['tuples']}
    mt = {tuple(t): d for t, g, d in r['tuples']}
    if set(st) != set(mt): setdiff.append(k)
    for c in CV:
        if c in r['curves'] and r['curves'][c]['sel'] != ses[k]['dims'][c]: dimdiff.append((k, c, r['curves'][c]['sel'], ses[k]['dims'][c]))
    for t in mt:
        for c in CV:
            a, b = (sdec(st[t][c]) if t in st else 'missing'), tdec(mt[t][c])
            if b == '?': undec.append((k, t, c, mt[t][c]))
            agree[(a, b)] += 1
            if a != b: dis.append((k, t, c, st.get(t, {}).get(c), mt[t][c]))
print('(B) my t10 reruns: %d distinct ratios (%s), of them in the pilot: %d; errors: %s' % (sum(byfile.values()), dict(byfile), inpilot, errs))
print('    ratios run twice (identical output?): %s' % dup)
print('    tuples %d, (tuple, curve) decisions %d; (session, t10) counts %s' % (sum(len(r['tuples']) for f, r in runs.values() if 'tuples' in r), sum(agree.values()), dict(agree)))
print('    residual sets differ at %d ratios %s; Sel2 dims differ: %d %s; undecided t10 verdicts: %d %s' % (len(setdiff), setdiff[:5], len(dimdiff), dimdiff[:5], len(undec), undec[:5]))
print('    DISAGREEMENTS: %d' % len(dis))
for x in dis: print('      ', x)
surv_t10 = sorted((k, tuple(t)) for k, (f, r) in runs.items() if 'tuples' in r for t, g, d in r['tuples'] if all(tdec(v) == 'N' for v in d.values()))
print('    tuples surviving all six curves in my reruns:', surv_t10)
gp = collections.Counter(g for k, (f, r) in runs.items() if 'tuples' in r for t, g, d in r['tuples'])
print('    least killing good prime of the residual tuples (t10, <= 5000):', sorted(gp.items(), key=lambda x: (x[0] is None, x[0] or 0)))
# distinct decisions: one per (ratio, curve, class), the class by the session's formulas (pipeline.py: 08 (1, mn d0 d8),
# LU (2C, 2 dL dU), mixed (2, k d_j d_k) with k = sqf(mn(m^2-n^2)), 2k for 0U and 8L); the four mixed curves are one
# curve in the session's model, so also counted with the four merged
import cypari2
_P = cypari2.Pari()
def _sqf(x): return int(_P.core(abs(int(x))))
def cls(m, n, c, t):
    A, B, C = m * m - n * n, 2 * m * n, m * m + n * n; k = _sqf(m * n * A)
    return {'08': (1, _sqf(m * n * t[0] * t[1])), 'LU': (_sqf(2 * C), _sqf(2 * t[2] * t[3])), '0L': (2, _sqf(k * t[0] * t[2])),
            '8U': (2, _sqf(k * t[1] * t[3])), '0U': (2, _sqf(2 * k * t[0] * t[3])), '8L': (2, _sqf(2 * k * t[1] * t[2]))}[c]
for label, merge in (('per curve', False), ('mixed curves merged', True)):
    dd = {}; clash = 0
    for k, (f, r) in runs.items():
        if 'tuples' not in r: continue
        st = {tuple(t): v for t, v in ses[k]['tuples']}
        for t, g, d in r['tuples']:
            for c in CV:
                key = (k, 'mix' if (merge and c not in ('08', 'LU')) else c, cls(k[0], k[1], c, t))
                val = (sdec(st[tuple(t)][c]), tdec(d[c]))
                if key in dd and dd[key] != val: clash += 1
                dd[key] = val
    print('    distinct (ratio, curve, class) decisions, %s: %d; (session, t10) %s; same class with two verdicts: %d'
          % (label, len(dd), dict(collections.Counter(dd.values())), clash))
# timings: t10_times.txt (moved out of the run records by z15_strip.py)
