#!/usr/bin/env python3
"""Which curves close the 19 census pairs of G1_FULL_SYSTEM.md Result 3 that "use L(E,1) != 0 (Kolyvagin) where
2-descent gives only [0, 2]".  Written by an AI model; not read by a person.  PARI-guided (identification, not proof).

The 19 pairs are not recorded anywhere with their curve.  They are recovered from the history of open_pairs_1000.txt:
(112,15) and (127,97) were closed by L-values already at 4ad9ae6 (not in that file's 174 open pairs, and named in
G1_FULL_SYSTEM.md at 4ad9ae6); 17 more left the file at d185af6 ("17 further closures by L-values"): 16 with a trivial
point and (700,529).  For each pair this script runs global.py's twelve curves (global.torsion_w: ellrank on the curve
and on its isogeny class, then the torsion test of global.analyse) and, for every curve with rank bounds [0, hi > 0],
the conductor, the root number and, when w = +1 and N <= 10^14, PARI's lfun(E, 1).  A pair is "closed by L" on a curve
when lfun(E,1) is clearly nonzero there and no torsion point gives a SOLUTION (global.classify_t).

usage: census_lvalue_identify.py [--write]      -> census_lvalue_identify.txt, census_lvalue_identify.json"""
import sys, os, json, time, importlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
G = importlib.import_module('global')            # scripts/g1_full/global.py (its name is a Python keyword)
pari = G.pari

OLD2 = [(112, 15), (127, 97)]                     # closed by L-values at 4ad9ae6
NEW17 = [(240, 161), (260, 69), (288, 175), (329, 191), (352, 135), (401, 79), (463, 113), (487, 217), (527, 336),
         (544, 33), (575, 48), (577, 511), (623, 527), (697, 696), (700, 529), (863, 191), (943, 576)]   # left the list at d185af6
PAIRS = OLD2 + NEW17
NMAX = 10**14

def model(poly):
    """global.torsion_w's Weierstrass model [0, c2, 0, c1 c3, c0 c3^2] of y^2 = poly(w)"""
    if len(poly) == 5: c = [poly[4], poly[3], poly[2], poly[1]]
    else: c = poly
    c0, c1, c2, c3 = c
    return [0, c2, 0, c1 * c3, c0 * c3 * c3]

def one(m, n, lines):
    rec = dict(pair=[m, n], curves=[])
    for name, poly in G.curves(m, n).items():
        lo, hi, ws, T = G.torsion_w(poly)
        verdicts = set()
        for w in ws:
            if w is None: verdicts.add('degenerate'); continue
            if w < 0 or not G.issq(w): continue
            t = G.sqrt_fr(w)
            for tt in {t, -t}: verdicts.add(G.classify_t(m, n, tt))
        c = dict(curve=name, rank=[lo, hi], torsion=T, solution=('SOLUTION' in verdicts), model=model(poly))
        if lo == 0 and hi > 0:
            E = pari.ellminimalmodel(pari.ellinit(c['model']))
            N = int(pari.ellglobalred(E)[0]); w = int(pari.ellrootno(E))
            c.update(N=N, w=w)
            if w == 1 and N <= NMAX:
                t0 = time.time(); L = float(pari.lfun(pari.lfuncreate(E), 1)); c.update(L=L, lfun_s=round(time.time() - t0, 1))
        rec['curves'].append(c)
    closers = [c for c in rec['curves'] if c['rank'][1] == 0 and not c['solution']]
    lcl = [c for c in rec['curves'] if c.get('L') is not None and abs(c['L']) > 1e-3 and not c['solution']]
    rec['closed_by_rank'] = [c['curve'] for c in closers]
    rec['closed_by_L'] = sorted(((c['N'], c['curve'], c['L']) for c in lcl))
    und = [(c['curve'], c['rank'], c.get('N'), c.get('w'), c.get('L')) for c in rec['curves'] if c['rank'][0] == 0 and c['rank'][1] > 0]
    lines.append('%d %d  rank-0 by ellrank: %s;  L-candidates [0,hi]: %s' % (m, n, rec['closed_by_rank'] or 'none',
                 '; '.join('%s %s N=%s w=%s L=%s' % (a, tuple(b), cN, cw, ('%.6f' % cL) if cL is not None else '-') for a, b, cN, cw, cL in und) or 'none'))
    lines.append('        closed by L (N <= 1e14, w = +1, lfun > 1e-3, no torsion solution): %s' % (
        ', '.join('%s N=%d L=%.6f' % (nm, N, L) for N, nm, L in rec['closed_by_L']) or 'NONE'))
    print('\n'.join(lines[-2:]), flush=True)
    return rec

if __name__ == '__main__':
    lines = ['census_lvalue_identify.py: the 19 L-value closures of G1_FULL_SYSTEM.md Result 3 (PARI %s; identification only)' % '.'.join(str(x) for x in pari.version()), '']
    out = []; t0 = time.time()
    sel = PAIRS
    if '--only' in sys.argv: sel = [PAIRS[int(i)] for i in sys.argv[sys.argv.index('--only') + 1].split(',')]
    for (m, n) in sel: out.append(one(m, n, lines))
    nL = sum(bool(r['closed_by_L']) for r in out); nR = sum(bool(r['closed_by_rank']) for r in out)
    lines.append(''); lines.append('pairs %d; closed by a rank-0 curve (ellrank): %d; closed by an L-value curve: %d; time %.0f s' % (len(out), nR, nL, time.time() - t0))
    print(lines[-1])
    if '--write' in sys.argv:
        open(os.path.join(HERE, 'census_lvalue_identify.txt'), 'w').write('\n'.join(lines) + '\n')
        json.dump(out, open(os.path.join(HERE, 'census_lvalue_identify.json'), 'w'), indent=0)
