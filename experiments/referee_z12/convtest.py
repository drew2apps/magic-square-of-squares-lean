import os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__))  # repository copy: paths relative to this folder
"""Task 2: the six assignments of (b_r, b_m, b_*) to the forms, by z12's own evaluator (zct), against rank_second."""
import json, sys, glob, collections
sys.path.insert(0, HERE)
import zct
out = open(sys.argv[1], 'w'); seeds = [int(x) for x in sys.argv[2].split(',')]
files = sys.argv[3:]
c = collections.Counter(); ex = collections.defaultdict(list)
for f in files:
    for l in open(f):
        r = json.loads(l)
        if 'error' in r: continue
        m, n = r['m'], r['n']
        A, C = (m * m - n * n, m * m + n * n) if r['cv'] == '08' else (8 * m * n, 4 * (m * m + n * n))
        sel = [tuple(b) for b in r['sel']]
        for row in r['rows']:
            if not row['inSel'] or not row['ct']: continue
            kind = 'tuple' if any(w[0] == 'tuple' for w in row['why']) else 'sel'
            dl = row['ap'][1]
            res = [zct.pairing_all(A, C, dl, sel, seed=s)[0] for s in seeds]
            for b, h in zip(sel, row['rs']):
                for pn in zct.PERMS:
                    vs = [rr[(b, pn)] for rr in res]
                    v = vs[0]
                    tag = 'nonconst' if any(isinstance(x, str) for x in vs) else ('seeddep' if len(set(vs)) > 1 else ('agree' if v == h else 'DISAGREE'))
                    c[(r['cv'], pn, tag, 'b_r=1' if b[0] == 1 else 'b_r!=1')] += 1
                    if tag == 'DISAGREE' and len(ex[pn]) < 3: ex[pn].append((m, n, r['cv'], tuple(row['a']), dl, b, h, v))
for k, v in sorted(c.items()): out.write('%s %d\n' % (k, v))
out.write('examples %s\n' % dict(ex)); out.close()
