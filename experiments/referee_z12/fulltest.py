import os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__))  # repository copy: paths relative to this folder
"""z12 task 1: full pairing matrix by rank_second (DFull.pairing) against ct08.pairing, on E_08 and E_LU.
a ranges over: the class of every residual tuple (both signs), and every element of Sel^2 with trivial e_r component;
b ranges over the basis of Sel^2 returned by full_selmer2 (and, for ct08, also over all of Sel^2 for linearity).
Results appended to fulltest.jsonl, one line per (ratio, curve)."""
import sys, os, json, time, random
sys.argv = [sys.argv[0]] + [a for a in sys.argv[1:]]
ARGS = sys.argv[1:]; sys.argv = sys.argv[:1]
T10 = _os.path.join(HERE, '..', 'conductor-t10')
sys.path.insert(0, T10)
import t10_pairct as T
RS = T.RS; RI = T.RI
sys.path.insert(0, _os.path.join(HERE, '..', 'research-pairing'))
os.chdir(_os.path.join(HERE, '..', 'conductor-t10'))
import ct08
from fractions import Fraction as Fr
sys.path.insert(0, HERE)
import zparse

def sqf(x):
    x = Fr(x); return RI.sqfree_part(x.numerator * x.denominator)
def span_all(sel):
    els = [(1, 1)]
    for b in sel:
        els = els + [(sqf(x[0] * b[0]), sqf(x[1] * b[1])) for x in els]
    return els      # element i = product of sel[k] over the bits k of i

def run(m, n, cv, out, seeds=(1,)):
    t0 = time.time()
    j, k = cv[0], cv[1]
    Pd = T.PairData(m, n, j, k)
    es = Pd.Pr.es
    assert Pd.role['r'] == 0 and Pd.role['*'] == 2, Pd.role
    if cv == '08': A, C = m * m - n * n, m * m + n * n
    else: A, C = 8 * m * n, 4 * (m * m + n * n)
    assert (es[1] - es[0], es[2] - es[0]) == (A * A, C * C), ('model', es, A, C)
    sel = Pd.sel
    allsel = span_all(sel)
    assert len(allsel) == 2 ** len(sel) and len(set(allsel)) == len(allsel)
    nall = 0
    R = zparse.load()[(m, n)]
    tuples = [t['tup'] for t in R['tuples']]
    aset = {}
    for t in tuples:
        for sg in (1, -1):
            c = Pd.cls_of_tuple(t, sg)
            aset.setdefault(c, []).append(('tuple', t, sg))
    for c in allsel:
        if c[0] == 1: aset.setdefault(c, []).append(('sel(1,d,d)',))
    # B = C^2 - A^2 = B0^2
    from math import isqrt
    B0 = isqrt(C * C - A * A); assert B0 * B0 == C * C - A * A
    img4 = (sqf(C * (C + B0)), sqf(B0 * (B0 + C)))   # image of a point of order 4 (x = C^2 + C B0)
    rec = dict(m=m, n=n, cv=cv, dim=len(sel), sel=sel, rows=[])
    for a, why in aset.items():
        inSel = Pd.vec(a) is not None and len(RS.f2_span(Pd.selv + [Pd.vec(a)])) == len(sel)
        row = dict(a=a, why=[list(map(str, w)) for w in why], inSel=inSel, rs=[], ct={}, ct_all={})
        if not inSel: rec['rows'].append(row); continue
        D = T.DFullX(Pd.esh, a)
        for b in sel:
            h, per = T.safe_pair(D, b)
            row['rs'].append(h)
        # ct08 needs a = (1, delta, delta): move a by the image of a point of order 4 if needed
        if a[0] == 1: ap = a; row['shift'] = None
        else:
            ap = (sqf(a[0] * img4[0]), sqf(a[1] * img4[1])); row['shift'] = img4
        row['ap'] = ap
        if ap[0] != 1: row['ct_note'] = 'not of the form (1,d,d) even after the order-4 shift'
        else:
            dl = ap[1]
            for sd in seeds:
                vals = []
                for b in sel:
                    v, loc = ct08.pairing(A, C, dl, b[0], sqf(b[0] * b[1]), seed=sd)
                    vals.append(v)
                row['ct'][sd] = vals
            # linearity over all of Sel^2 (seed 1 only; tuple classes and at most 8 other classes)
            if why[0][0] != 'tuple':
                nall += 1
                if nall > 8: rec['rows'].append(row); continue
            va = []
            for b in allsel:
                v, loc = ct08.pairing(A, C, dl, b[0], sqf(b[0] * b[1]), seed=1)
                va.append(v)
            row['ct_all'] = dict(b=allsel, v=va)
        rec['rows'].append(row)
    rec['time'] = time.time() - t0
    with open(out, 'a') as f: f.write(json.dumps(rec) + '\n')
    return rec

if __name__ == '__main__':
    out = ARGS[0]; seeds = tuple(int(x) for x in ARGS[1].split(','))
    for a in ARGS[2:]:
        m, n, cv = a.split(',')
        try:
            r = run(int(m), int(n), cv, out, seeds)
            print(m, n, cv, 'dim', r['dim'], 'rows', len(r['rows']), 'time %.1f' % r['time'], flush=True)
        except Exception as ex:
            import traceback; traceback.print_exc()
            with open(out, 'a') as f: f.write(json.dumps(dict(m=int(m), n=int(n), cv=cv, error=repr(ex))) + '\n')
