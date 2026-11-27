"""G4, second pass (EXPERIMENT, written by an AI model; not read by a person; nothing here changes the census).
The first pass (g4_search.py) skips a curve when PARI's ellrank gives rank bounds lo = hi but fewer than lo independent
points.  PARI's lower bound lo can count a 2-Selmer element that no point found explains (see ?? ellrank: r > #L), so
on such a curve the rank is not known to be lo.  Here each such curve is searched on the subgroup that the points
found generate (with the point of t = +-1 when it is independent), exactly as in the first pass: points of height
<= K^2 max h(gens), mapped to t and tested.  This is a partial search: the missing direction is not searched.
usage: g4_pass2.py LOGS... --out FILE [--curves E+,E-]   (reads the first pass's g4_<i>_<N>.json or log_*.jsonl; K = 8; resumes FILE;
       worker G4_I of G4_N from the environment)"""
import sys, os, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import g4_search as S
pari, issq, sqrt_fr, G = S.pari, S.issq, S.sqrt_fr, S.G
K = 8

def todo(logs, only=None):
    out = []
    for f in logs:
        recs = json.load(open(f)) if f.endswith('.json') else [json.loads(l) for l in open(f) if l.startswith('{')]
        for d in recs:
            names = [k for k, v in d['curves'].items() if not v.get('searched') and v['rank'][0] == v['rank'][1] <= S.RMAX
                     and (only is None or k in only)]
            if names: out.append((tuple(d['ratio']), names))
    return sorted(out)

def search(m, n, names):
    rec = dict(ratio=[m, n], curves={}, hits=[], solutions=[])
    for name, E, to_t, Xs in S.curve_list(m, n):
        if name not in names: continue
        kp = []
        for X in Xs:
            y2 = X**3 + S.fr(E[1])*X*X + S.fr(E[3])*X + S.fr(E[4])
            if issq(y2):
                P = pari('[%s, %s]' % (X, sqrt_fr(y2)))
                if pari.ellorder(E, P) == 0: kp.append(P)
        kv = pari('[' + ','.join(str(P) for P in kp) + ']') if kp else pari('[]')
        rk = pari.ellrank(E, 2, kv); lo, hi = int(rk[0]), int(rk[1]); gens = list(rk[3])
        for P in kp:
            cand = gens + [P]; Hm = pari.ellheightmatrix(E, cand); d = float(pari.matdet(Hm)); sc = 1.0
            for i in range(len(cand)): sc *= float(Hm[i, i])
            if d > 1e-9 * sc: gens = cand
        info = dict(rank=[lo, hi], gens=len(gens))
        if not gens:
            info['searched'] = False; rec['curves'][name] = info; continue
        gens = list(pari.ellsaturation(E, gens, 100))
        hmax = max(float(pari.ellheight(E, P)) for P in gens)
        pts, vecs = S.points(E, gens, K*K*hmax)
        info.update(searched='partial' if len(gens) < lo else True, points=len(pts), H=round(K*K*hmax, 3))
        ts = set()
        for P in pts:
            for t in to_t(P): ts.add(t)
        for t in sorted((x for x in ts if x != 'inf'), key=lambda x: (abs(x.numerator) + x.denominator, x)):
            f = S.forms(m, n, t); k = sum(issq(v) for v in f)
            if k >= 4:
                c = G.classify_t(m, n, t)
                rec['hits'].append(dict(curve=name, t=str(t), squares=k, which=[int(issq(v)) for v in f], verdict=c))
                if c == 'SOLUTION': rec['solutions'].append(str(t))
        found = [t for t in (Fr(1), Fr(-1)) if t in ts]
        info['planted'] = dict(carried=len(kp) > 0, found=[str(t) for t in found])
        if kp: assert found, ('planted test fails', m, n, name)
        rec['curves'][name] = info
    return rec

if __name__ == '__main__':
    a = sys.argv[1:]; out = a[a.index('--out') + 1]
    logs = [x for x in a if (x.endswith('.jsonl') or x.endswith('.json')) and x != out]
    only = a[a.index('--curves') + 1].split(',') if '--curves' in a else None
    jobs = todo(logs, only)
    i, N = (int(os.environ.get('G4_I', 0)), int(os.environ.get('G4_N', 1)))
    with open(out, 'a') as fo:
        done = set()
        if os.path.exists(out):
            for l in open(out):
                if l.startswith('{'): done.add(tuple(json.loads(l)['ratio']))
        for (m, n), names in jobs[i::N]:
            if (m, n) in done: continue
            r = search(m, n, names); fo.write(json.dumps(r) + '\n'); fo.flush()
