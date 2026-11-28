import os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__))  # repository copy: paths relative to this folder
"""z12's own parser of the t10 outputs: per ratio, the header of every pair curve, all residual tuples (per-tuple
lines), and the kill lines."""
import re, glob, os
T10 = _os.path.join(HERE, '..', 'conductor-t10')
def load():
    out = {}
    for f in sorted(glob.glob(os.path.join(T10, 't10_*_[ab].txt'))):
        r = None; cv = None
        for l in open(f):
            mm = re.match(r'ratio \((\d+),(\d+)\): (\d+) residual', l)
            if mm:
                r = (int(mm[1]), int(mm[2])); out[r] = dict(file=os.path.basename(f), curves={}, tuples=[], nres=int(mm[3])); cv = None; continue
            mm = re.match(r'  E_(\w\w): base t0 = (\S+); e = \(([-\d, ]+)\); e_r = e(\d), e\* = e(\d).*dim Sel2 = (\d+), point images (\d+)', l)
            if mm:
                cv = mm[1]; out[r]['curves'][cv] = dict(e=tuple(int(x) for x in mm[3].split(',')), er=int(mm[4]), es=int(mm[5]), dim=int(mm[6]), pim=int(mm[7]), kills=[]); continue
            mm = re.match(r'    tuple \(([-\d, ]+)\): class \(([-\d]+), ([-\d]+)\): (\S+) ?(.*)', l)
            if mm and cv:
                tup = tuple(int(x) for x in mm[1].split(','))
                k = re.match(r'b = \(([-\d]+), ([-\d]+)\), local -1 at \[(.*)\]', mm[5])
                out[r]['curves'][cv]['kills'].append(dict(tup=tup, cls=(int(mm[2]), int(mm[3])), verdict=mm[4], b=(int(k[1]), int(k[2])) if k else None, places=k[3] if k else None))
                continue
            mm = re.match(r'    \(([-\d, ]+)\)  good prime (\S+)  (.*)  => (.*)$', l.rstrip('\n'))
            if mm:
                tup = tuple(int(x) for x in mm[1].split(','))
                ver = dict(x.split(':', 1) for x in mm[3].split())
                out[r]['tuples'].append(dict(tup=tup, good=mm[2], ver=ver))
    return out
if __name__ == '__main__':
    R = load()
    print('ratios', len(R), 'residual tuples', sum(len(v['tuples']) for v in R.values()), 'nres check', all(len(v['tuples']) == v['nres'] for v in R.values()))
    for cv in ('08', 'LU'):
        kl = [(r, k['cls'], k['b']) for r, v in R.items() for k in v['curves'].get(cv, {}).get('kills', [])]
        print(cv, 'kill lines', len(kl), 'distinct (ratio,class,b)', len(set(kl)), 'distinct (ratio,class)', len(set((a, b) for a, b, c in kl)),
              'ratios with kill lines', len(set(a for a, b, c in kl)), 'verdicts', set(k['verdict'] for v in R.values() for k in v['curves'].get(cv, {}).get('kills', [])))
    allk = sum(len(c['kills']) for v in R.values() for c in v['curves'].values())
    print('all kill lines', allk)
