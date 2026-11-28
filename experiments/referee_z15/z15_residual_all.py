"""z15: residual tuples of every ratio of the two G2 lists that keep a tuple, with the repository's code only
(tuple_descent.magic, exact_local.allowed_at at every bad prime, exact_local.join, as t10_pairct.residual does, without
the good-prime search), compared with the session's out_*.jsonl.  Written by an AI model (reviewer z15), not read by
a person.  usage: python3 z15_residual_all.py > z15_residual_all.txt   (writes z15_residual_all.jsonl)"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
os.environ.setdefault('MSS_REPO', REPO)
sys.path.insert(0, os.path.join(REPO, 'experiments', 'roadmap-g', 'g1'))
import exact_local as X
TD = X.TD
import cypari2
PARI = cypari2.Pari(); PARI.allocatemem(5 * 10 ** 8, 2 * 10 ** 9)
def fact(x):
    x = abs(int(x))
    if x <= 1: return {}
    f = PARI.factor(x); out = {int(q): int(e) for q, e in zip(f[0], f[1])}
    for q in out: assert PARI.isprime(q)
    return out
TD.factorint = fact
def residual(m, n):
    forms, sup, known = TD.magic(m, n)
    divs = []
    for S in sup:
        ds = [1]
        for q in S: ds += [x * q for x in ds]
        divs.append(ds)
    bad = set()
    for x in (2 * m * n * (m * m - n * n), m * m + n * n, m ** 8 + 14 * m ** 4 * n ** 4 + n ** 8): bad |= set(fact(x))
    loc = []
    for p in sorted(bad):
        L, al = X.allowed_at(forms, p)
        if len(al) < len(L.labels()) ** len(forms): loc.append((p, L, al))
    return sorted(X.join(divs, loc)), known
def read_list(path):
    out = []
    for l in open(path):
        l = l.strip()
        if not l or l[0] == '#' or 'left' not in l: continue
        out.append(tuple(map(int, l.split()[0].split(','))))
    return out
G2 = os.path.join(REPO, 'experiments', 'roadmap-g', 'g2')
L = read_list(os.path.join(G2, 'g2_997_all.txt')) + read_list(os.path.join(G2, 'g2_15000_need_good.txt'))
ses = {}
for f in ('out_997.jsonl', 'out_15000_0.jsonl', 'out_15000_1.jsonl'):
    for l in open(os.path.join(HERE, 'session_data', f)):
        r = json.loads(l); ses[(r['m'], r['n'])] = sorted(tuple(t) for t, d in r['tuples'])
print('z15_residual_all.py: written by an AI model (reviewer z15), not read by a person.')
print('ratios in the G2 lists with tuples left:', len(L), '; ratios in the session output:', len(ses), '; same set:', set(L) == set(ses))
nt = 0; diff = []; known_n = 0; t0 = time.process_time()
with open(os.path.join(HERE, 'z15_residual_all.jsonl'), 'w') as fo:
    for m, n in sorted(L):
        T, known = residual(m, n)
        nt += len(T); known_n += len(known)
        if T != ses.get((m, n)): diff.append((m, n))
        fo.write(json.dumps(dict(m=m, n=n, tuples=[list(t) for t in T], known=[list(t) for t in known])) + '\n')
print('residual tuples (repository code):', nt, '; known tuples:', known_n, '; ratios where the session differs:', len(diff), diff[:10])
sys.stderr.write('cpu %.0f s\n' % (time.process_time() - t0))
