"""I8: the Cassels-Tate pairing test of conductor-t10 on every ratio that needs a good prime (EXPERIMENT, written by an
AI model; not read by a person; nothing here changes the census; tier D).
For each ratio of i2_needgood_list.txt (G2's 937, H2's 2,807 and the 132 of H2's hunts only), each residual tuple tau
(the tuples allowed at every bad prime) is tested on the six pair curves C_jk: y^2 = F_j F_k by t10_pairct.py
(experiments/conductor-t10; its PairData, unchanged).  The residual tuples and their killing good primes are read from
I2's residual file (i2_residual.jsonl, or the .part file of the run while it is running), each tuple checked here to be
allowed at every bad prime and not allowed at its killing prime; for a ratio not yet in that file t10's own residual()
computes them (source "t10").
Per tuple:  "pairing" (killed on some pair curve: outside Sel^2, or <b, c(tau)> = -1; the curves are listed),
else "rank 0 of E_jk" (the ratio is excluded by a pair curve of rank 0: dim Sel^2 = 2 with full rational 2-torsion,
the torsion group complete (as many points as the gcd of #E(F_p) over 25 good odd primes, counted exactly) and every
point over t in {0, 1, -1, infinity}), else "neither" (or "undecided" when a local point was not found).
Output: one JSON object per line, in the order of the list:
  {"m", "n", "source", "tuples": [{"tau", "kill", "pairs": {"08": ..., ...}, "verdict"}], "rank0": {jk: {...}}}
usage: i8_pairct.py i N [--max K] [--write]    worker i of N (i8_<i>of<N>.jsonl; .part while running); at most K new
                                                ratios per process (exit 3 if more remain: run it again)
       i8_pairct.py --merge N [--write]         -> i8_pairct.jsonl, then the summary -> i8_summary.txt
       i8_pairct.py --compare50 [--write]       50 ratios (seed 1) recomputed by t10's residual() -> i8_compare50.txt
       i8_pairct.py --negative [--write]        -> i8_negative.txt
Folders: MSS_REPO (default ../../..), T10_DIR (default MSS_REPO/experiments/conductor-t10), I2_DIR (default ../i2)."""
import sys, os, json, random
from fractions import Fraction as Fr
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..', '..'))
T10 = os.environ.get('T10_DIR', os.path.join(REPO, 'experiments', 'conductor-t10'))
I2 = os.environ.get('I2_DIR', os.path.join(HERE, '..', 'i2'))
os.environ['MSS_REPO'] = REPO
argv = sys.argv; sys.argv = [argv[0]]; sys.path.insert(0, T10)
import t10_pairct as T
sys.argv = argv
T.PARI.default('nbthreads', 1)
LOG = []
T.say = lambda s='': LOG.append(s)          # t10's detailed lines go to the .log file, not the diffed output
PDS = []
_PairData = T.PairData
class _PD(_PairData):
    def __init__(s, *a, **k):
        _PairData.__init__(s, *a, **k); PDS.append(s)
T.PairData = _PD

def read_list(path):
    out = []
    for l in open(path):
        l = l.strip()
        if l and l[0] != '#': m, n = map(int, l.split()[0].split(',')); out.append((m, n))
    return out

def residual_records():
    recs = {}
    for f in ('i2_residual.jsonl', 'i2_residual_0of1.jsonl', 'i2_residual_0of1.jsonl.part'):
        p = os.path.join(I2, f)
        if os.path.exists(p):
            for l in open(p):
                try: o = json.loads(l)
                except ValueError: continue
                if 'residual' in o: recs[(o['m'], o['n'])] = o
    return recs

def check_record(o):
    """each tuple allowed at every bad prime and not at its killing prime (the record's completeness is I2's own check)."""
    m, n = o['m'], o['n']
    forms, sup, known = T.TD.magic(m, n)
    bad = T.badprimes(m, n); cache = {}
    def al(p):
        if p not in cache: cache[p] = T.X.allowed_at(forms, p)
        return cache[p]
    for t, k in zip(o['residual'], o['kill']):
        for p in bad:
            L, a = al(p)
            if tuple(L.dlabel(d) for d in t) not in a: return False, 'tuple %s not allowed at the bad prime %d' % (t, p)
        if k is not None:
            if k in bad: return False, 'killing prime %d is bad' % k
            L, a = al(k)
            if tuple(L.dlabel(d) for d in t) in a: return False, 'tuple %s allowed at its killing prime %d' % (t, k)
    return True, ''

def npoints(A, p):
    a2, a4, a6 = [x % p for x in A]; c = 1
    for x in range(p):
        f = (x * x * x + a2 * x * x + a4 * x + a6) % p
        c += 1 if f == 0 else (2 if pow(f, (p - 1) // 2, p) == 1 else 0)
    return c

def rank0(Pd, drop=0):
    """for a pair curve with dim Sel^2 = 2: the complete torsion group and the t of its points on C_jk."""
    if len(Pd.sel) != 2: return None
    Pr = Pd.Pr; es = Pr.es
    diffs = [es[0] - es[1], es[0] - es[2], es[1] - es[2]]
    B = 0; ps = []; p = 3
    while len(ps) < 25:
        if T.PARI.isprime(p) and all(d % p for d in diffs): ps.append(p); B = gcd(B, npoints(Pr.A, p))
        p += 2
    E = T.PARI.ellinit([0, Pr.A[0], 0, Pr.A[1], Pr.A[2]])
    gens = [(Fr(str(P[0])), Fr(str(P[1]))) for P in T.PARI.elltors(E)[2] if len(P) == 2]   # a finder; checked
    gens += [(Fr(e), Fr(0)) for e in es]
    for P in gens: assert Pr.on(P), 'a torsion point is not on the curve'
    grp = [None]
    while True:
        new = [Pr.add(P, g) for P in grp for g in gens]
        new = [Q for Q in new if Q not in grp]
        if not new: break
        for Q in new:
            if Q not in grp: grp.append(Q)
        assert len(grp) <= B, 'more rational torsion points than #E(F_p) allows'
    grp = grp[:len(grp) - drop]
    ts = []
    for P in grp:
        uv = Pr.point_uv(P) if P is not None else None
        if uv is None: t0 = Pr.t0; uv = (1, 0) if t0 == 'inf' else (Fr(t0).numerator, Fr(t0).denominator)
        ts.append('%d/%d' % uv if uv[1] else 'inf')
    complete = len(grp) == B
    triv = all((lambda u, v: u * v * (u * u - v * v) == 0)(*([int(x) for x in t.split('/')] if t != 'inf' else (1, 0))) for t in ts)
    return {'torsion': len(grp), 'bound': B, 'complete': complete, 't': sorted(set(ts)), 'excludes': complete and triv}

def analyse(m, n, rec=None):
    del PDS[:]
    if rec is not None:
        forms, sup, known = T.TD.magic(m, n)
        left = [tuple(t) for t in rec['residual']]; kill = dict(zip(left, rec['kill']))
        _, _, res, _ = T.analyse_ratio(m, n, tuples=(left, list(known)))
        src = 'i2_residual'
    else:
        left, kill, res, known = T.analyse_ratio(m, n); src = 't10'
    r0 = {}
    for Pd in PDS:
        x = rank0(Pd)
        if x: r0[Pd.j + Pd.k] = x
    excl = sorted(jk for jk, x in r0.items() if x['excludes'])
    tups = []
    for t in left:
        if t in known: continue
        ks = [jk for jk, v in res[t].items() if v in ('killed', 'outside Sel2')]
        if ks: ver = 'pairing'
        elif excl: ver = 'rank 0 of E_' + ','.join(excl)
        elif any(v == 'undecided' for v in res[t].values()): ver = 'undecided'
        else: ver = 'neither'
        tups.append({'tau': list(t), 'kill': kill.get(t), 'pairs': res[t], 'verdict': ver, 'killed_by': ks})
    return {'m': m, 'n': n, 'source': src, 'tuples': tups, 'rank0': r0}

def line(o): return json.dumps(o, separators=(', ', ': '), sort_keys=True)

def summary(rows):
    out = ['I8 summary: the Cassels-Tate pairing test of conductor-t10 on the ratios that need a good prime. Written by an AI model; not read by a person. Experiment, not counted; tier D.']
    c = {}; per = {}; nr = 0; err = []; allk = 0; neither = []; src = {}
    for o in rows:
        if 'tuples' not in o: err.append((o['m'], o['n'], o.get('status'))); continue
        nr += 1; src[o['source']] = src.get(o['source'], 0) + 1
        allk += all(t['verdict'] == 'pairing' for t in o['tuples'])
        for t in o['tuples']:
            v = t['verdict'] if not t['verdict'].startswith('rank 0') else 'rank 0'
            c[v] = c.get(v, 0) + 1
            for jk in t['killed_by']: per[jk] = per.get(jk, 0) + 1
            if v in ('neither', 'undecided'): neither.append((o['m'], o['n'], t['tau'], t['kill'], v))
    out.append('ratios done: %d (residual tuples from %s); not done: %d' % (nr, ', '.join('%s %d' % kv for kv in sorted(src.items())), len(err)))
    out.append('residual tuples: %d; pairing-killed %d; ratio excluded by a rank-0 pair curve %d; neither %d; undecided %d'
               % (sum(c.values()), c.get('pairing', 0), c.get('rank 0', 0), c.get('neither', 0), c.get('undecided', 0)))
    out.append('ratios all of whose tuples are pairing-killed: %d of %d' % (allk, nr))
    out.append('kills per pair curve: %s' % ' '.join('%s:%d' % (jk, per.get(jk, 0)) for jk in ('08', '0L', '0U', '8L', '8U', 'LU')))
    for m, n, t, k, v in neither: out.append('  %s: (%d,%d) %s killing good prime %s' % (v.upper(), m, n, t, k))
    for o in rows:
        for t in o.get('tuples', []):
            if t['verdict'].startswith('rank 0'): out.append('  rank 0: (%d,%d) %s by %s' % (o['m'], o['n'], t['tau'], t['verdict']))
    for e in err: out.append('  not done: (%d,%d) %s' % e)
    return out

def compare50():
    recs = residual_records(); keys = read_list(os.path.join(I2, 'i2_needgood_list.txt'))[:2900]   # fixed: done by I2 already
    assert all(k in recs for k in keys)
    out = ['I8: the residual tuples of i2_residual against t10_pairct.residual() on 50 of the first 2,900 ratios of i2_needgood_list.txt, drawn with seed 1. Written by an AI model; not read by a person.']
    agree = 0
    for m, n in sorted(random.Random(1).sample(keys, 50)):
        left, kill, bad, known = T.residual(m, n)
        o = recs[(m, n)]; mine = {tuple(t): k for t, k in zip(o['residual'], o['kill'])}
        theirs = {t: kill.get(t) for t in left}
        ok = mine == theirs; agree += ok
        out.append('(%d,%d): %d tuples; %s' % (m, n, len(mine), 'same tuples and killing primes' if ok else 'DIFFER'))
    out.append('agree: %d of 50' % agree)
    return out

def negative():
    out = ['I8 negative and positive controls. Written by an AI model; not read by a person.']
    recs = residual_records()
    o = recs[(15001, 2081)]
    ok, _ = check_record(o); out.append('true residual record of (15001,2081): %s' % ('accepted' if ok else 'REFUSED')); assert ok
    f = json.loads(json.dumps(o)); f['residual'][0][3] *= 3
    ok, why = check_record(f); out.append('forged residual tuple: %s (%s)' % ('refused' if not ok else 'ACCEPTED', why)); assert not ok
    bad = T.badprimes(15001, 2081); q = max(int(p) for p in T.PARI.primes([3, o['kill'][0] - 1]) if int(p) not in bad)
    f = json.loads(json.dumps(o)); f['kill'][0] = q
    ok, why = check_record(f); out.append('forged killing prime %d (a good prime below the true one): %s (%s)' % (q, 'refused' if not ok else 'ACCEPTED', why)); assert not ok
    r = analyse(5623, 1017)
    x = r['rank0'].get('08')
    out.append('positive control (5623,1017) (t10: E_08 of rank 0): E_08 torsion %d of bound %d, points over t in %s: %s'
               % (x['torsion'], x['bound'], x['t'], 'excludes the ratio' if x['excludes'] else 'DOES NOT EXCLUDE'))
    assert x['excludes'] and x['torsion'] == 8
    vs = sorted(set(t['verdict'] for t in r['tuples']))
    out.append('  its tuples: %s' % vs)
    Pd = [p for p in PDS if p.j + p.k == '08'][0]
    y = rank0(Pd, drop=1)
    out.append('forged: the torsion of E_08 with one point left out: %s' % ('refused (not complete)' if not y['excludes'] else 'ACCEPTED')); assert not y['excludes']
    Pd = [p for p in PDS if p.j + p.k == '0L'][0]
    out.append('forged: a rank-0 claim on E_0L (dim Sel2 = %d): %s' % (len(Pd.sel), 'refused' if rank0(Pd) is None else 'ACCEPTED')); assert rank0(Pd) is None
    r = analyse(4457, 287)
    out.append('positive control (4457,287) (t10: both tuples pairing-killed by 0L,0U,8L,8U): %s'
               % '; '.join('%s %s by %s' % (t['tau'], t['verdict'], ','.join(t['killed_by'])) for t in r['tuples']))
    assert all(t['verdict'] == 'pairing' and t['killed_by'] == ['0L', '0U', '8L', '8U'] for t in r['tuples'])
    return out

if __name__ == '__main__':
    W = '--write' in argv; a = [x for x in argv[1:] if not x.startswith('--')]
    def arg(k, d=None): return argv[argv.index(k) + 1] if k in argv else d
    if '--negative' in argv or '--compare50' in argv:
        out = negative() if '--negative' in argv else compare50()
        f = os.path.join(HERE, 'i8_negative.txt' if '--negative' in argv else 'i8_compare50.txt'); s = '\n'.join(out) + '\n'
        if W: open(f, 'w').write(s)
        else:
            sys.stdout.write(s)
            if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS')
    elif '--merge' in argv:
        N = int(a[0]); files = [[l for l in open(os.path.join(HERE, 'i8_%dof%d.jsonl' % (i, N)))] for i in range(N)]
        rows = []; j = 0
        while any(files): rows.append(files[j % N].pop(0)); j += 1
        s = ''.join(rows); s2 = '\n'.join(summary([json.loads(l) for l in rows])) + '\n'
        if W: open(os.path.join(HERE, 'i8_pairct.jsonl'), 'w').write(s); open(os.path.join(HERE, 'i8_summary.txt'), 'w').write(s2)
        else: sys.stdout.write(s2)
    else:
        i, N = int(a[0]), int(a[1]); K = int(arg('--max', '1000000'))
        lst = [r for j, r in enumerate(read_list(os.path.join(I2, 'i2_needgood_list.txt'))) if j % N == i]
        f = os.path.join(HERE, 'i8_%dof%d.jsonl' % (i, N)); part = f + '.part'
        done = [l.rstrip('\n') for l in open(part)] if W and os.path.exists(part) else []
        pf = open(part, 'a') if W else None; logf = open(f[:-6] + '.log', 'a') if W else None
        recs = residual_records(); out = list(done); new = 0
        for m, n in lst[len(done):]:
            if new >= K: sys.exit(3)
            rec = recs.get((m, n))
            try:
                if rec is not None:
                    ok, why = check_record(rec); assert ok, why
                s = line(analyse(m, n, rec))
            except MemoryError: s = line({'m': m, 'n': n, 'status': 'over the per-process memory limit, not done'})
            except AssertionError as ex: s = line({'m': m, 'n': n, 'status': 'assertion: %s' % str(ex)[:160]})
            except Exception as ex: s = line({'m': m, 'n': n, 'status': 'error: %s' % type(ex).__name__})
            out.append(s); new += 1
            if pf: pf.write(s + '\n'); pf.flush()
            if logf: logf.write('\n'.join(LOG) + '\n'); logf.flush(); del LOG[:]
        s = '\n'.join(out) + '\n'
        if W:
            open(f, 'w').write(s)
            if os.path.exists(part): os.remove(part)
        else: sys.stdout.write(s)
