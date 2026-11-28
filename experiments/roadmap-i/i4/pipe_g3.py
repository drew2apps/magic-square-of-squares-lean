"""G3, the stages after the second descent for one chunk (lo < m <= lo + 1000): the rank-0 quotients (the repository's
global.py) and Demjanenko-Manin (../../census2000/dm_pari.py), each on the open ratios together with their partners
(a ratio closes if it or its partner closes, as in CENSUS9000.md).  Same programs as pipe11k.py.  EXPERIMENT, written
by an AI model; not read by a person.   usage: pipe_g3.py LO"""
import json, sys, os, subprocess, math
HERE = os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..', '..'))
G = os.path.join(REPO, 'scripts', 'g1_full', 'global.py')
def partner(m, n):
    a, b = m + n, abs(m - n)
    if a % 2 == 0 and b % 2 == 0: a, b = a // 2, b // 2
    g = math.gcd(a, b); a, b = a // g, b // g
    return [max(a, b), min(a, b)]
lo = int(sys.argv[1]); tag = str(lo)
op = []; allp = []
for k in range(2):
    for l in open('sd_%s_%d.jsonl' % (tag, k)):
        if l.startswith('{'):
            o = json.loads(l); allp.append(o['pair'])
            if o['status'] == 'open': op.append(o['pair'])
pm = {str(tuple(p)): partner(*p) for p in op}
U = sorted(set(tuple(p) for p in op) | set(tuple(v) for v in pm.values()))
procs = []
for k in range(2):
    json.dump({'open': [list(u) for u in U[k::2]]}, open('r0_%s_in%d.json' % (tag, k), 'w'))
    procs.append(subprocess.Popen('nice python3 %s r0_%s_in%d.json r0_%s_out%d.json > r0_%s_%d.txt 2>&1' % (G, tag, k, tag, k, tag, k), shell=True))
for p in procs: p.wait()
closed = set()
for k in range(2):
    for key, lst in json.load(open('r0_%s_out%d.json' % (tag, k))).items():
        if any(c.get('closes') for c in lst): closed.add(tuple(int(t) for t in key.split(',')))
left = [p for p in op if tuple(p) not in closed and tuple(pm[str(tuple(p))]) not in closed]
U = sorted(set(tuple(p) for p in left) | set(tuple(pm[str(tuple(p))]) for p in left))
procs = []
for k in range(2):
    json.dump({'open': [list(u) for u in U[k::2]]}, open('dm_%s_in%d.json' % (tag, k), 'w'))
    procs.append(subprocess.Popen('nice python3 ../../census2000/dm_pari.py dm_%s_in%d.json dm_%s_%d.jsonl > dm_%s_%d.log 2>&1' % (tag, k, tag, k, tag, k), shell=True))
for p in procs: p.wait()
dmc = set(); alerts = []
for k in range(2):
    for l in open('dm_%s_%d.jsonl' % (tag, k)):
        if l.startswith('{'):
            o = json.loads(l)
            if o.get('status') == 'closed': dmc.add(tuple(o['pair']))
            if 'alert' in str(o.get('status')): alerts.append(o['pair'])
left2 = [p for p in left if tuple(p) not in dmc and tuple(pm[str(tuple(p))]) not in dmc]
res = dict(chunk=[lo, lo + 1000], ratios=len(allp), after_second_descent=len(op), after_rank0=len(left), after_dm=len(left2),
           open=left2, dm_alerts=alerts)
json.dump(res, open('stage_%s.json' % tag, 'w')); print(json.dumps(res))
