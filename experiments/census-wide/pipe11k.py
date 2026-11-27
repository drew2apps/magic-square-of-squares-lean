"""Stages after the second descent for 9000 < m <= 11000 (same programs as the 7000 run).  EXPERIMENT."""
import json, sys, subprocess, math
G = __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..', '..', 'scripts', 'g1_full', 'global.py')
def partner(m, n):
    a, b = m + n, abs(m - n)
    if a % 2 == 0 and b % 2 == 0: a, b = a // 2, b // 2
    g = math.gcd(a, b); a, b = a // g, b // g
    return [max(a, b), min(a, b)]
stage = sys.argv[1]
if stage == 'rank0':
    op = []
    for f in ['s11_0.jsonl','s11_1.jsonl']:
        for l in open(f):
            if l.startswith('{'):
                o = json.loads(l)
                if o['status'] == 'open': op.append(o['pair'])
    pm = {str(tuple(p)): partner(*p) for p in op}
    json.dump({'pairs': op, 'partner': pm}, open('open11k_map.json', 'w'))
    U = sorted(set(tuple(p) for p in op) | set(tuple(v) for v in pm.values()))
    print(len(op), 'open,', len(U), 'with partners')
    procs = []
    for k in range(2):
        json.dump({'open': [list(u) for u in U[k::2]]}, open(f'g11k_in{k}.json', 'w'))
        procs.append(subprocess.Popen(f'python3 {G} g11k_in{k}.json g11k_out{k}.json > g11k{k}.txt 2>&1', shell=True))
    for p in procs: p.wait()
    closed = set()
    for k in range(2):
        for key, lst in json.load(open(f'g11k_out{k}.json')).items():
            if any(c.get('closes') for c in lst): closed.add(tuple(int(t) for t in key.split(',')))
    left = [p for p in op if tuple(p) not in closed and tuple(pm[str(tuple(p))]) not in closed]
    json.dump(left, open('after_rank0_11k.json', 'w')); print('after rank 0:', len(left))
if stage == 'dm':
    left = json.load(open('after_rank0_11k.json')); pm = json.load(open('open11k_map.json'))['partner']
    U = sorted(set(tuple(p) for p in left) | set(tuple(pm[str(tuple(p))]) for p in left))
    procs = []
    for k in range(2):
        json.dump({'open': [list(u) for u in U[k::2]]}, open(f'dm11k_in{k}.json', 'w'))
        procs.append(subprocess.Popen(f'python3 ../census2000/dm_pari.py dm11k_in{k}.json dm11k{k}.jsonl > dm11k{k}.log 2>&1', shell=True))
    for p in procs: p.wait()
    closed = set(); alerts = []
    for k in range(2):
        for l in open(f'dm11k{k}.jsonl'):
            if l.startswith('{'):
                o = json.loads(l)
                if o.get('status') == 'closed': closed.add(tuple(o['pair']))
                if 'alert' in str(o.get('status')): alerts.append(o)
    left2 = [p for p in left if tuple(p) not in closed and tuple(pm[str(tuple(p))]) not in closed]
    json.dump(left2, open('after_dm_11k.json', 'w')); print('after DM:', len(left2), left2, 'alerts', len(alerts))
